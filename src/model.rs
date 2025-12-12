use crate::analysis::calculate_ressens;
use crate::parameters::Parameters;
use anyhow::Result;
use std::f64::consts::PI;
use std::fs::File;
use std::io::{BufWriter, Write};

pub struct Model {
    pub params: Parameters,
    tapetal_pigment: f64,
    shielding_pigment: f64,
    pub ommatidial_angle: f64, // Public for analysis
    number_of_facets: usize,
    rhabdom_radius: f64,
    incidence_ommatidial_angle: f64,
    critical_angle: f64,
    refracted_ommatidial_angle: f64,
    old_rhabdom_length: f64,
}

impl Model {
    pub fn new(mut params: Parameters) -> Self {
        // Run post-processing on parameters (e.g., clamp blur circle)
        params.post_process();

        let old_rhabdom_length = params.rhabdom_length;
        let circumference_of_eye = PI * params.eye_diameter;
        let aperture_radius = params.aperture_diameter / 2.0;
        let eye_radius = params.eye_diameter / 2.0;

        let distance_to_aperture = (eye_radius.powi(2) - aperture_radius.powi(2)).sqrt();
        let angle_at_center = (aperture_radius / distance_to_aperture).atan().to_degrees();
        let aperture_arc = (angle_at_center / 360.0) * circumference_of_eye;
        let ommatidial_angle = (params.facet_width / circumference_of_eye) * 360.0;
        let rhabdom_radius = params.rhabdom_width / 2.0;

        let number_of_facets = (aperture_arc / params.facet_width).round() as usize;

        let snells_law =
            (params.cytoplasm_refractive_index / params.rhabdom_refractive_index).asin();
        let critical_angle = 90.0 - snells_law.to_degrees();

        Model {
            params,
            tapetal_pigment: 0.0,
            shielding_pigment: 0.0,
            ommatidial_angle,
            number_of_facets,
            rhabdom_radius,
            incidence_ommatidial_angle: 0.0,
            critical_angle,
            refracted_ommatidial_angle: 0.0,
            old_rhabdom_length,
        }
    }

    pub fn run_simulation(&mut self) -> Result<()> {
        let p = &self.params;
        println!("--- Running simulation for {} ---", p.species_name);
        println!("Calculating pathlengths for {}...", p.species_name);

        let pathlengths_filename = format!("{}_pathlengths.csv", p.species_name);
        let debug_filename = format!("{}_debug.csv", p.species_name);

        let pathlengths_file = File::create(&pathlengths_filename)?;
        let mut pathlengths_writer = BufWriter::new(pathlengths_file);

        let debug_file = File::create(&debug_filename)?;
        let mut debug_writer = BufWriter::new(debug_file);

        self.shielding_pigment = 0.0;
        self.tapetal_pigment = 0.0;
        let increment_amount = p.rhabdom_length / 10.0;

        loop {
            println!(
                "P: {:.2}, T: {:.2}",
                self.shielding_pigment, self.tapetal_pigment
            );
            writeln!(
                pathlengths_writer,
                "{:.6}\n{:.6}",
                self.shielding_pigment, self.tapetal_pigment
            )?;
            writeln!(
                debug_writer,
                "P: {:.2}, T: {:.2}",
                self.shielding_pigment, self.tapetal_pigment
            )?;

            self.incidence_ommatidial_angle = 0.0;
            self.refracted_ommatidial_angle = 0.0;
            let mut current_facet = 0;

            while current_facet < self.number_of_facets {
                let mut row_data: Vec<String> = Vec::new();

                // Account for refraction at the cornea
                if self.incidence_ommatidial_angle > 0.0 && self.incidence_ommatidial_angle <= 15.0
                {
                    self.refracted_ommatidial_angle =
                        (self.incidence_ommatidial_angle * 0.9494) + 0.004667;
                } else if self.incidence_ommatidial_angle > 15.0
                    && self.incidence_ommatidial_angle <= 35.0
                {
                    self.refracted_ommatidial_angle =
                        (self.incidence_ommatidial_angle * 0.9407) + 0.1648;
                } else if self.incidence_ommatidial_angle > 35.0
                    && self.incidence_ommatidial_angle <= 50.0
                {
                    self.refracted_ommatidial_angle =
                        (self.incidence_ommatidial_angle * 0.9196) + 0.8676;
                } else if self.incidence_ommatidial_angle > 50.0
                    && self.incidence_ommatidial_angle <= 60.0
                {
                    self.refracted_ommatidial_angle =
                        (self.incidence_ommatidial_angle * 0.8677) + 3.38;
                } else if self.incidence_ommatidial_angle > 60.0 {
                    println!("UNREAL ANGLE AT CORNEA");
                    row_data.push("UNREAL ANGLE AT CORNEA".to_string());
                }

                // Light loss at cone due to angle of incidence
                let mut facet_num = if self.refracted_ommatidial_angle == 0.0 {
                    1.0
                } else {
                    let cc =
                        p.facet_width / self.refracted_ommatidial_angle.to_radians().tan().abs();
                    let fw = if cc > p.facet_width * 2.0 {
                        self.incidence_ommatidial_angle.to_radians().cos() * p.facet_width
                    } else {
                        let ll = (2.0 * cc) - (2.0 * p.facet_width);
                        self.incidence_ommatidial_angle.to_radians().sin() * ll
                    };
                    fw / p.facet_width
                };
                if facet_num > 1.0 {
                    facet_num = 1.0;
                }

                // --- Path Calculation Logic ---
                let rhabdom_length = self.old_rhabdom_length;

                if self.incidence_ommatidial_angle == 0.0 {
                    // CASE 4: Perpendicular ray
                    let val = if self.tapetal_pigment == 0.0 || self.shielding_pigment > 0.0 {
                        rhabdom_length * facet_num
                    } else {
                        (rhabdom_length * 2.0) * facet_num
                    };
                    row_data.push(format!("{:.6}", val));
                } else {
                    let y = self.rhabdom_radius
                        / self.refracted_ommatidial_angle.to_radians().tan().abs();

                    if y >= rhabdom_length {
                        // CASE 3: Bounce off base
                        let mx = (rhabdom_length.powi(2) + self.rhabdom_radius.powi(2)).sqrt();
                        let x = if y == rhabdom_length {
                            mx
                        } else {
                            rhabdom_length
                                / self.refracted_ommatidial_angle.to_radians().cos().abs()
                        };

                        let v = if x > self.old_rhabdom_length {
                            x
                        } else {
                            self.old_rhabdom_length
                        };

                        let val = if self.tapetal_pigment == 0.0 || self.shielding_pigment > 0.0 {
                            x * facet_num
                        } else {
                            (x + v) * facet_num
                        };
                        row_data.push(format!("{:.6}", val));
                    } else if y > (rhabdom_length - self.shielding_pigment)
                        || y > (rhabdom_length - self.tapetal_pigment)
                        || self.refracted_ommatidial_angle < self.critical_angle
                    {
                        // CASE 2: Reflection from edge
                        let x = self.rhabdom_radius
                            / self.refracted_ommatidial_angle.to_radians().sin().abs();
                        let mut z = (rhabdom_length - y)
                            / self.refracted_ommatidial_angle.to_radians().cos().abs();
                        if z > x {
                            z = x;
                        }

                        let v = if (x + z) > self.old_rhabdom_length {
                            x + z
                        } else {
                            self.old_rhabdom_length
                        };

                        let mut val = if self.tapetal_pigment == 0.0 {
                            (x + z) * facet_num
                        } else {
                            (x + z + v) * facet_num
                        };

                        if self.shielding_pigment > 0.0 {
                            val = (x + z) * facet_num;
                        }
                        if self.shielding_pigment > (rhabdom_length - y) {
                            val = x * facet_num;
                        }
                        row_data.push(format!("{:.6}", val));
                    } else {
                        // CASE 1: No reflection
                        let boa = self.refracted_ommatidial_angle;
                        let x = self.rhabdom_radius / boa.to_radians().sin().abs();
                        row_data.push(format!("{:.6}", x * facet_num));
                    }
                }

                current_facet += 1;
                self.incidence_ommatidial_angle += self.ommatidial_angle;

                // Account for blur circle
                if p.blur_circle_extent > 0.0 {
                    let fd = (self.number_of_facets as f64) / p.blur_circle_extent;
                    let mut nx = 0;
                    for _ in 0..p.blur_circle_extent as i32 {
                        nx += 1;
                        if (current_facet as f64) > (fd * nx as f64) {
                            self.refracted_ommatidial_angle += self.ommatidial_angle;
                            row_data.push("0".to_string());
                        }
                    }
                }

                row_data.push("998".to_string());
                writeln!(pathlengths_writer, "{}", row_data.join(","))?;
            }

            writeln!(pathlengths_writer, "999")?;

            // Pigment increment logic
            if self.tapetal_pigment >= self.old_rhabdom_length
                && self.shielding_pigment >= self.old_rhabdom_length
            {
                break;
            } else if self.tapetal_pigment >= self.old_rhabdom_length {
                self.tapetal_pigment = 0.0;
                self.shielding_pigment += increment_amount;
            } else {
                self.tapetal_pigment += increment_amount;
            }
        }

        // Calculate results
        calculate_ressens(&p.species_name, self.ommatidial_angle)?;

        println!("--- Finished simulation for {} ---\n", p.species_name);

        Ok(())
    }
}
