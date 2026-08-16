use crate::analysis::calculate_ressens;
use crate::parameters::Parameters;
use anyhow::Result;
use std::f64::consts::PI;
use std::fs::File;
use std::io::{BufWriter, Write};

pub struct Model {
    pub params: Parameters,
    pub circumference_of_eye: f64,
    pub aperture_arc: f64,
    pub ommatidial_angle: f64,
    pub number_of_facets: usize,
    pub rhabdom_radius: f64,
    pub critical_angle: f64,
    pub old_rhabdom_length: f64,
    pub debug_mode: bool,
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
            circumference_of_eye,
            aperture_arc,
            ommatidial_angle,
            number_of_facets,
            rhabdom_radius,
            critical_angle,
            old_rhabdom_length,
            debug_mode: false,
        }
    }

    pub fn run_simulation(&mut self) -> Result<()> {
        let p = &self.params;
        println!("--- Running simulation for {} ---", p.species_name);
        println!("Calculating pathlengths for {}...", p.species_name);

        let pathlengths_filename = format!("{}_pathlengths.csv", p.species_name);
        let pathlengths_file = File::create(&pathlengths_filename)?;
        let mut pathlengths_writer = BufWriter::new(pathlengths_file);

        let mut debug_writer = if self.debug_mode {
            let debug_filename = format!("{}_debug.csv", p.species_name);
            let debug_file = File::create(&debug_filename)?;
            Some(BufWriter::new(debug_file))
        } else {
            None
        };

        let increment_amount = p.rhabdom_length / 10.0;

        for s_step in 0..=10 {
            let shielding_pigment = s_step as f64 * increment_amount;
            for t_step in 0..=10 {
                let tapetal_pigment = t_step as f64 * increment_amount;

                println!("P: {:.2}, T: {:.2}", shielding_pigment, tapetal_pigment);
                writeln!(
                    pathlengths_writer,
                    "{:.6}\n{:.6}",
                    shielding_pigment, tapetal_pigment
                )?;
                if let Some(ref mut dw) = debug_writer {
                    writeln!(dw, "P: {:.2}, T: {:.2}", shielding_pigment, tapetal_pigment)?;
                }

                for current_facet in 0..self.number_of_facets {
                    let incidence_ommatidial_angle = current_facet as f64 * self.ommatidial_angle;

                    // Account for refraction at the cornea
                    let mut refracted_ommatidial_angle = 0.0;
                    if incidence_ommatidial_angle > 0.0 && incidence_ommatidial_angle <= 15.0 {
                        refracted_ommatidial_angle =
                            (incidence_ommatidial_angle * 0.9494) + 0.004667;
                    } else if incidence_ommatidial_angle > 15.0
                        && incidence_ommatidial_angle <= 35.0
                    {
                        refracted_ommatidial_angle = (incidence_ommatidial_angle * 0.9407) + 0.1648;
                    } else if incidence_ommatidial_angle > 35.0
                        && incidence_ommatidial_angle <= 50.0
                    {
                        refracted_ommatidial_angle = (incidence_ommatidial_angle * 0.9196) + 0.8676;
                    } else if incidence_ommatidial_angle > 50.0
                        && incidence_ommatidial_angle <= 60.0
                    {
                        refracted_ommatidial_angle = (incidence_ommatidial_angle * 0.8677) + 3.38;
                    } else if incidence_ommatidial_angle > 60.0 {
                        println!("UNREAL ANGLE AT CORNEA");
                        refracted_ommatidial_angle = 60.0;
                    }

                    // Light loss at cone due to angle of incidence
                    let facet_num = if refracted_ommatidial_angle == 0.0 {
                        1.0
                    } else {
                        let cc =
                            p.facet_width / refracted_ommatidial_angle.to_radians().tan().abs();
                        let fw = if cc > p.facet_width * 2.0 {
                            incidence_ommatidial_angle.to_radians().cos() * p.facet_width
                        } else {
                            let ll = (2.0 * cc) - (2.0 * p.facet_width);
                            incidence_ommatidial_angle.to_radians().sin() * ll
                        };
                        fw / p.facet_width
                    };
                    let facet_num = facet_num.clamp(0.0, 1.0);

                    let mut row_data: Vec<String> = Vec::new();
                    let mut boa = refracted_ommatidial_angle;

                    if p.blur_circle_extent > 0.0 {
                        let fd = (self.number_of_facets as f64) / p.blur_circle_extent;
                        for i in 1..=p.blur_circle_extent as usize {
                            if (current_facet as f64) > (fd * i as f64) {
                                boa += self.ommatidial_angle;
                                row_data.push("0".to_string());
                            }
                        }
                    }

                    // Ray tracing through rhabdom array
                    let mut rhabdom_length = self.old_rhabdom_length;
                    let mut cz = 0;

                    loop {
                        // Account for tapered/pointy rhabdom shape at proximal entrance
                        if boa > self.critical_angle && cz == 0 {
                            boa -= p.proximal_rhabdom_angle;
                        }

                        if incidence_ommatidial_angle == 0.0 {
                            // CASE 4: Perpendicular ray
                            let val = if tapetal_pigment == 0.0 || shielding_pigment > 0.0 {
                                rhabdom_length * facet_num
                            } else {
                                (rhabdom_length * 2.0) * facet_num
                            };
                            row_data.push(format!("{:.6}", val));
                            break;
                        }

                        let y = self.rhabdom_radius / boa.to_radians().tan().abs();

                        if y >= rhabdom_length {
                            // CASE 3: Bounce off base
                            let mx = (rhabdom_length.powi(2) + self.rhabdom_radius.powi(2)).sqrt();
                            let x = if y == rhabdom_length {
                                mx
                            } else {
                                rhabdom_length / boa.to_radians().cos().abs()
                            };

                            let v = self.old_rhabdom_length / boa.to_radians().cos().abs();

                            let val = if tapetal_pigment == 0.0 || shielding_pigment > 0.0 {
                                x * facet_num
                            } else {
                                (x + v) * facet_num
                            };
                            row_data.push(format!("{:.6}", val));
                            break;
                        } else if y > (rhabdom_length - shielding_pigment)
                            || y > (rhabdom_length - tapetal_pigment)
                            || boa < self.critical_angle
                        {
                            // CASE 2: Reflection from edge
                            let x = self.rhabdom_radius / boa.to_radians().sin().abs();
                            let z = (rhabdom_length - y) / boa.to_radians().cos().abs();

                            let v = self.old_rhabdom_length / boa.to_radians().cos().abs();

                            let mut val = if tapetal_pigment == 0.0 {
                                (x + z) * facet_num
                            } else {
                                (x + z + v) * facet_num
                            };

                            if shielding_pigment > 0.0 {
                                val = (x + z) * facet_num;
                            }
                            if shielding_pigment > (rhabdom_length - y) {
                                val = x * facet_num;
                            }
                            row_data.push(format!("{:.6}", val));
                            break;
                        } else {
                            // CASE 1: No reflection (ray passes through rhabdom wall into adjacent rhabdom)
                            let x = self.rhabdom_radius / boa.to_radians().sin().abs();
                            row_data.push(format!("{:.6}", x * facet_num));
                            rhabdom_length -= y;
                            boa += self.ommatidial_angle;
                            cz = 1;
                            if rhabdom_length <= tapetal_pigment
                                || rhabdom_length <= shielding_pigment
                            {
                                break;
                            }
                        }
                    }

                    // Write row
                    writeln!(pathlengths_writer, "{}", row_data.join(","))?;
                }

                writeln!(pathlengths_writer, "999")?;
            }
        }

        drop(pathlengths_writer);
        drop(debug_writer);

        // Calculate results
        calculate_ressens(&p.species_name, self.ommatidial_angle)?;

        println!("--- Finished simulation for {} ---\n", p.species_name);

        Ok(())
    }
}
