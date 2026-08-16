use crate::analysis::{BlockSummary, accumulate, calculate_ressens, summarise_block};
use crate::parameters::Parameters;
use anyhow::{Result, bail};
use std::f64::consts::PI;
use std::fs::File;
use std::io::{BufWriter, Write};

/// Labels the columns of the raw geometry output. Every row carries its own keys, so
/// the file is a plain rectangular CSV with no positional state and no block terminator.
pub const PATHLENGTHS_HEADER: &str = "block,shielding_um,tapetal_um,facet,rhabdom,pathlength_um";

/// Number of migration positions sampled for each pigment, from fully retracted (0)
/// to fully covering the rhabdom (rhabdom_length).
pub const PIGMENT_STEPS: usize = 11;

/// The largest angle to the rhabdom axis at which a ray can still advance towards
/// the proximal end. At or beyond 90 degrees the ray travels perpendicular to (or
/// back along) the axis and is treated as lost.
pub const MAX_PROPAGATION_ANGLE: f64 = 90.0;

pub struct Model {
    pub params: Parameters,
    pub circumference_of_eye: f64,
    pub aperture_arc: f64,
    pub ommatidial_angle: f64,
    pub number_of_facets: usize,
    pub rhabdom_radius: f64,
    pub critical_angle: f64,
    pub debug_mode: bool,
}

/// The outcome of tracing one ray through the rhabdom array.
pub struct TraceResult {
    /// Pathlengths through each successive rhabdom the ray enters, in micrometres of
    /// raw geometry. Facet transmission is NOT folded in here; it is a flux factor
    /// applied when the absorbed intensity is computed.
    pub pathlengths: Vec<f64>,
    /// Which branch ended the trace, for the debug output.
    pub terminal_case: &'static str,
    /// The largest angle to the rhabdom axis reached during the trace. Every segment
    /// covers some axial depth d at an angle no greater than this, so the total path
    /// is bounded by 2*rhabdom_length/cos(max_angle).
    pub max_angle: f64,
    /// Set when the ray stopped propagating towards the proximal end.
    pub lost: bool,
}

/// Applies the empirical corneal refraction regression to an angle of incidence,
/// in degrees. Returns NaN for angles the regression does not cover.
pub fn refracted_angle(incidence: f64) -> f64 {
    if incidence <= 0.0 {
        0.0
    } else if incidence <= 15.0 {
        incidence * 0.9494 + 0.004667
    } else if incidence <= 35.0 {
        incidence * 0.9407 + 0.1648
    } else if incidence <= 50.0 {
        incidence * 0.9196 + 0.8676
    } else if incidence <= 60.0 {
        incidence * 0.8677 + 3.38
    } else {
        f64::NAN
    }
}

impl Model {
    pub fn new(params: Parameters) -> Result<Self> {
        params.validate()?;

        let circumference_of_eye = PI * params.eye_diameter;
        let aperture_radius = params.aperture_diameter / 2.0;
        let eye_radius = params.eye_diameter / 2.0;

        let distance_to_aperture = (eye_radius.powi(2) - aperture_radius.powi(2)).sqrt();
        let angle_at_center = (aperture_radius / distance_to_aperture).atan().to_degrees();
        let aperture_arc = (angle_at_center / 360.0) * circumference_of_eye;
        let ommatidial_angle = (params.facet_width / circumference_of_eye) * 360.0;
        let rhabdom_radius = params.rhabdom_width / 2.0;

        let number_of_facets = (aperture_arc / params.facet_width).round() as usize;

        // boa is measured from the rhabdom axis, so the angle at the wall normal is
        // (90 - boa) and light is guided while boa < critical_angle.
        let snells_law = (params.cytoplasm_refractive_index / params.rhabdom_refractive_index)
            .asin()
            .to_degrees();
        let critical_angle = 90.0 - snells_law;

        if number_of_facets < 1 {
            bail!(
                "eyeshine patch spans no facets (aperture arc {:.2} um / facet width {:.2} um); check aperture and facet dimensions",
                aperture_arc,
                params.facet_width
            );
        }
        // The blur circle spreads the light gathered across the eyeshine patch over
        // blur_circle_extent rhabdoms. With more blur steps than facets, the mapping
        // from facet to rhabdom offset leaves offsets with no contributing facet at
        // all, producing a comb-shaped profile whose half maximum is a binning
        // artefact rather than a measurement.
        if params.blur_circle_extent > number_of_facets as f64 {
            bail!(
                "blur circle extent ({} rhabdoms) exceeds the {} facets across the eyeshine patch; the resulting profile would contain rhabdom offsets that receive no light",
                params.blur_circle_extent,
                number_of_facets
            );
        }

        Ok(Model {
            params,
            circumference_of_eye,
            aperture_arc,
            ommatidial_angle,
            number_of_facets,
            rhabdom_radius,
            critical_angle,
            debug_mode: false,
        })
    }

    /// The fraction of light a facet admits at the given angle of incidence, relative
    /// to a facet viewed normally. This is a flux factor in [0, 1] and is applied to
    /// the absorbed intensity, not to the geometric path length.
    pub fn facet_transmission(&self, facet_index: usize) -> f64 {
        let incidence = facet_index as f64 * self.ommatidial_angle;
        let refracted = refracted_angle(incidence);
        if refracted.is_nan() {
            return 0.0;
        }
        if refracted == 0.0 {
            return 1.0;
        }
        let cc = self.params.facet_width / refracted.to_radians().tan().abs();
        let fw = if cc > self.params.facet_width * 2.0 {
            incidence.to_radians().cos() * self.params.facet_width
        } else {
            let ll = (2.0 * cc) - (2.0 * self.params.facet_width);
            incidence.to_radians().sin() * ll
        };
        (fw / self.params.facet_width).clamp(0.0, 1.0)
    }

    /// The displacement, in rhabdoms, of the image formed by the facet at the given
    /// radial position in the eyeshine patch.
    ///
    /// The blur circle spans blur_circle_extent rhabdoms, so the outermost facet is
    /// displaced by (blur_circle_extent - 1) and the central facet by zero. The
    /// offset is continuous in the facet index: quantising it to whole rhabdoms (as
    /// earlier versions did, via a chain of `facet > fd*i` tests) aliased the facets
    /// unevenly across the available offsets and cut notches into the profile at
    /// offsets where fd*i landed on an exact integer.
    pub fn blur_offset(&self, facet_index: usize) -> f64 {
        if self.number_of_facets <= 1 {
            return 0.0;
        }
        facet_index as f64 * (self.params.blur_circle_extent - 1.0)
            / (self.number_of_facets - 1) as f64
    }

    /// Follows a single ray from the given facet through the rhabdom array for one
    /// combination of pigment positions.
    pub fn trace_ray(&self, facet_index: usize, shielding: f64, tapetal: f64) -> TraceResult {
        let p = &self.params;
        let mut result = TraceResult {
            pathlengths: Vec::new(),
            terminal_case: "",
            max_angle: 0.0,
            lost: false,
        };

        // Angle to the rhabdom axis on entry: corneal refraction plus the blur-circle
        // displacement, which tilts the ray by one ommatidial angle per rhabdom offset.
        let mut boa = refracted_angle(facet_index as f64 * self.ommatidial_angle)
            + self.blur_offset(facet_index) * self.ommatidial_angle;

        let mut rhabdom_length = p.rhabdom_length;
        let mut cz = 0;

        loop {
            // Tapered ("pointy") rhabdom tip widens the acceptance angle on first entry.
            if boa > self.critical_angle && cz == 0 {
                boa -= p.proximal_rhabdom_angle;
                if boa < 0.0 {
                    boa = 0.0;
                }
            }

            // A ray at 90 degrees or more to the axis cannot advance towards the
            // proximal end. Earlier versions took the absolute value of tan and cos,
            // which silently folded such rays back and produced path lengths many
            // times the rhabdom length.
            if boa >= MAX_PROPAGATION_ANGLE || boa.is_nan() {
                result.terminal_case = "lost";
                result.lost = true;
                return result;
            }
            if boa > result.max_angle {
                result.max_angle = boa;
            }

            if facet_index == 0 {
                // CASE 4: axial ray. Equivalent to case 3 at boa = 0, kept explicit.
                let val = if tapetal > 0.0 && shielding == 0.0 {
                    rhabdom_length * 2.0
                } else {
                    rhabdom_length
                };
                result.pathlengths.push(val);
                result.terminal_case = "C4";
                return result;
            }

            let sin = boa.to_radians().sin();
            let cos = boa.to_radians().cos();
            let tan = boa.to_radians().tan();

            // Axial distance travelled before the ray meets the rhabdom wall.
            let y = self.rhabdom_radius / tan;

            if y >= rhabdom_length {
                // CASE 3: the ray reaches the base without meeting the wall and is
                // reflected by the tapetum at the base.
                let x = rhabdom_length / cos;
                let v = p.rhabdom_length / cos;
                let val = if tapetal > 0.0 && shielding == 0.0 {
                    x + v
                } else {
                    x
                };
                result.pathlengths.push(val);
                result.terminal_case = "C3";
                return result;
            } else if y > (rhabdom_length - shielding)
                || y > (rhabdom_length - tapetal)
                || boa < self.critical_angle
            {
                // CASE 2: the ray is reflected at the wall, either by total internal
                // reflection or by the tapetal mirror.
                let x = self.rhabdom_radius / sin;

                // A guided ray never leaves the rhabdom, so the proximal screening
                // pigment - which lies in the cytoplasm outside the rhabdom - cannot
                // absorb it. An unguided ray exits through the wall and is absorbed
                // where the pigment starts.
                let guided = boa < self.critical_angle;
                let axial = if shielding > 0.0 && !guided {
                    (rhabdom_length - shielding - y).max(0.0)
                } else {
                    rhabdom_length - y
                };
                let z = axial / cos;
                let v = p.rhabdom_length / cos;

                let val = if tapetal > 0.0 && shielding == 0.0 {
                    x + z + v
                } else {
                    x + z
                };
                result.pathlengths.push(val);
                result.terminal_case = "C2";
                return result;
            } else {
                // CASE 1: no reflection. The ray crosses the wall into the adjacent
                // rhabdom, and the inter-rhabdom angle steps by one ommatidial angle.
                result.pathlengths.push(self.rhabdom_radius / sin);
                rhabdom_length -= y;
                boa += self.ommatidial_angle;
                cz = 1;
                if rhabdom_length <= tapetal || rhabdom_length <= shielding {
                    result.terminal_case = "C1";
                    return result;
                }
            }
        }
    }

    pub fn run_simulation(&mut self) -> Result<()> {
        let p = &self.params;
        println!("--- Running simulation for {} ---", p.species_name);
        println!(
            "{} facets across the eyeshine patch, ommatidial angle {:.4} deg, critical angle {:.4} deg",
            self.number_of_facets, self.ommatidial_angle, self.critical_angle
        );
        println!("Calculating pathlengths for {}...", p.species_name);

        let pathlengths_filename = format!("{}_pathlengths.csv", p.species_name);
        let pathlengths_file = File::create(&pathlengths_filename)?;
        let mut pathlengths_writer = BufWriter::new(pathlengths_file);

        let mut debug_writer = if self.debug_mode {
            let debug_file = File::create(format!("{}_debug.csv", p.species_name))?;
            let mut w = BufWriter::new(debug_file);
            writeln!(
                w,
                "block,shielding_um,tapetal_um,facet,incidence_deg,refracted_deg,blur_offset_rhabdoms,entry_boa_deg,facet_transmission,terminal_case,rhabdoms_entered,pathlengths_um"
            )?;
            Some(w)
        } else {
            None
        };

        let increment_amount = p.rhabdom_length / 10.0;
        let mut lost_rays = 0usize;
        let mut block = 0usize;
        let mut summaries: Vec<BlockSummary> = Vec::with_capacity(PIGMENT_STEPS * PIGMENT_STEPS);

        writeln!(pathlengths_writer, "{}", PATHLENGTHS_HEADER)?;

        for s_step in 0..PIGMENT_STEPS {
            let shielding = s_step as f64 * increment_amount;
            for t_step in 0..PIGMENT_STEPS {
                let tapetal = t_step as f64 * increment_amount;

                // Area-weighted absorbed light at each rhabdom offset from the optic axis.
                let mut profile: Vec<f64> = Vec::new();

                for facet in 0..self.number_of_facets {
                    let trace = self.trace_ray(facet, shielding, tapetal);
                    if trace.lost {
                        lost_rays += 1;
                    }
                    accumulate(self, &mut profile, facet, &trace.pathlengths);

                    if trace.pathlengths.is_empty() {
                        // A lost ray absorbs nothing, but the facet still belongs in
                        // the record, so emit an explicit zero for it.
                        writeln!(
                            pathlengths_writer,
                            "{},{:.6},{:.6},{},0,0.000000",
                            block, shielding, tapetal, facet
                        )?;
                    }
                    for (rhabdom, v) in trace.pathlengths.iter().enumerate() {
                        writeln!(
                            pathlengths_writer,
                            "{},{:.6},{:.6},{},{},{:.6}",
                            block, shielding, tapetal, facet, rhabdom, v
                        )?;
                    }

                    if let Some(ref mut w) = debug_writer {
                        let parts: Vec<String> = trace
                            .pathlengths
                            .iter()
                            .map(|v| format!("{:.6}", v))
                            .collect();
                        let incidence = facet as f64 * self.ommatidial_angle;
                        writeln!(
                            w,
                            "{},{:.4},{:.4},{},{:.4},{:.4},{:.4},{:.4},{:.6},{},{},{}",
                            block,
                            shielding,
                            tapetal,
                            facet,
                            incidence,
                            refracted_angle(incidence),
                            self.blur_offset(facet),
                            refracted_angle(incidence)
                                + self.blur_offset(facet) * self.ommatidial_angle,
                            self.facet_transmission(facet),
                            trace.terminal_case,
                            trace.pathlengths.len(),
                            parts.join(" ")
                        )?;
                    }
                }

                summaries.push(summarise_block(self, &profile));
                block += 1;
            }
        }

        drop(pathlengths_writer);
        drop(debug_writer);

        if lost_rays > 0 {
            println!(
                "WARNING: {} of {} rays exceeded 90 degrees to the rhabdom axis and were discarded.",
                lost_rays,
                PIGMENT_STEPS * PIGMENT_STEPS * self.number_of_facets
            );
        }

        calculate_ressens(self, &summaries)?;

        println!(
            "--- Finished simulation for {} ---\n",
            self.params.species_name
        );

        Ok(())
    }
}
