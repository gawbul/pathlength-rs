use anyhow::{Result, bail};
use serde::Deserialize;

#[derive(Debug, Clone, Deserialize)]
pub struct Parameters {
    pub species_name: String,
    pub rhabdom_length: f64,
    pub rhabdom_width: f64,
    pub eye_diameter: f64,
    pub facet_width: f64,
    pub aperture_diameter: f64,
    pub cytoplasm_refractive_index: f64,
    pub rhabdom_refractive_index: f64,
    pub blur_circle_extent: f64,
    pub proximal_rhabdom_angle: f64,
}

impl Parameters {
    /// Rejects inputs that would produce NaNs or nonsensical optics.
    ///
    /// Without these checks an aperture wider than the eye yields `sqrt` of a
    /// negative number and an empty simulation, while a cytoplasm refractive index
    /// at or above the rhabdom's yields a NaN critical angle that silently makes
    /// every total-internal-reflection test evaluate to false.
    pub fn validate(&self) -> Result<()> {
        if self.species_name.trim().is_empty() {
            bail!("species name is required");
        }
        for (name, value) in [
            ("rhabdom length", self.rhabdom_length),
            ("rhabdom width", self.rhabdom_width),
            ("eye diameter", self.eye_diameter),
            ("facet width", self.facet_width),
            ("aperture diameter", self.aperture_diameter),
        ] {
            // Spelled out rather than negated so that NaN is rejected explicitly.
            if value.is_nan() || value <= 0.0 {
                bail!("{} must be greater than 0 um, got {}", name, value);
            }
        }
        if self.aperture_diameter >= self.eye_diameter {
            bail!(
                "aperture diameter ({} um) must be smaller than eye diameter ({} um)",
                self.aperture_diameter,
                self.eye_diameter
            );
        }
        if self.cytoplasm_refractive_index <= 1.0 {
            bail!(
                "cytoplasm refractive index must be greater than 1.0, got {}",
                self.cytoplasm_refractive_index
            );
        }
        if self.rhabdom_refractive_index <= self.cytoplasm_refractive_index {
            bail!(
                "rhabdom refractive index ({}) must exceed cytoplasm refractive index ({}) for total internal reflection",
                self.rhabdom_refractive_index,
                self.cytoplasm_refractive_index
            );
        }
        if self.blur_circle_extent < 1.0 {
            bail!(
                "blur circle extent must be at least 1 rhabdom, got {}",
                self.blur_circle_extent
            );
        }
        if self.proximal_rhabdom_angle < 0.0 {
            bail!(
                "proximal rhabdom angle must not be negative, got {}",
                self.proximal_rhabdom_angle
            );
        }
        Ok(())
    }
}
