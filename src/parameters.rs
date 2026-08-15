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
    /// Validates and corrects parameters after loading
    pub fn post_process(&mut self) {
        // Ensure blur circle is at least 1.0
        if self.blur_circle_extent < 1.0 {
            self.blur_circle_extent = 1.0;
        }
    }
}
