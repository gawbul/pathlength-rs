use crate::model::{Model, PIGMENT_STEPS};
use anyhow::{Context, Result, bail};
use std::f64::consts::PI;
use std::fs::File;
use std::io::{BufWriter, Write};

/// The rhabdom absorption coefficient in um^-1, used in the Beer-Lambert absorbance
/// 1 - exp(-k*L). Reported values for crustacean rhabdoms span roughly 0.0067 to
/// 0.01 um^-1.
pub const ABSORPTION_COEFFICIENT: f64 = 0.01;

/// The area, in squared facet widths, of the annulus of rhabdoms lying at the given
/// whole-rhabdom offset from the optic axis. This is the same measure used to weight
/// the contributing facets, so dividing an area-weighted total by it yields a radial
/// area density.
pub fn ring_area(offset: usize) -> f64 {
    if offset == 0 {
        return PI * 0.25;
    }
    let outer = PI * (offset as f64 + 0.5).powi(2);
    let inner = PI * (offset as f64 - 0.5).powi(2);
    outer - inner
}

/// Adds an amount at the given offset, growing the profile as required. Earlier
/// versions used a fixed 21-element array and silently discarded everything beyond
/// it, which lost up to a quarter of the absorbed light for widely blurred eyes.
pub fn deposit(dst: &mut Vec<f64>, offset: usize, amount: f64) {
    if amount == 0.0 {
        return;
    }
    if dst.len() <= offset {
        dst.resize(offset + 1, 0.0);
    }
    dst[offset] += amount;
}

/// The resolution and sensitivity derived from one pigment block.
#[derive(Debug, Clone, Copy)]
pub struct BlockSummary {
    /// The acceptance angle: the full width at half maximum of the angular sensitivity
    /// function, in degrees. NaN when the profile carries no light or is annular, in
    /// which case there is no acceptance angle to report.
    pub fwhm_degrees: f64,
    /// Percentage of incident light absorbed, averaged over the eyeshine patch (0-100).
    pub sensitivity_percent: f64,
    /// The rhabdom offset carrying the most light.
    pub peak_offset: usize,
    /// Set when the profile dips below half its maximum on the optic axis, so the
    /// light forms a ring rather than a central spot.
    pub annular: bool,
}

/// Converts one block's area-weighted absorption profile into resolution and
/// sensitivity.
pub fn summarise_block(model: &Model, rhabdoms: &[f64]) -> BlockSummary {
    let mut out = BlockSummary {
        fwhm_degrees: f64::NAN,
        sensitivity_percent: 0.0,
        peak_offset: 0,
        annular: false,
    };

    // Sensitivity: the area-weighted mean of the absorbed percentage over the
    // eyeshine patch. The facet weights telescope to exactly pi*(N-0.5)^2, so
    // dividing by that area makes this a true weighted mean in the range 0-100.
    let total: f64 = rhabdoms.iter().sum();
    let patch_area = PI * (model.number_of_facets as f64 - 0.5).powi(2);
    if patch_area > 0.0 {
        out.sensitivity_percent = total / patch_area;
    }

    if rhabdoms.is_empty() {
        return out;
    }

    // Point spread function: light per unit area at each rhabdom offset. The
    // contributing facets are weighted by their source annulus, so the light arriving
    // at an offset must be divided by the annulus it is spread over to recover an
    // intensity. Without this the profile rises monotonically with offset simply
    // because outer annuli contain more ommatidia.
    //
    // The profile ends at the outermost rhabdom that receives any light, so the next
    // offset out is genuinely dark. Including that zero captures the falling edge of
    // the blur circle, which is where a top-hat profile crosses its half maximum.
    let mut psf = vec![0.0f64; rhabdoms.len() + 1];
    for (j, &v) in rhabdoms.iter().enumerate() {
        psf[j] = v / ring_area(j);
    }

    let mut peak = 0usize;
    for (j, &v) in psf.iter().enumerate() {
        if v > psf[peak] {
            peak = j;
        }
    }
    out.peak_offset = peak;
    if psf[peak] <= 0.0 {
        return out;
    }
    let half = psf[peak] / 2.0;

    // A profile that is already below half its maximum on the optic axis is annular:
    // the light forms a ring, and the region above half maximum is a band that does
    // not contain the axis. There is no acceptance angle to report, so the width is
    // left undefined rather than substituting the ring's thickness for it.
    if psf[0] < half {
        out.annular = true;
        return out;
    }

    // The angular sensitivity function is even about the optic axis - offset j stands
    // for both +j and -j - so its full width at half maximum is twice the radius at
    // which it first falls below half. That radius is measured from the axis, not from
    // the peak: on a flat-topped profile whose maximum sits slightly off-axis,
    // measuring from the peak would understate the width by the peak's own offset.
    for i in 0..psf.len() - 1 {
        if psf[i] >= half && psf[i + 1] < half {
            let frac = (psf[i] - half) / (psf[i] - psf[i + 1]);
            out.fwhm_degrees = 2.0 * (i as f64 + frac) * model.ommatidial_angle;
            break;
        }
    }
    out
}

/// Adds one facet's traced ray into the area-weighted absorption profile, which
/// records how much light reaches each whole-rhabdom offset from the optic axis.
pub fn accumulate(model: &Model, profile: &mut Vec<f64>, facet_index: usize, pathlengths: &[f64]) {
    // Light gathered by this facet, and the rhabdom offset its image lands on.
    let transmission = model.facet_transmission(facet_index);
    let source_area = ring_area(facet_index);
    let offset = model.blur_offset(facet_index);
    let base = offset.floor() as usize;
    let frac = offset - base as f64;

    let mut tot = 0.0f64;
    for (rhabdom, &pathlength) in pathlengths.iter().enumerate() {
        if pathlength <= 0.0 {
            continue;
        }
        // Fraction of the light still travelling that this rhabdom absorbs.
        let absorbed = (1.0 - tot) * (1.0 - (-ABSORPTION_COEFFICIENT * pathlength).exp());
        tot += absorbed;
        // Facet transmission attenuates the flux entering the eye; it does not shorten
        // the geometric path, so it multiplies the absorbed intensity rather than the
        // exponent.
        let weighted = 100.0 * transmission * absorbed * source_area;

        // The blur displacement is continuous, so split the light between the two
        // rhabdom offsets that bracket it.
        deposit(profile, base + rhabdom, weighted * (1.0 - frac));
        if frac > 0.0 {
            deposit(profile, base + rhabdom + 1, weighted * frac);
        }
    }
}

/// Writes the resolution and sensitivity matrices for the pigment states accumulated
/// during the simulation.
pub fn calculate_ressens(model: &Model, summaries: &[BlockSummary]) -> Result<()> {
    let species_name = &model.params.species_name;
    println!("INFO: Calculating resolution and sensitivity...");

    if summaries.len() != PIGMENT_STEPS * PIGMENT_STEPS {
        bail!(
            "expected {} pigment states, got {}",
            PIGMENT_STEPS * PIGMENT_STEPS,
            summaries.len()
        );
    }

    write_summary_matrix(
        &format!("{}_summary_res.csv", species_name),
        summaries,
        |b| b.fwhm_degrees,
    )?;
    write_summary_matrix(
        &format!("{}_summary_sen.csv", species_name),
        summaries,
        |b| b.sensitivity_percent,
    )?;

    let annular = summaries.iter().filter(|b| b.annular).count();
    let dark = summaries
        .iter()
        .filter(|b| !b.annular && b.fwhm_degrees.is_nan())
        .count();
    if annular > 0 {
        println!(
            "WARNING: {} of {} pigment states have an annular profile, with the light forming a ring rather than a central spot; they have no acceptance angle and are reported as NaN.",
            annular,
            summaries.len()
        );
    }
    if dark > 0 {
        println!(
            "WARNING: {} of {} pigment states absorb no light; their resolution is reported as NaN.",
            dark,
            summaries.len()
        );
    }

    Ok(())
}

/// Writes an 11x11 matrix with shielding pigment position varying down the rows and
/// tapetal pigment position across the columns.
fn write_summary_matrix<F>(filename: &str, summaries: &[BlockSummary], value: F) -> Result<()>
where
    F: Fn(BlockSummary) -> f64,
{
    let file = File::create(filename).with_context(|| format!("Failed to create {}", filename))?;
    let mut writer = BufWriter::new(file);

    for row in 0..PIGMENT_STEPS {
        let cells: Vec<String> = (0..PIGMENT_STEPS)
            .map(|col| format!("{:.4}", value(summaries[row * PIGMENT_STEPS + col])))
            .collect();
        writeln!(writer, "{}", cells.join(","))?;
    }
    writer.flush()?;
    Ok(())
}
