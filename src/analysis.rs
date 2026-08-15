use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

pub fn calculate_ressens(species_name: &str, ommatidial_angle: f64) -> Result<()> {
    println!("INFO: Calculating resolution and sensitivity...");

    let pathlengths_filename = format!("{}_pathlengths.csv", species_name);
    let res_filename = format!("{}_summary_res.csv", species_name);
    let sens_filename = format!("{}_summary_sen.csv", species_name);

    let file = File::open(&pathlengths_filename)
        .with_context(|| format!("Failed to open {}", pathlengths_filename))?;
    let reader = BufReader::new(file);

    let mut res_file = File::create(&res_filename)?;
    let mut sens_file = File::create(&sens_filename)?;

    let mut rhabdoms = vec![0.0f64; 21];
    let mut intensities = vec![0.0f64; 21];
    let mut matrix_sens: Vec<String> = Vec::new();
    let mut matrix_res: Vec<String> = Vec::new();

    let mut facet: f64 = 0.0;
    let mut arem = 0.0;
    let mut cc = 0;
    let mut dd = 0;
    let mut header_count = 0;

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();

        if line.is_empty() {
            continue;
        }

        if line == "999" {
            // End of block, summarize
            let mut sens = 0.0;
            for val in &rhabdoms {
                sens += val;
            }

            // Find true peak of the intensity profile
            let mut max_intensity = intensities[0];
            for &v in &intensities {
                if v > max_intensity {
                    max_intensity = v;
                }
            }
            let halfway_point = max_intensity / 2.0;

            let mut optic_axis = 0.0;
            let mut xz = intensities[0];
            let mut yy = intensities[1];

            // Scan full array to find threshold crossing
            let mut found_crossing = false;
            for i in 0..(intensities.len() - 1) {
                if intensities[i] >= halfway_point && intensities[i + 1] < halfway_point {
                    xz = intensities[i];
                    yy = intensities[i + 1];
                    optic_axis = ommatidial_angle * (i as f64);
                    found_crossing = true;
                    break;
                }
            }

            let frac = if !found_crossing {
                // If it never drops below half-max (extreme blur), cap at the maximum measured angle
                optic_axis = ommatidial_angle * ((intensities.len() - 1) as f64);
                0.0
            } else {
                let diff = xz - yy;
                let hwp = xz - halfway_point;
                if diff > 0.0 {
                    let mut f = hwp / diff;
                    // strictly bound frac between 0.0 and 1.0 to guarantee correct interpolation
                    f = f.clamp(0.0, 1.0);
                    f
                } else {
                    0.0
                }
            };

            let oab = frac * ommatidial_angle;
            let res = oab + optic_axis;

            if cc == 0 && dd > 0 {
                writeln!(sens_file, "{}", matrix_sens.join(","))?;
                writeln!(res_file, "{}", matrix_res.join(","))?;
                matrix_sens.clear();
                matrix_res.clear();
            }

            if arem > 0.0 {
                matrix_sens.push(format!("{}", (sens / arem) as i64));
            } else {
                matrix_sens.push("0".to_string());
            }

            matrix_res.push(format!("{}", (res * 200.0) as i64));

            cc += 1;
            if cc == 11 {
                dd += 1;
                cc = 0;
            }

            // Reset for next block
            rhabdoms.fill(0.0);
            intensities.fill(0.0);
            facet = 0.0;
            header_count = 0;
        } else if header_count < 2 {
            // Shielding or Tapetal pigment header line
            header_count += 1;
        } else {
            // Facet pathlength row
            let parts: Vec<&str> = line.split(',').collect();
            let mut tot = 0.0;

            let area = std::f64::consts::PI * (facet + 0.5).powi(2);
            let mut inci = std::f64::consts::PI * (facet - 0.5).powi(2);
            if facet == 0.0 {
                inci = 0.0;
            }
            let torus = area - inci;
            if area > arem {
                arem = area;
            }

            for (rhabdom_idx, part) in parts.iter().enumerate() {
                let pathlength: f64 = part.parse().unwrap_or(0.0);
                let absorbance = if pathlength > 0.0 {
                    1.0 - (-0.01 * pathlength).exp()
                } else {
                    0.0
                };

                let mut bx = 0.0;
                if rhabdom_idx == 0 && absorbance > 0.0 {
                    bx = 100.0 * absorbance;
                } else if rhabdom_idx > 0 && absorbance > 0.0 {
                    bx = 100.0 * ((1.0 - tot) * absorbance);
                }

                if absorbance == 0.0 {
                    bx = 0.0;
                }

                tot += bx / 100.0;

                if rhabdom_idx < rhabdoms.len() {
                    intensities[rhabdom_idx] += bx; // Intensity profile (no torus weighting)
                    rhabdoms[rhabdom_idx] += bx * torus; // Total Sensitivity (area weighted)
                }
            }
            facet += 1.0;
        }
    }

    // Write final line
    if !matrix_sens.is_empty() {
        writeln!(sens_file, "{}", matrix_sens.join(","))?;
        writeln!(res_file, "{}", matrix_res.join(","))?;
    }

    Ok(())
}
