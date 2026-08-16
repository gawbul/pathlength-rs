use pathlength_rs::analysis::{ABSORPTION_COEFFICIENT, deposit, ring_area, summarise_block};
use pathlength_rs::model::{MAX_PROPAGATION_ANGLE, Model, PATHLENGTHS_HEADER, PIGMENT_STEPS};
use pathlength_rs::parameters::Parameters;
use std::collections::HashMap;
use std::f64::consts::PI;
use std::fs;
use std::path::Path;

/// The reference parameter set used throughout the tests: Nephrops norvegicus,
/// flat lateral measurements.
fn nephrops_flat_lateral(name: &str) -> Parameters {
    Parameters {
        species_name: name.to_string(),
        rhabdom_length: 180.0,
        rhabdom_width: 25.0,
        eye_diameter: 7800.0,
        facet_width: 50.0,
        aperture_diameter: 3200.0,
        cytoplasm_refractive_index: 1.34,
        rhabdom_refractive_index: 1.37,
        blur_circle_extent: 18.0,
        proximal_rhabdom_angle: 0.0,
    }
}

fn read_matrix(filename: &str) -> Vec<Vec<f64>> {
    let content = fs::read_to_string(filename)
        .unwrap_or_else(|e| panic!("Failed to read {}: {}", filename, e));
    let rows: Vec<Vec<f64>> = content
        .trim()
        .split('\n')
        .map(|line| {
            line.split(',')
                .map(|f| f.trim().parse::<f64>().unwrap_or(f64::NAN))
                .collect()
        })
        .collect();
    assert_eq!(
        rows.len(),
        PIGMENT_STEPS,
        "{}: unexpected row count",
        filename
    );
    for row in &rows {
        assert_eq!(
            row.len(),
            PIGMENT_STEPS,
            "{}: unexpected column count",
            filename
        );
    }
    rows
}

fn cleanup(species: &str) {
    for suffix in ["pathlengths", "summary_res", "summary_sen", "debug"] {
        let _ = fs::remove_file(format!("{}_{}.csv", species, suffix));
    }
}

#[test]
fn test_initial_calculations() {
    let mut params = nephrops_flat_lateral("test_eye");
    params.proximal_rhabdom_angle = 12.5;
    let model = Model::new(params).expect("reference parameters must be valid");

    let expected_circumference = PI * 7800.0;
    assert!((model.circumference_of_eye - expected_circumference).abs() < 1e-6);

    let expected_ommatidial_angle = (50.0 / expected_circumference) * 360.0;
    assert!((model.ommatidial_angle - expected_ommatidial_angle).abs() < 1e-6);

    // boa is measured from the rhabdom axis, so the angle at the wall normal is
    // (90 - boa) and light is guided while boa < critical_angle.
    let expected_critical = 90.0 - (1.34 / 1.37f64).asin().to_degrees();
    assert!((model.critical_angle - expected_critical).abs() < 1e-6);
    assert_eq!(model.number_of_facets, 33);
}

/// Covers inputs that previously produced NaNs which then silently disabled the
/// total-internal-reflection test or emptied the simulation, in both cases without
/// any warning to the user.
#[test]
fn test_rejects_unphysical_parameters() {
    /// A named mutation of the reference parameters, and the text its rejection
    /// message must contain.
    type InvalidCase = (&'static str, Box<dyn Fn(&mut Parameters)>, &'static str);

    let cases: Vec<InvalidCase> = vec![
        (
            "cytoplasm index exceeds rhabdom",
            Box::new(|p: &mut Parameters| {
                p.cytoplasm_refractive_index = 1.40;
                p.rhabdom_refractive_index = 1.37;
            }),
            "total internal reflection",
        ),
        (
            "equal refractive indices",
            Box::new(|p: &mut Parameters| {
                p.rhabdom_refractive_index = p.cytoplasm_refractive_index;
            }),
            "total internal reflection",
        ),
        (
            "aperture exceeds eye",
            Box::new(|p: &mut Parameters| p.aperture_diameter = 9000.0),
            "aperture diameter",
        ),
        (
            "zero rhabdom length",
            Box::new(|p: &mut Parameters| p.rhabdom_length = 0.0),
            "rhabdom length",
        ),
        (
            "blur circle below one",
            Box::new(|p: &mut Parameters| p.blur_circle_extent = 0.0),
            "blur circle extent",
        ),
        // Rust's f64 parser accepts "NaN" and "inf", and every ordered comparison
        // against NaN is false, so these used to slip past every range check and
        // produce plausible-looking output: a NaN cytoplasm index gave a NaN critical
        // angle and 83% sensitivity with no warning at all.
        (
            "NaN blur circle",
            Box::new(|p: &mut Parameters| p.blur_circle_extent = f64::NAN),
            "must be a finite number",
        ),
        (
            "NaN cytoplasm index",
            Box::new(|p: &mut Parameters| p.cytoplasm_refractive_index = f64::NAN),
            "must be a finite number",
        ),
        (
            "infinite rhabdom index",
            Box::new(|p: &mut Parameters| p.rhabdom_refractive_index = f64::INFINITY),
            "must be a finite number",
        ),
        (
            "infinite proximal angle",
            Box::new(|p: &mut Parameters| p.proximal_rhabdom_angle = f64::INFINITY),
            "must be a finite number",
        ),
        (
            "negative infinite eye diameter",
            Box::new(|p: &mut Parameters| p.eye_diameter = f64::NEG_INFINITY),
            "must be a finite number",
        ),
        // astacodes ships with an 18-rhabdom blur circle but only 7 facets across
        // the eyeshine patch, leaving 11 rhabdom offsets receiving no light at all.
        (
            "blur circle exceeds facets",
            Box::new(|p: &mut Parameters| {
                p.rhabdom_length = 84.0;
                p.rhabdom_width = 16.0;
                p.eye_diameter = 890.0;
                p.facet_width = 32.0;
                p.aperture_diameter = 445.0;
                p.blur_circle_extent = 18.0;
            }),
            "exceeds the 7 facets",
        ),
    ];

    for (name, mutate, want) in cases {
        let mut params = nephrops_flat_lateral("test_invalid");
        mutate(&mut params);
        match Model::new(params) {
            Ok(m) => panic!(
                "{}: expected an error, got a model with {} facets",
                name, m.number_of_facets
            ),
            Err(err) => assert!(
                err.to_string().contains(want),
                "{}: expected an error mentioning {:?}, got: {}",
                name,
                want,
                err
            ),
        }
    }
}

/// Guards the blur-circle mapping. The previous `facet > fd*i` formulation aliased
/// facets unevenly onto whole rhabdom offsets and skipped an offset entirely wherever
/// fd*i landed on an exact integer, cutting a notch into the profile that the
/// half-maximum search then locked onto.
#[test]
fn test_blur_offset_spans_extent_evenly() {
    let model = Model::new(nephrops_flat_lateral("test_blur")).unwrap();

    assert_eq!(
        model.blur_offset(0),
        0.0,
        "the central facet must be undisplaced"
    );
    let want = model.params.blur_circle_extent - 1.0;
    assert!(
        (model.blur_offset(model.number_of_facets - 1) - want).abs() < 1e-9,
        "the outermost facet must sit at offset {}",
        want
    );

    // The mapping must be strictly increasing with a constant step, so no offset is
    // starved of contributing facets and no facet index sits on a tie boundary.
    let step = model.blur_offset(1) - model.blur_offset(0);
    assert!(step > 0.0, "expected a positive blur step, got {}", step);
    for facet in 1..model.number_of_facets {
        let delta = model.blur_offset(facet) - model.blur_offset(facet - 1);
        assert!(
            (delta - step).abs() < 1e-9,
            "facet {}: expected a uniform blur step of {}, got {}",
            facet,
            step,
            delta
        );
    }

    // A single-rhabdom blur circle means no displacement at all.
    let mut unblurred = nephrops_flat_lateral("test_noblur");
    unblurred.blur_circle_extent = 1.0;
    let plain = Model::new(unblurred).unwrap();
    for facet in 0..plain.number_of_facets {
        assert_eq!(plain.blur_offset(facet), 0.0);
    }
}

/// Guards against the ray tracer folding rays back on themselves. Taking |tan| and
/// |cos| of an angle past 90 degrees used to yield path lengths many times the
/// rhabdom length for a ray that cannot in fact advance towards the proximal end.
#[test]
fn test_rays_stay_within_physical_geometry() {
    for pra in [0.0, 12.5] {
        let mut params = nephrops_flat_lateral("test_geom");
        params.proximal_rhabdom_angle = pra;
        let rhabdom_length = params.rhabdom_length;
        let model = Model::new(params).unwrap();
        let increment = rhabdom_length / 10.0;

        for s_step in 0..PIGMENT_STEPS {
            for t_step in 0..PIGMENT_STEPS {
                for facet in 0..model.number_of_facets {
                    let trace = model.trace_ray(
                        facet,
                        s_step as f64 * increment,
                        t_step as f64 * increment,
                    );
                    assert!(
                        !trace.lost,
                        "facet {} lost the ray at ({},{})",
                        facet, s_step, t_step
                    );
                    assert!(
                        trace.max_angle >= 0.0 && trace.max_angle < MAX_PROPAGATION_ANGLE,
                        "facet {} reached {} deg to the rhabdom axis",
                        facet,
                        trace.max_angle
                    );

                    // Every segment covers some axial depth at an angle no greater
                    // than max_angle, and the ray traverses the rhabdom at most twice
                    // (down, then back off the tapetum), so the total path is bounded
                    // by 2*rhabdom_length/cos(max_angle).
                    let limit = 2.0 * rhabdom_length / trace.max_angle.to_radians().cos() + 1e-9;
                    let total: f64 = trace.pathlengths.iter().sum();
                    for v in &trace.pathlengths {
                        assert!(v.is_finite() && *v >= 0.0, "invalid pathlength {}", v);
                    }
                    assert!(
                        total <= limit,
                        "facet {} at ({},{}) traced {:.1} um, exceeding the {:.1} um bound for a maximum angle of {:.2} deg",
                        facet,
                        s_step,
                        t_step,
                        total,
                        limit,
                        trace.max_angle
                    );
                }
            }
        }
    }
}

/// Confirms that facet transmission is no longer folded into the geometry. The axial
/// ray must traverse exactly the rhabdom length, and twice that when the tapetum
/// reflects it back with no screening pigment in the way.
#[test]
fn test_pathlengths_are_raw_geometry() {
    let params = nephrops_flat_lateral("test_geometry");
    let rhabdom_length = params.rhabdom_length;
    let model = Model::new(params).unwrap();

    let single = model.trace_ray(0, 0.0, 0.0);
    assert_eq!(single.pathlengths.len(), 1);
    assert!((single.pathlengths[0] - rhabdom_length).abs() < 1e-9);

    let reflected = model.trace_ray(0, 0.0, rhabdom_length);
    assert_eq!(reflected.pathlengths.len(), 1);
    assert!((reflected.pathlengths[0] - 2.0 * rhabdom_length).abs() < 1e-9);
}

/// The proximal screening pigment lies in the cytoplasm outside the rhabdom, so it
/// can only absorb light that has already crossed the wall. A totally internally
/// reflected ray stays inside and must not be truncated by it.
#[test]
fn test_guided_ray_ignores_screening_pigment() {
    let mut params = nephrops_flat_lateral("test_guided");
    // A narrow blur circle keeps the entry angle below the critical angle.
    params.blur_circle_extent = 1.0;
    let rhabdom_length = params.rhabdom_length;
    let model = Model::new(params).unwrap();

    let facet = 1;
    let boa = pathlength_rs::model::refracted_angle(facet as f64 * model.ommatidial_angle)
        + model.blur_offset(facet) * model.ommatidial_angle;
    assert!(
        boa < model.critical_angle,
        "test needs a guided ray: boa {} is not below the critical angle {}",
        boa,
        model.critical_angle
    );

    let unscreened = model.trace_ray(facet, 0.0, 0.0);
    let screened = model.trace_ray(facet, rhabdom_length, 0.0);
    assert_eq!(unscreened.pathlengths.len(), 1);
    assert_eq!(screened.pathlengths.len(), 1);
    assert!(
        (unscreened.pathlengths[0] - screened.pathlengths[0]).abs() < 1e-9,
        "screening pigment truncated a guided ray: {} um without, {} um with",
        unscreened.pathlengths[0],
        screened.pathlengths[0]
    );
}

/// Covers the accumulator that used to be a fixed 21-element array, silently
/// discarding every rhabdom past the twenty-first.
#[test]
fn test_deposit_grows_beyond_fixed_array() {
    let mut acc: Vec<f64> = Vec::new();
    deposit(&mut acc, 40, 3.5);
    assert_eq!(acc.len(), 41);
    assert_eq!(acc[40], 3.5);
    deposit(&mut acc, 40, 1.5);
    assert_eq!(acc[40], 5.0);

    // A zero deposit must not extend the profile with meaningless trailing offsets.
    let before = acc.len();
    deposit(&mut acc, 99, 0.0);
    assert_eq!(acc.len(), before);
}

/// Checks the half-maximum interpolation against profiles whose full width at half
/// maximum is known exactly. The stored profile is area-weighted, so each offset is
/// scaled by its ring area and the summary divides it back out.
#[test]
fn test_summarise_block_resolution() {
    let mut params = nephrops_flat_lateral("test_summarise");
    params.blur_circle_extent = 1.0;
    let model = Model::new(params).unwrap();
    let omm = model.ommatidial_angle;

    let cases: Vec<(&str, Vec<f64>, f64)> = vec![
        // Falls from the peak straight to zero: the half maximum lies midway across
        // the first step, so the half width is 0.5 ommatidial angles.
        ("single step", vec![1.0, 0.0], 2.0 * 0.5 * omm),
        // Sits exactly on the half maximum at offset 1, giving a half width of 1.
        (
            "exact half at unit offset",
            vec![1.0, 0.5, 0.0],
            2.0 * 1.0 * omm,
        ),
        // A wider profile crossing half maximum midway between offsets 2 and 3.
        ("wider profile", vec![1.0, 1.0, 1.0, 0.0], 2.0 * 2.5 * omm),
        // A flat top-hat is measured at its own edge: rhabdoms beyond the outermost
        // illuminated one are genuinely dark.
        (
            "top hat measured at its edge",
            vec![1.0, 1.0, 1.0],
            2.0 * 2.5 * omm,
        ),
    ];

    for (name, psf, want) in cases {
        let weighted: Vec<f64> = psf
            .iter()
            .enumerate()
            .map(|(j, v)| v * ring_area(j))
            .collect();
        let got = summarise_block(&model, &weighted);
        assert!(
            (got.fwhm_degrees - want).abs() < 1e-9,
            "{}: expected FWHM {:.6} deg, got {:.6}",
            name,
            want,
            got.fwhm_degrees
        );
        assert_eq!(got.peak_offset, 0, "{}: expected an on-axis peak", name);
    }

    let empty = summarise_block(&model, &[]);
    assert!(empty.fwhm_degrees.is_nan());
    assert_eq!(empty.sensitivity_percent, 0.0);

    // The angular sensitivity function is even about the optic axis, so its width is
    // measured from the axis. Measuring from the peak understates a flat-topped
    // profile by the peak's own offset: here the light is above half maximum out to
    // radius 2.828, so the full width is 5.657 ommatidial angles, not 3.657.
    let flat_top: Vec<f64> = [0.99, 1.0, 0.98, 0.4, 0.0]
        .iter()
        .enumerate()
        .map(|(j, v)| v * ring_area(j))
        .collect();
    let got = summarise_block(&model, &flat_top);
    assert_eq!(got.peak_offset, 1);
    assert!(
        !got.annular,
        "a profile at maximum on the axis is not annular"
    );
    let want = 2.0 * 2.8275862068965516 * omm;
    assert!(
        (got.fwhm_degrees - want).abs() < 1e-9,
        "expected FWHM {:.6} deg measured from the axis, got {:.6}",
        want,
        got.fwhm_degrees
    );

    // A profile that dips below half maximum on the axis is a ring. Its supra-half
    // region does not contain the axis, so there is no acceptance angle: reporting the
    // ring's thickness instead would read as an implausibly sharp eye.
    let annular: Vec<f64> = [0.2, 0.6, 1.0, 0.6, 0.2, 0.0]
        .iter()
        .enumerate()
        .map(|(j, v)| v * ring_area(j))
        .collect();
    let got = summarise_block(&model, &annular);
    assert!(got.annular, "expected the profile to be flagged as annular");
    assert!(
        got.fwhm_degrees.is_nan(),
        "expected an undefined FWHM for an annular profile, got {}",
        got.fwhm_degrees
    );
    // Sensitivity is independent of the resolution classification.
    assert!(got.sensitivity_percent > 0.0);
}

#[test]
fn test_run_simulation_produces_well_formed_blocks() {
    let mut model_flat = Model::new(nephrops_flat_lateral("test_flat")).unwrap();
    let mut pointy_params = nephrops_flat_lateral("test_pointy");
    pointy_params.proximal_rhabdom_angle = 12.5;
    let mut model_pointy = Model::new(pointy_params).unwrap();

    let facets = model_flat.number_of_facets;
    model_flat.run_simulation().unwrap();
    model_pointy.run_simulation().unwrap();

    let flat = fs::read_to_string("test_flat_pathlengths.csv").unwrap();
    let pointy = fs::read_to_string("test_pointy_pathlengths.csv").unwrap();

    // The file is a plain rectangular CSV: a header, then one row per rhabdom, each
    // carrying its own keys. There is no block terminator and no positional state.
    let mut lines = flat.lines().map(str::trim).filter(|l| !l.is_empty());
    assert_eq!(lines.next(), Some(PATHLENGTHS_HEADER));
    assert!(
        !flat.contains("\n999\n"),
        "expected no 999 block terminator in the output"
    );

    // Every (block, facet) pair must appear exactly once as a group, with rhabdom
    // indices running from zero.
    let mut seen: HashMap<(usize, usize), usize> = HashMap::new();
    for line in lines {
        let fields: Vec<&str> = line.split(',').collect();
        assert_eq!(fields.len(), 6, "expected 6 fields, got {:?}", line);
        for field in &fields {
            field
                .parse::<f64>()
                .unwrap_or_else(|_| panic!("non-numeric field {:?}", field));
        }
        let block: usize = fields[0].parse().unwrap();
        let facet: usize = fields[3].parse().unwrap();
        let rhabdom: usize = fields[4].parse().unwrap();
        let next = seen.entry((block, facet)).or_insert(0);
        assert_eq!(
            rhabdom, *next,
            "block {} facet {}: rhabdom indices must run from zero",
            block, facet
        );
        *next += 1;
    }
    assert_eq!(seen.len(), PIGMENT_STEPS * PIGMENT_STEPS * facets);
    assert_ne!(
        flat, pointy,
        "the proximal rhabdom angle must change the traced pathlengths"
    );

    cleanup("test_flat");
    cleanup("test_pointy");
}

#[test]
fn test_debug_flag_output() {
    let mut model_no_debug = Model::new(nephrops_flat_lateral("test_nodebug")).unwrap();
    model_no_debug.debug_mode = false;
    model_no_debug.run_simulation().unwrap();
    assert!(
        !Path::new("test_nodebug_debug.csv").exists(),
        "expected no debug file when debug_mode is false"
    );
    cleanup("test_nodebug");

    let mut model_debug = Model::new(nephrops_flat_lateral("test_debug")).unwrap();
    let facets = model_debug.number_of_facets;
    model_debug.debug_mode = true;
    model_debug.run_simulation().unwrap();

    let debug = fs::read_to_string("test_debug_debug.csv").unwrap();
    let lines: Vec<&str> = debug.lines().filter(|l| !l.trim().is_empty()).collect();
    // One header plus one row per traced ray, rather than the block headings alone
    // that earlier versions emitted.
    assert_eq!(lines.len(), 1 + PIGMENT_STEPS * PIGMENT_STEPS * facets);
    assert!(lines[0].starts_with("block,shielding_um,tapetal_um,facet,"));
    assert_eq!(lines[1].split(',').count(), 12);

    cleanup("test_debug");
}

/// Runs the reference eye end to end and checks the resolution matrix is fully
/// defined. Earlier versions passed a `value > 0` check even when every cell was the
/// fabricated no-crossing fallback.
#[test]
fn test_summary_matrices_are_usable() {
    let mut model = Model::new(nephrops_flat_lateral("test_matrix")).unwrap();
    model.run_simulation().unwrap();

    let res = read_matrix("test_matrix_summary_res.csv");
    let sens = read_matrix("test_matrix_summary_sen.csv");

    for row in 0..PIGMENT_STEPS {
        for col in 0..PIGMENT_STEPS {
            assert!(
                res[row][col].is_finite() && res[row][col] > 0.0,
                "resolution [{}][{}] = {}, expected a positive acceptance angle",
                row,
                col,
                res[row][col]
            );
            assert!(
                (0.0..=100.0).contains(&sens[row][col]),
                "sensitivity [{}][{}] = {}, expected a percentage in 0-100",
                row,
                col,
                sens[row][col]
            );
        }
    }

    // Migrating the screening pigment across the whole rhabdom must not raise
    // sensitivity: a light-adapted eye absorbs less than a dark-adapted one.
    assert!(sens[PIGMENT_STEPS - 1][0] <= sens[0][0]);
    // Extending the tapetum must not lower sensitivity: reflected light is absorbed twice.
    assert!(sens[0][PIGMENT_STEPS - 1] >= sens[0][0]);

    cleanup("test_matrix");
}

/// The whole aperture is a single axial ray, so the absorbed percentage is exactly
/// the Beer-Lambert absorbance over the rhabdom length and the area weight cancels
/// against the patch area.
#[test]
fn test_single_facet_sensitivity_is_beer_lambert() {
    let params = Parameters {
        species_name: "test_single".to_string(),
        rhabdom_length: 100.0,
        rhabdom_width: 20.0,
        eye_diameter: 1000.0,
        facet_width: 20.0,
        aperture_diameter: 40.0,
        cytoplasm_refractive_index: 1.34,
        rhabdom_refractive_index: 1.37,
        blur_circle_extent: 1.0,
        proximal_rhabdom_angle: 0.0,
    };
    let mut model = Model::new(params).unwrap();
    assert_eq!(model.number_of_facets, 1, "test needs a single-facet patch");
    let omm = model.ommatidial_angle;
    model.run_simulation().unwrap();

    let sens = read_matrix("test_single_summary_sen.csv");
    let want = 100.0 * (1.0 - (-ABSORPTION_COEFFICIENT * 100.0f64).exp());
    assert!(
        (sens[0][0] - want).abs() < 1e-4,
        "expected {:.4}% absorbed, got {:.4}%",
        want,
        sens[0][0]
    );

    // All the light lands on the axial rhabdom and the next offset out is dark, so
    // the half maximum falls midway across that single step.
    let res = read_matrix("test_single_summary_res.csv");
    assert!(
        (res[0][0] - omm).abs() < 1e-4,
        "expected an acceptance angle of {:.4} deg, got {:.4}",
        omm,
        res[0][0]
    );

    cleanup("test_single");
}
