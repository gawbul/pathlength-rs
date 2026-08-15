use pathlength_rs::model::Model;
use pathlength_rs::parameters::Parameters;
use std::f64::consts::PI;
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[test]
fn test_initial_calculations() {
    let params = Parameters {
        species_name: "test_eye".to_string(),
        rhabdom_length: 180.0,
        rhabdom_width: 25.0,
        eye_diameter: 7800.0,
        facet_width: 50.0,
        aperture_diameter: 3200.0,
        cytoplasm_refractive_index: 1.34,
        rhabdom_refractive_index: 1.37,
        blur_circle_extent: 18.0,
        proximal_rhabdom_angle: 12.5,
    };

    let model = Model::new(params);

    let expected_circumference = PI * 7800.0;
    assert!((model.circumference_of_eye - expected_circumference).abs() < 1e-6);

    let expected_ommatidial_angle = (50.0 / expected_circumference) * 360.0;
    assert!((model.ommatidial_angle - expected_ommatidial_angle).abs() < 1e-6);

    let expected_snell = (1.34 / 1.37f64).asin().to_degrees();
    let expected_critical = 90.0 - expected_snell;
    assert!((model.critical_angle - expected_critical).abs() < 1e-6);
    assert!(model.number_of_facets > 0);
}

#[test]
fn test_run_model_blur_circle_and_pointy_rhabdom() {
    let params_flat = Parameters {
        species_name: "test_flat".to_string(),
        rhabdom_length: 180.0,
        rhabdom_width: 25.0,
        eye_diameter: 7800.0,
        facet_width: 50.0,
        aperture_diameter: 3200.0,
        cytoplasm_refractive_index: 1.34,
        rhabdom_refractive_index: 1.37,
        blur_circle_extent: 18.0,
        proximal_rhabdom_angle: 0.0,
    };

    let params_pointy = Parameters {
        species_name: "test_pointy".to_string(),
        rhabdom_length: 180.0,
        rhabdom_width: 25.0,
        eye_diameter: 7800.0,
        facet_width: 50.0,
        aperture_diameter: 3200.0,
        cytoplasm_refractive_index: 1.34,
        rhabdom_refractive_index: 1.37,
        blur_circle_extent: 18.0,
        proximal_rhabdom_angle: 12.5,
    };

    let mut model_flat = Model::new(params_flat);
    let mut model_pointy = Model::new(params_pointy);

    model_flat.run_simulation().unwrap();
    model_pointy.run_simulation().unwrap();

    let flat_pathlengths = "test_flat_pathlengths.csv";
    let pointy_pathlengths = "test_pointy_pathlengths.csv";

    // Read flat file
    let file_flat = fs::File::open(flat_pathlengths).unwrap();
    let lines_flat: Vec<String> = BufReader::new(file_flat)
        .lines()
        .map(|l| l.unwrap())
        .collect();

    let mut found_leading_zero = false;
    let mut block_count = 0;
    for line in &lines_flat {
        if line == "999" {
            block_count += 1;
        }
        if line.starts_with("0,") {
            found_leading_zero = true;
        }
    }

    assert!(
        found_leading_zero,
        "Expected leading zeros in blur circle facets"
    );
    assert_eq!(block_count, 121, "Expected 121 pigment blocks (11x11)");

    // Read pointy file
    let file_pointy = fs::File::open(pointy_pathlengths).unwrap();
    let lines_pointy: Vec<String> = BufReader::new(file_pointy)
        .lines()
        .map(|l| l.unwrap())
        .collect();

    let mut has_difference = false;
    for (f, p) in lines_flat.iter().zip(lines_pointy.iter()) {
        if f != p {
            has_difference = true;
            break;
        }
    }

    assert!(
        has_difference,
        "Expected differences between flat and pointy rhabdoms"
    );

    // Clean up test files
    let _ = fs::remove_file("test_flat_pathlengths.csv");
    let _ = fs::remove_file("test_flat_summary_res.csv");
    let _ = fs::remove_file("test_flat_summary_sen.csv");
    let _ = fs::remove_file("test_pointy_pathlengths.csv");
    let _ = fs::remove_file("test_pointy_summary_res.csv");
    let _ = fs::remove_file("test_pointy_summary_sen.csv");
}

#[test]
fn test_debug_flag_output() {
    let params_no_debug = Parameters {
        species_name: "test_nodebug".to_string(),
        rhabdom_length: 100.0,
        rhabdom_width: 20.0,
        eye_diameter: 1000.0,
        facet_width: 20.0,
        aperture_diameter: 500.0,
        cytoplasm_refractive_index: 1.34,
        rhabdom_refractive_index: 1.37,
        blur_circle_extent: 1.0,
        proximal_rhabdom_angle: 0.0,
    };
    let mut model_no_debug = Model::new(params_no_debug);
    model_no_debug.debug_mode = false;
    model_no_debug.run_simulation().unwrap();

    let _ = fs::remove_file("test_nodebug_pathlengths.csv");
    let _ = fs::remove_file("test_nodebug_summary_res.csv");
    let _ = fs::remove_file("test_nodebug_summary_sen.csv");
    assert!(
        !Path::new("test_nodebug_debug.csv").exists(),
        "Expected test_nodebug_debug.csv to NOT exist when debug_mode is false"
    );

    let params_debug = Parameters {
        species_name: "test_debug".to_string(),
        rhabdom_length: 100.0,
        rhabdom_width: 20.0,
        eye_diameter: 1000.0,
        facet_width: 20.0,
        aperture_diameter: 500.0,
        cytoplasm_refractive_index: 1.34,
        rhabdom_refractive_index: 1.37,
        blur_circle_extent: 1.0,
        proximal_rhabdom_angle: 0.0,
    };
    let mut model_debug = Model::new(params_debug);
    model_debug.debug_mode = true;
    model_debug.run_simulation().unwrap();

    let _ = fs::remove_file("test_debug_pathlengths.csv");
    let _ = fs::remove_file("test_debug_summary_res.csv");
    let _ = fs::remove_file("test_debug_summary_sen.csv");
    assert!(
        Path::new("test_debug_debug.csv").exists(),
        "Expected test_debug_debug.csv to exist when debug_mode is true"
    );
    let _ = fs::remove_file("test_debug_debug.csv");
}

#[test]
fn test_calculate_ressens_full_matrix() {
    let params = Parameters {
        species_name: "test_matrix".to_string(),
        rhabdom_length: 127.0,
        rhabdom_width: 15.8,
        eye_diameter: 2480.0,
        facet_width: 22.5,
        aperture_diameter: 870.0,
        cytoplasm_refractive_index: 1.34,
        rhabdom_refractive_index: 1.37,
        blur_circle_extent: 1.0,
        proximal_rhabdom_angle: 0.0,
    };

    let mut model = Model::new(params);
    model.run_simulation().unwrap();

    let res_content = fs::read_to_string("test_matrix_summary_res.csv").unwrap();
    let res_lines: Vec<&str> = res_content.trim().split('\n').collect();
    assert_eq!(
        res_lines.len(),
        11,
        "Expected 11 rows in resolution summary file"
    );

    for (row_idx, line) in res_lines.iter().enumerate() {
        let parts: Vec<&str> = line.split(',').collect();
        assert_eq!(parts.len(), 11, "Row {}: expected 11 columns", row_idx);
        for (col_idx, val_str) in parts.iter().enumerate() {
            let val: i64 = val_str.trim().parse().expect("Valid integer");
            assert!(
                val > 0,
                "Row {}, Col {}: expected positive resolution value, got {}",
                row_idx,
                col_idx,
                val
            );
        }
    }

    // Clean up
    let _ = fs::remove_file("test_matrix_pathlengths.csv");
    let _ = fs::remove_file("test_matrix_summary_res.csv");
    let _ = fs::remove_file("test_matrix_summary_sen.csv");
}

#[test]
fn test_calculate_ressens_nephrops_wide_blur_circle() {
    let params = Parameters {
        species_name: "test_nephrops_bce18".to_string(),
        rhabdom_length: 180.0,
        rhabdom_width: 25.0,
        eye_diameter: 7800.0,
        facet_width: 50.0,
        aperture_diameter: 3200.0,
        cytoplasm_refractive_index: 1.34,
        rhabdom_refractive_index: 1.37,
        blur_circle_extent: 18.0,
        proximal_rhabdom_angle: 0.0,
    };

    let mut model = Model::new(params);
    model.run_simulation().unwrap();

    let res_content = fs::read_to_string("test_nephrops_bce18_summary_res.csv").unwrap();
    let res_lines: Vec<&str> = res_content.trim().split('\n').collect();
    assert_eq!(
        res_lines.len(),
        11,
        "Expected 11 rows in resolution summary file"
    );

    for (row_idx, line) in res_lines.iter().enumerate() {
        let parts: Vec<&str> = line.split(',').collect();
        assert_eq!(parts.len(), 11, "Row {}: expected 11 columns", row_idx);
        for (col_idx, val_str) in parts.iter().enumerate() {
            let val: i64 = val_str.trim().parse().expect("Valid integer");
            assert!(
                val > 0,
                "Row {}, Col {}: expected positive resolution value, got {}",
                row_idx,
                col_idx,
                val
            );
        }
    }

    // Clean up
    let _ = fs::remove_file("test_nephrops_bce18_pathlengths.csv");
    let _ = fs::remove_file("test_nephrops_bce18_summary_res.csv");
    let _ = fs::remove_file("test_nephrops_bce18_summary_sen.csv");
}
