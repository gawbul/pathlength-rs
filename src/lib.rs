pub mod analysis;
pub mod model;
pub mod parameters;

use anyhow::{Context, Result};
use csv::ReaderBuilder;
use parameters::Parameters;
use std::fs::File;
use std::path::PathBuf;

pub fn parse_input_parameters(filename: &PathBuf) -> Result<Vec<Parameters>> {
    let file = File::open(filename)
        .with_context(|| format!("Failed to open parameter file {:?}", filename))?;
    let mut rdr = ReaderBuilder::new().has_headers(false).from_reader(file);

    let mut params_list = Vec::new();

    for result in rdr.deserialize() {
        let record: Parameters = result.context("Failed to deserialize CSV record")?;
        params_list.push(record);
    }

    if params_list.is_empty() {
        anyhow::bail!("No valid parameter data found in {:?}", filename);
    }

    Ok(params_list)
}
