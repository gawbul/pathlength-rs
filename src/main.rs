mod analysis;
mod model;
mod parameters;

use anyhow::{Context, Result};
use clap::Parser;
use csv::ReaderBuilder;
use model::Model;
use parameters::Parameters;
use std::fs::File;
use std::path::PathBuf;

const VERSION: &str = "0.6.0";

#[derive(Parser)]
#[command(name = "pathlength")]
#[command(version = VERSION)]
#[command(
    about = "Calculates resolution and sensitivity in reflective superposition compound eyes"
)]
struct Cli {
    /// Path to a parameter file (CSV format). (Required)
    #[arg(short = 'f')]
    filename: Option<PathBuf>,

    /// Show the program citation.
    #[arg(short = 'c')]
    citation: bool,

    /// Show the program license.
    #[arg(short = 'l')]
    license: bool,

    /// Show program version.
    #[arg(short = 'v')]
    show_version: bool,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    if cli.license {
        print_license();
        return Ok(());
    }

    if cli.citation {
        print_citation();
        return Ok(());
    }

    if cli.show_version {
        println!("pathlength version {}", VERSION);
        return Ok(());
    }

    let filename = match cli.filename {
        Some(f) => f,
        None => {
            use clap::CommandFactory;
            let mut cmd = Cli::command();
            cmd.print_help()?;
            println!("\nError: No parameter file supplied. Use the -f flag to specify a file.");
            std::process::exit(1);
        }
    };

    println!("Parsing input parameters from {:?}...", filename);
    let params_list = parse_input_parameters(&filename)?;

    for params in params_list {
        let mut model = Model::new(params);
        model.run_simulation()?;
    }

    println!("All simulations complete.");

    Ok(())
}

fn parse_input_parameters(filename: &PathBuf) -> Result<Vec<Parameters>> {
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

fn print_citation() {
    println!(
        r#"Gaten, E., Moss, S., Johnson, M. 2013. The Reniform Reflecting Superposition Compound Eyes of Nephrops Norvegicus:
Optics, Susceptibility to Light-Induced Damage, Electrophysiology and a Ray Tracing Model. In: M. L. Johnson and M. P. Johnson, ed(s).
Advances in Marine Biology: The Ecology and Biology of Nephrops norvegicus. Oxford: Academic Press, 107:148.
"#
    );
}

fn print_license() {
    println!(
        r#"pathlength - calculates resolution and sensitivity in reflective superposition compound eyes.

Copyright (C) 2025 Dr Stephen P Moss

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>
"#
    );
}
