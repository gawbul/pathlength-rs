use anyhow::Result;
use clap::Parser;
use pathlength_rs::model::Model;
use pathlength_rs::parse_input_parameters;
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

    /// Generate debug CSV output file.
    #[arg(short = 'd')]
    debug: bool,

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

    let total = params_list.len();
    let mut failed = 0usize;
    for params in params_list {
        let species_name = params.species_name.clone();
        let mut model = match Model::new(params) {
            Ok(model) => model,
            Err(err) => {
                eprintln!("Skipping {}: {}", species_name, err);
                failed += 1;
                continue;
            }
        };
        model.debug_mode = cli.debug;
        if let Err(err) = model.run_simulation() {
            eprintln!("Simulation for {} failed: {}", species_name, err);
            failed += 1;
        }
    }

    if failed > 0 {
        anyhow::bail!(
            "{} of {} parameter sets could not be simulated.",
            failed,
            total
        );
    }
    println!("All simulations complete.");

    Ok(())
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
