# Superposition Eye Pathlength Program

PathLength implements a ray tracing model to calculate the resolution and sensitivity of reflective superposition compound eyes.

Original QBASIC version by Dr Magnus L Johnson and Genevre Parker, 1995

Golang rewrite by Dr Stephen P Moss, 2025

Rust rewrite by Dr Stephen P Moss, 2025

Author: Dr Stephen P Moss

Website: [https://www.gawbul.io](https://www.gawbul.io)

Email: gawbul@gmail.com

## Install Rust compiler

### macOS / Linux / Windows

The recommended way to install Rust is via `rustup`.

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Follow the on-screen instructions to complete the installation. You may need to restart your terminal or run `source $HOME/.cargo/env`.

## Checkout source code

```bash
# Create and change to projects directory
mkdir -p ~/projects
cd ~/projects

# Clone the GitHub repository
git clone git@github.com:gawbul/pathlength-rs.git
```

## Usage

### Run the program from source

```bash
cd ~/projects/pathlength-rs
cargo run --release -- -f example_data/acanthephyra_parameters.txt
```

### Run the test suite

```bash
cd ~/projects/pathlength-rs
cargo test
```

### Build the program

To create an optimized release binary:

```bash
cd ~/projects/pathlength-rs
cargo build --release
```

The binary will be located at `target/release/pathlength`.

### Display program usage

```bash
./target/release/pathlength -h
```

Outputs:

```bash
Calculates resolution and sensitivity in reflective superposition compound eyes

Usage: pathlength [OPTIONS]

Options:
  -f <FILENAME>      Path to a parameter file (CSV format). (Required)
  -d                 Generate debug CSV output file.
  -c                 Show the program citation.
  -l                 Show the program license.
  -v                 Show program version.
  -h, --help         Print help
  -V, --version      Print version
```

### Display citation information

```bash
./target/release/pathlength -c
```

Outputs:

```bash
Gaten, E., Moss, S., Johnson, M. 2013. The Reniform Reflecting Superposition Compound Eyes of Nephrops Norvegicus:
Optics, Susceptibility to Light-Induced Damage, Electrophysiology and a Ray Tracing Model. In: M. L. Johnson and M. P. Johnson, ed(s).
Advances in Marine Biology: The Ecology and Biology of Nephrops norvegicus. Oxford: Academic Press, 107:148.
```

### Display license information

```bash
./target/release/pathlength -l
```

### Display program version

```bash
./target/release/pathlength -v
```

Outputs:

```bash
pathlength version 0.6.0
```

## Run the program

```bash
./target/release/pathlength -f example_data/acanthephyra_parameters.txt
```

Outputs:

```bash
Parsing input parameters from "example_data/acanthephyra_parameters.txt"...
--- Running simulation for acanthephyra ---
Calculating pathlengths for acanthephyra...
P: 0.00, T: 0.00
P: 0.00, T: 12.70
...
INFO: Calculating resolution and sensitivity...
--- Finished simulation for acanthephyra ---
All simulations complete.
```

### Run with debug output

To generate optional `{species}_debug.csv` files containing simulation debug logs (e.g. pigment migration steps):

```bash
./target/release/pathlength -f example_data/acanthephyra_parameters.txt -d
```

## Required parameters

A CSV format file is required as input to the program. You can provide multiple lines for separate runs of the model. The format should be as follows:

```text
nephropsfl,180,25,7800,50,3200,1.34,1.37,18,0
nephropspl,180,25,7800,50,3200,1.34,1.37,18,12.5
nephropsfa,180,25,6760,50,3060,1.34,1.37,10,0
nephropspa,180,25,6760,50,3060,1.34,1.37,10,12.5
```

Each row is comprised of the following fields:

```text
genus	= A prefix for the output filenames e.g. organism genus name (lowercase alphanumeric only)
180 	= Rhabdom Length
25 	= Rhabdom Width
7800 	= Eye Diameter
50 	= Facet Width
3200	= Aperture Diameter
1.34	= Cytoplasm Refractive Index
1.37	= Rhabdom Refractive Index
18	= Blur Circle Extent
0	= Proximal Rhabdom Angle (used to create pointy-ended rhabdoms)
```

*NB: The genus name is NOT case sensitive. It is always converted to lowercase and should be unique to avoid filename conflicts.*

## Output files

The following output files are created:

* `genus_pathlengths.csv` - Raw pathlength data for each facet/pigment combination
* `genus_summary_res.csv` - Calculated resolution summary
* `genus_summary_sen.csv` - Calculated sensitivity summary
* `genus_debug.csv` - (Optional) Debug information from simulation, enabled with `-d`

The pathlengths file contains multiple rows for each facet with the various combinations of tapetal and shielding pigment lengths in the adjacent columns and then multiple columns with the pathlengths.

The summary resolution file contains the calculated resolution values.

The summary sensitivity file contains the calculated sensitivity values.

## Citation

If you use this program, please cite:

> Gaten, E., Moss, S., Johnson, M. 2013. The Reniform Reflecting Superposition Compound Eyes of Nephrops Norvegicus: Optics, Susceptibility to Light-Induced Damage, Electrophysiology and a Ray Tracing Model. In: M. L. Johnson and M. P. Johnson, ed(s). Advances in Marine Biology: The Ecology and Biology of Nephrops norvegicus. Oxford: Academic Press, 107:148.
