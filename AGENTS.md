# AGENTS.md

This file provides guidance to AI coding agents when working with code in this repository.

## Project Overview

PathLength is a ray tracing model that calculates the resolution and sensitivity of reflective superposition compound eyes (found in crustaceans like Nephrops norvegicus). It's a Rust rewrite of an original QBASIC program from 1995, ported from a Python version.

## Development Commands

### Build
```bash
cargo build --release
```
This creates an optimized executable named `pathlength` in `target/release/`.

### Run from source
```bash
cargo run -- -f example_data/nephrops_parameters.txt
```
Note: The program requires a `-f` flag with a CSV parameter file path.

### Test
```bash
cargo test
```

### Lint and Format
```bash
cargo fmt --check
cargo clippy -- -D warnings
```

### Display help
```bash
cargo run -- -h
```

## Development Workflow

**IMPORTANT: Always run tests before committing changes.**

This project uses pre-commit hooks to automatically enforce code quality. The hooks will:
- Format Rust code with `cargo fmt`
- Run type checks with `cargo check`
- Run lint analysis with `cargo clippy`
- Run all tests with `cargo test`
- Check for common issues (merge conflicts, large files, etc.)

### Commit Message Format

Use conventional commit format for all commit messages:

```
<type>(<scope>): <description>

[optional body]

[optional footer]
```

**Types:**
- `feat`: New feature
- `fix`: Bug fix
- `refactor`: Code refactoring (no functional changes)
- `test`: Adding or updating tests
- `docs`: Documentation changes
- `chore`: Maintenance tasks (dependencies, build config, etc.)
- `perf`: Performance improvements

## Code Architecture

### Core Components

The codebase is organized as a library and CLI binary:

1. **`src/main.rs`** - Entry point with CLI argument parsing via `clap`
   - Handles flags: `-f` (file), `-d` (debug), `-h` (help), `-v` (version), `-c` (citation), `-l` (license)
   - Parses parameter file and orchestrates simulation runs

2. **`src/lib.rs`** - Module exports and `parse_input_parameters`

3. **`src/parameters.rs`** - `Parameters` struct holding 10 eye-specific configuration values

4. **`src/model.rs`** - Simulation engine and ray tracing logic
   - `Model` struct: Contains calculated parameters and simulation state
   - `Model::new()`: Computes derived values (eye radius, critical angle via Snell's law, facet counts, etc.)
   - `Model::run_simulation()`: Ray tracing loop iterating over 11×11 pigment steps and facets
     - Four main cases: perpendicular ray, no reflection, edge reflection, base bounce
     - Prepend blur circle shifts
     - Pointy-ended rhabdom acceptance angle offset

5. **`src/analysis.rs`** - Post-processes pathlengths to calculate resolution and sensitivity matrix summaries

6. **`tests/model_tests.rs`** - Test suite covering initialization, ray tracing, blur circles, pointy rhabdoms, and summary matrix output

## Output Files

For each parameter set with species name "X", generates:
- `X_pathlengths.csv` - Raw pathlength data for each facet/pigment combination
- `X_summary_res.csv` - Calculated resolution values
- `X_summary_sen.csv` - Calculated sensitivity values
- `X_debug.csv` - (Optional, generated with `-d`) Debug information from simulation

## Citation

When using this program, cite:
Gaten, E., Moss, S., Johnson, M. 2013. The Reniform Reflecting Superposition Compound Eyes of Nephrops Norvegicus: Optics, Susceptibility to Light-Induced Damage, Electrophysiology and a Ray Tracing Model. In: M. L. Johnson and M. P. Johnson, ed(s). Advances in Marine Biology: The Ecology and Biology of Nephrops norvegicus. Oxford: Academic Press, 107:148.
