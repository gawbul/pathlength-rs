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

3. **`src/parameters.rs`** - `Parameters` struct holding 10 eye-specific configuration
   values, plus `validate()`, which rejects any eye that is not physically realisable
   rather than letting NaNs propagate into the results

4. **`src/model.rs`** - Simulation engine and ray tracing logic
   - `Model::new()`: Validates the parameters and computes the derived optics
     (ommatidial angle, facet count, critical angle via Snell's law)
   - `refracted_angle()`: Empirical corneal refraction regression
   - `Model::facet_transmission()`: Light a facet admits at a given incidence, as a
     flux factor in [0, 1]. Applied to absorbed intensity, never to path length
   - `Model::blur_offset()`: Continuous displacement in rhabdoms of the image formed
     by a facet, spanning 0 to (blur_circle_extent - 1) across the eyeshine patch
   - `Model::trace_ray()`: Follows one ray through the rhabdom array. Four cases:
     axial ray, no reflection, wall reflection, base bounce. Returns raw geometry in
     µm, the terminating case, and the largest angle reached
   - `Model::run_simulation()`: Iterates the 121 pigment states and writes the output.
     The absorption profile is accumulated as the rays are traced rather than by
     reading the file back, so the summary does not depend on the output format

5. **`src/analysis.rs`** - Summary
   - `ring_area()`: Area of the annulus of rhabdoms at a given offset. Used both to
     weight contributing facets and to normalise the light arriving at an offset
   - `deposit()`: Accumulates into a growable profile, so no rhabdom is silently dropped
   - `summarise_block()`: Turns one block's area-weighted profile into an acceptance
     angle (FWHM, degrees) and a sensitivity (percent absorbed)
   - `accumulate()`: Adds one facet's traced ray into the area-weighted absorption
     profile
   - `calculate_ressens()`: Writes the summary matrices from the per-state summaries
     accumulated during the simulation

6. **`tests/model_tests.rs`** - Test suite covering parameter validation, the blur
   circle mapping, ray geometry bounds, raw-geometry path lengths, guided rays against
   the screening pigment, profile accumulation, half-maximum interpolation, and the
   summary matrices

## Output Files

For each parameter set with species name "X", generates:
- `X_pathlengths.csv` - Raw ray geometry, in µm: a rectangular CSV with the header
  `block,shielding_um,tapetal_um,facet,rhabdom,pathlength_um` and one row per rhabdom
  entered. No block terminator and no positional state
- `X_summary_res.csv` - Acceptance angle (FWHM of the point spread function), degrees
- `X_summary_sen.csv` - Incident light absorbed, percent (0-100)
- `X_debug.csv` - (Optional, generated with `-d`) One row per traced ray

Both summary files are 11x11: **rows vary the shielding (proximal screening) pigment,
columns vary the tapetal (reflecting) pigment**. A resolution cell reading `NaN` means
the profile never falls to half its maximum, so the acceptance angle is undefined.

## Key Implementation Notes

- Ray angles (`boa`) are measured from the rhabdom axis, so the angle at the wall
  normal is (90 - boa) and light is guided while `boa < critical_angle`, where
  `critical_angle = 90 - asin(n_cytoplasm / n_rhabdom)`
- The blur circle displaces each facet's image by a **continuous** offset. Light is
  split between the two whole-rhabdom offsets that bracket it. Quantising this aliases
  the facets unevenly and cuts notches into the profile that the half-maximum search
  then locks onto
- The blur circle extent may not exceed the facet count across the eyeshine patch
- Rays reaching 90 degrees or more to the rhabdom axis are discarded with a warning,
  rather than being folded back by taking absolute values of `tan` and `cos`
- The absorption coefficient is fixed at 0.01 µm⁻¹ and tapetal reflectance at 1.0

### Known limitation

A trace terminates at the first reflection rather than following the reflected ray
onward. Because an unreflected ray continues leaking into adjacent rhabdoms and
absorbing there, extending the tapetum can occasionally shorten the total absorbing
path and lower the reported sensitivity, by up to about 6 percentage points. A tapetal
mirror can only add path length in reality, so this is a limitation of the 1995 case
structure rather than a property of the eye.

## Citation

When using this program, cite:
Gaten, E., Moss, S., Johnson, M. 2013. The Reniform Reflecting Superposition Compound Eyes of Nephrops Norvegicus: Optics, Susceptibility to Light-Induced Damage, Electrophysiology and a Ray Tracing Model. In: M. L. Johnson and M. P. Johnson, ed(s). Advances in Marine Biology: The Ecology and Biology of Nephrops norvegicus. Oxford: Academic Press, 107:148.
