# Project Context: PathLength

## Overview
PathLength is a scientific simulation tool written in Go (Golang) that models the optics of **reflective superposition compound eyes**. It uses ray tracing to calculate two key optical properties:
1.  **Resolution**
2.  **Sensitivity**

This project is a modern rewrite (v0.6.0) of an original QBASIC program from 1995. It is designed to take biological parameters of an eye (e.g., *Nephrops norvegicus*) and output detailed optical performance metrics.

## Architecture
The application is a command-line interface (CLI) tool structured as follows:

*   **Entry Point (`pathlength.go`):** Handles command-line flags (`-f`, `-c`, `-l`, `-v`), prints help/license info, and orchestrates the simulation loop for each parameter set found in the input file.
*   **Data Handling (`csv.go`):**
    *   **Input:** Parses the input CSV file containing biological parameters (genus, rhabdom dimensions, refractive indices, etc.).
    *   **Output:** Implements `calculateRessens` to process raw pathlength data and generate summary CSV files for resolution and sensitivity.
*   **Core Logic (`model.go`):** Contains the `Model` struct and the `runModel()` method which performs the actual ray tracing computations (calculating light paths through the eye's geometry).

## Usage & Development

### Prerequisites
*   **Go:** Version 1.25 or later.
*   **Dependencies:** `github.com/stretchr/testify` (for testing).

### Building and Running
1.  **Run directly:**
    ```bash
    go run . -f example_data/acanthephyra_parameters.txt
    ```
2.  **Build binary:**
    ```bash
    go build
    ./pathlength -f example_data/acanthephyra_parameters.txt
    ```

### Testing
The project uses the standard Go testing framework with `testify` assertions.
```bash
go test -v
```

### Input Data Format
The input must be a CSV file with 10 fields per row:
1.  **Genus:** Species name (used for output filenames).
2.  **Rhabdom Length:** (microns)
3.  **Rhabdom Width:** (microns)
4.  **Eye Diameter:** (microns)
5.  **Facet Width:** (microns)
6.  **Aperture Diameter:** (microns)
7.  **Cytoplasm Refractive Index:** (float)
8.  **Rhabdom Refractive Index:** (float)
9.  **Blur Circle Extent:** (integer)
10. **Proximal Rhabdom Angle:** (degrees)

### Output Files
For each species defined in the input, the program generates:
*   `{genus}_debug.csv`: Debugging information.
*   `{genus}_pathlengths.csv`: Raw ray tracing path length data.
*   `{genus}_summary_res.csv`: Calculated resolution summary.
*   `{genus}_summary_sen.csv`: Calculated sensitivity summary.
