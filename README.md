/* ==============================================================================
File:        ExpansionScan.cpp
Version:     1.3.0
Author:      Luciano Ristori

Description:
------------
ExpansionScan analyzes thermal expansion effects by comparing two sets of
measured 3D points (e.g. warm vs cold scans) acquired with a CMM.

For each matched point, the program:
  • Computes raw displacements (warm → cold)
  • Fits a global expansion model
  • Computes residual displacements after the fit
  • Produces ROOT plots (histograms and XY scatter plots)
  • Optionally writes a detailed CSV summary

The program is designed to be run from a dataset-specific *results directory*.

Input and output path handling:
-------------------------------
• Input files:
    - Absolute paths are used as given
    - Relative paths are interpreted relative to $CMM_DATA if set,
      otherwise relative to the current working directory

• Output files (ROOT, PNG, CSV):
    - Absolute paths are used as given
    - Relative paths are interpreted relative to $CMM_RESULTS/ExpansionScan
      if CMM_RESULTS is set
    - Otherwise, relative paths refer to the current working directory

Command-line usage:
-------------------
ExpansionScan <warm.csv> <cold.csv> [output.root] [options]

Positional arguments:
  warm.csv        CSV file with warm measurements
  cold.csv        CSV file with cold measurements
  output.root     Optional ROOT output file
                  (default: ExpansionScan.root)

Options:
  --labels            Draw point labels on XY scatter plots
  --dT <value>        Specify temperature difference (ΔT)
  --csv <file.csv>    Write a CSV summary of per-point results

CSV output:
-----------
When --csv is specified, a CSV file is written with one row per point,
containing:

  label,
  warm_x_mm, warm_y_mm, warm_z_mm,
  cold_x_mm, cold_y_mm, cold_z_mm,
  dx_before_um, dy_before_um, dr_before_um,
  dx_after_um,  dy_after_um,  dr_after_um

Formatting conventions:
  • Coordinates in millimeters are written with 3 decimal digits
  • Displacements in micrometers are written as integers

Notes:
------
• Plotting behavior and --labels semantics are unchanged from previous versions
• CSV generation does not alter any analysis or plotting logic
• Executables and build artifacts are not tracked in git

=============================================================================== */


/* ============================== README.md ===================================

# ExpansionScan

**ExpansionScan** is a C++/ROOT tool for analyzing thermal expansion effects
by comparing two sets of 3D point measurements (e.g. warm vs cold CMM scans).

The program computes raw displacements, fits a global expansion model, evaluates
residuals after the fit, and produces both graphical output and optional CSV
summaries.

## Features

- Reads two CSV point files with identical point ordering
- Computes warm → cold displacements
- Fits a global expansion model
- Computes residual displacements after the fit
- Produces ROOT histograms and XY scatter plots
- Optional point labeling on scatter plots
- Optional CSV output with per-point results

## Build

Requires:
- C++17
- ROOT (built with C++17 support)

Build using:
    make

## Usage

ExpansionScan <warm.csv> <cold.csv> [output.root] [options]

Options:
  --labels            Draw point labels on XY scatter plots
  --dT <value>        Specify temperature difference ΔT
  --csv <file.csv>    Write a detailed CSV summary

## CSV Output Format

Columns:
  label,
  warm_x_mm, warm_y_mm, warm_z_mm,
  cold_x_mm, cold_y_mm, cold_z_mm,
  dx_before_um, dy_before_um, dr_before_um,
  dx_after_um,  dy_after_um,  dr_after_um

Formatting:
- mm values: 3 decimal digits
- µm values: integers

## Directory Model

Recommended layout:

  ~/code/CMM/ExpansionScan
  ~/Dropbox/CMM/data
  ~/Dropbox/CMM/results/ExpansionScan/<run>

Environment variables:
  CMM_DATA     root directory for input data
  CMM_RESULTS  root directory for output results

## Versioning

v1.3.0
- Added CSV output with before/after fit displacements
- Defined input/output directory model
- Added support for CMM_DATA and CMM_RESULTS
- No changes to plotting or analysis logic

=============================================================================== */
