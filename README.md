# ExpansionScan

v1.2.1

ExpansionScan is a C++ / ROOT utility to analyze in-plane thermal expansion
(or contraction) of a flat object measured with a Coordinate Measurement
Machine (CMM).

Given two point files corresponding to warm and cold measurements of
the same object, the program:

- computes point-by-point displacements in the XY plane
- fits a global 2D similarity transform (translation + rotation + scale)
- evaluates residuals after removing the global expansion
- produces diagnostic plots
- optionally writes a detailed CSV summary

--------------------------------------------------------------------------

Features

- 2D similarity fit in the XY plane (translation, rotation, uniform scale)
- Displacements computed BEFORE and AFTER the fit
- Histograms of dX, dY, and R (before and after)
- XY scatter plots with color-coded displacement magnitude
- Optional numeric labels at each point
- ROOT file output and PNG snapshots
- Optional CSV output with stable formatting

--------------------------------------------------------------------------

Input format

Input files are parsed using the shared Points.h / Points.cpp module.

Each line must begin with a non-numeric label, followed by coordinates.
At least three coordinates are required:

    <label>  X  Y  Z

Example:

    P01  123.456  78.901  2.345

- Coordinates are interpreted as millimeters
- Z is read but NOT used in the fit

The warm and cold files must contain the same points in the same order.

--------------------------------------------------------------------------

Usage

    ./ExpansionScan [options] warm.txt cold.txt [out.root]

Positional arguments (order is fixed):

- warm.txt   : warm measurement file
- cold.txt   : cold measurement file
- out.root   : optional ROOT output file
               default: ExpansionScan.root

Options (may appear anywhere on the command line):

- --labels
    Draw numeric displacement labels on scatter plots

- --dT <value>
    Temperature difference (Twarm − Tcold)
    Used only to print the thermal expansion coefficient

- --csv <file>
    Write per-point CSV summary to <file>

Unknown options are rejected with an error.

--------------------------------------------------------------------------

Output

ROOT and PNG output (always produced):

- ROOT file containing all histograms and canvases
- PNG snapshots of:
  - displacement histograms (before / after fit)
  - XY scatter plots with arrows and color scale

CSV output (optional):

When --csv <file> is specified, a CSV file is written with one row per point.

Columns:

    label,
    x_warm_mm, y_warm_mm, z_warm_mm,
    x_cold_mm, y_cold_mm, z_cold_mm,
    dx_before_um, dy_before_um, r_before_um,
    dx_after_um,  dy_after_um,  r_after_um

Units and formatting:

- Coordinates are in millimeters, written with 3 decimal digits
- Displacements are in micrometers, written as integers
- CSV output is purely additive and does not affect plots or ROOT output

--------------------------------------------------------------------------

Definitions

For each point i:

Warm coordinates:
    (xw, yw)

Cold coordinates:
    (xc, yc)

Before fit:

    dX = xc − xw
    dY = yc − yw
    R  = sqrt(dX^2 + dY^2)

After fit:

Warm coordinates are first transformed using the best-fit similarity transform:

    (xw, yw) → (xw’, yw’)

Residuals are then computed:

    dX’ = xc − xw’
    dY’ = yc − yw’
    R’  = sqrt(dX’^2 + dY’^2)

--------------------------------------------------------------------------

Build

Requirements:

- C++17 compiler
- CERN ROOT (tested with Homebrew ROOT on macOS)

Typical build:

    make clean
    make

--------------------------------------------------------------------------

Notes

- Only the first min(Nwarm, Ncold) points are used if file lengths differ
- All internal computations use millimeters
- Arrow lengths in scatter plots are normalized for readability
- The interactive ROOT session starts after all computations and outputs

--------------------------------------------------------------------------

Versioning

v1.2 introduces:

- CSV output with before/after fit displacements
- Stable numeric formatting
- Updated documentation and usage

--------------------------------------------------------------------------

Author

Luciano Ristori
