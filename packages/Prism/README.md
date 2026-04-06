# Prism

**Prism** provides support for [GraphPad Prism](https://www.graphpad.com/) `.prism` files
(version 10+), including file preview with navigation, data import, and dose-response
curve visualization via the [Curves](../Curves) package.

## Features

- **File preview**: tree-based navigation across data sheets, analyses, and graphs
- **Data import**: all data sheets are converted to Datagrok DataFrames
- **All data formats**: y_replicates, y_single, y_sd, y_se, y_cv, y_sd_n, y_se_n, y_cv_n,
  y_plus_minus, y_high_low
- **Curve fitting**: XY dose-response data is automatically converted to fitted curves,
  rendered by the Curves package
- **Analysis results**: nonlinear fits, survival analyses, and other Prism analysis
  results are extracted and displayed as tables
- **Automatic visualization**: data sheets get appropriate chart viewers (scatter plots,
  bar charts) based on their type

## File format

The `.prism` format (introduced in GraphPad Prism 10) is an open ZIP archive containing:

- `document.json` — top-level manifest linking sheets, analyses, and graphs
- `data/sheets/` — per-sheet metadata (format, column definitions, replicate count)
- `data/tables/` — raw data as CSV files
- `data/sets/` — dataset metadata (column titles)
- `analyses/` — analysis parameters and result references

All data, parameters, and results are stored in industry-standard formats (CSV, JSON),
making them accessible without GraphPad Prism.

## Usage

### Preview

Single-click a `.prism` file in the file browser to see a tree-based preview:

- **Data Sheets** — browse all tables with their format labels
- **Analyses** — view curve fit results, statistics, and other analysis outputs

Clicking a sheet in the tree shows the corresponding grid and chart.

### Import

Double-click a `.prism` file (or drag-and-drop) to import all data as DataFrames:

- Each data sheet becomes a separate DataFrame
- XY dose-response sheets also produce a "Prism Curves" DataFrame with a `fit`
  column, automatically rendered by the Curves package
- Analysis results are imported as additional DataFrames

### Curves integration

When XY dose-response data is detected, the plugin creates cells with the `fit` semantic
type. If the [Curves](../Curves) package is installed, these cells are rendered as
interactive fitted curves with:

- Sigmoid (4PL) curve fitting
- IC50 droplines
- Click-to-toggle outlier support
- Multi-curve overlay

## Data format mapping

| Prism Format   | DataFrame Columns                          |
|----------------|--------------------------------------------|
| y_single       | 1 value column per dataset                 |
| y_replicates   | N replicate columns (Title_1, Title_2, ...) |
| y_sd           | Mean, SD                                   |
| y_se           | Mean, SEM                                  |
| y_cv           | Mean, %CV                                  |
| y_sd_n         | Mean, SD, N                                |
| y_se_n         | Mean, SEM, N                               |
| y_cv_n         | Mean, %CV, N                               |
| y_plus_minus   | Mean, +Error, -Error                       |
| y_high_low     | Mean, Upper Limit, Lower Limit             |

## Demo files

The package includes two sample files in `files/`:

- **demo_dataset.prism** — comprehensive file with 14 data sheets covering all
  formats, 2 analyses (nonlinear fit + survival), and 17 graphs
- **no_analysis.prism** — simple file with 2 data sheets and no analyses

See also:

- [Curves](../Curves)
- [Packages](../../help/develop/develop.md#packages)
- [JavaScript API](../../help/develop/packages/js-api.md)
