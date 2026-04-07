# Minitab

**Minitab** provides support for [Minitab](https://www.minitab.com/) file formats,
including file preview, data import, and workspace integration.

## Supported File Types

| Extension | Name                | Support                              |
| --------- | ------------------- | ------------------------------------ |
| `.mwx`    | Minitab Worksheet   | Import, preview                      |
| `.mpx`    | Minitab Project     | Import (all worksheets), preview     |

## Features

- **MWX import**: single worksheet converted to a Datagrok DataFrame
- **MWX preview**: grid view with worksheet name and dimensions
- **MPX import**: all worksheets in the project are imported as separate DataFrames
- **MPX preview**: project summary (creator, comments, worksheet list, session history)
  with an "Open" button to import all worksheets into the workspace

## Usage

### Preview

Single-click an `.mwx` or `.mpx` file in the file browser to see a preview:

- **MWX** — shows the worksheet data in a grid
- **MPX** — shows project metadata, a table of worksheets with row/column counts,
  and collapsible session history

### Import

Double-click a file (or drag-and-drop) to import:

- **MWX** — opens one DataFrame
- **MPX** — opens one DataFrame per worksheet in the project

## File Format Details

Both `.mwx` and `.mpx` are ZIP archives containing JSON:

- **MWX**: `sheet_metadata_20.json` + `sheets/0/sheet.json`
- **MPX**: `project_metadata_20.json` + `sheets/N/sheet.json` (one per worksheet)

Column types are mapped as follows:

| Minitab Type | Datagrok Type |
| ------------ | ------------- |
| Numeric      | float64       |
| Text         | string        |
| Date/Time    | datetime      |

Missing values (`*` in Minitab) are converted to nulls.

## Demo Files

The package includes sample files in `files/`:

- **ChocolatePreferences.MWX** — 400 rows, 2 text columns
- **PaintHardness.MWX** — 24 rows, 1 text + 3 numeric columns
- **SoftwareSales.MWX** — 26 rows, 3 numeric columns
- **ANOVA-SALES.mpx** — 2 worksheets (stacked and unstacked sales data)
- **SPC-I-MR.mpx** — individuals and moving range chart data
- **SPC-XBAR-R.mpx** — 2 worksheets for Xbar-R chart
- **SPC - P CHART.mpx** — proportion defective data
- **SPC - U CHART.mpx** — defects per unit data
- **TEST FOR EQUAL VARIANCES.mpx** — variance test data

See also:

- [Packages](../../help/develop/develop.md#packages)
- [JavaScript API](../../help/develop/packages/js-api.md)
