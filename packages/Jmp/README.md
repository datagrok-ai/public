# JMP Reader

**JMP Reader** opens [JMP](https://www.jmp.com/)® data tables (`.jmp`) in Datagrok. The file is
parsed in the browser, so you don't need a JMP installation or any server component.

JMP® is a registered trademark of JMP Statistical Discovery LLC, a subsidiary of SAS Institute Inc.
This plugin is not affiliated with, sponsored by or endorsed by either company.

## Supported File Types

| Extension | Name           | Support          |
| --------- | -------------- | ---------------- |
| `.jmp`    | JMP data table | Import, preview  |

## Usage

- **From disk:** drop a `.jmp` file onto Datagrok, or use File > Open > File... It opens as a
  table, named after the file, with its saved analyses as viewers (see below).
- **From the Browse tree:** click a `.jmp` file to preview it as its table with its saved
  analyses, as a CSV file previews; double-click to keep it open. The **About this JMP file** icon
  in the ribbon shows the JMP version that saved it, its size, when it was saved, the columns left
  out and each saved analysis with what it opens as or why it doesn't.

If some columns can't be read, the rest of the table still opens. A warning names the columns left
out, and the summary lists them under **Columns not imported**.
A file that isn't a readable JMP table says why instead of opening.

## What gets imported

| JMP                                   | Datagrok                         |
| ------------------------------------- | -------------------------------- |
| Numeric                               | `double`                         |
| 1/2/4-byte integer                    | `int`                            |
| Character                             | `string`                         |
| Date, date-time                       | `datetime`                       |
| Duration                              | `double` (seconds, units `s`)    |
| Latitude, longitude                   | `double` with the GIS semantic types |
| Currency                              | `double`, with the currency code as its units |
| Row state column                      | two columns, `<name> marker` and `<name> color`, drawn with the JMP markers and colors |
| Missing values                        | nulls                            |
| Column with value labels              | `string` of the labels, as JMP shows them; the values are kept in `jmp.valueLabels` |

Column properties come across where Datagrok has a counterpart: **Notes** become the column's
description, **Units** its units, and a custom **Value Order** its category order. A **Formula** is
kept as the `jmp.formula` tag (the computed values are imported). **Spec Limits** are kept as
`jmp.specLimits` and drawn as LSL, USL and Target lines on every chart that plots the column (hidden
when JMP doesn't show them). Every property is also kept as JSL text in `jmp.properties`.

The JMP version, build and save time are kept as table tags (`jmp.version`, `jmp.buildString`,
`jmp.saveTime`). Columns hidden in JMP are tagged `jmp.hidden` and hidden in the grid.

Not imported yet: which rows are excluded, hidden or selected; value colors; formulas as
calculated columns.

## Saved analyses

JMP keeps a table's analyses as saved scripts. Those the package understands open as Datagrok
viewers, titled after the script:

| JMP                       | Datagrok                                   |
| ------------------------- | ------------------------------------------ |
| Distribution              | Histogram, or a bar chart of counts         |
| Bivariate, Fit Y by X     | Scatter plot, with the regression line for a line fit |
| Oneway                    | Box plot by group                          |
| Contingency               | Stacked bar chart of counts                |
| Graph Builder             | Scatter plot, line chart, bar chart, box plot, histogram or pie chart, one per panel, by its first element, with its title, axis ranges, log scales and reference lines; a trellis plot of it when the graph is grouped or paged |
| Tabulate                  | Pivot table                                |
| Multivariate              | Correlation plot (Pearson, or Spearman when only that was asked for) and a matrix plot |
| Scatterplot Matrix        | Matrix plot of scatter plots, coloured by its Group |
| Matched Pairs             | Scatter plot of one Y against its pair, on the line of no difference; a matrix plot for several Ys |
| Chart                     | Bar chart (or pie chart) of each Y's statistic by the X categories |
| Fit Curve                 | Scatter plot of the data it fits, coloured by its Group; a straight-line fit is drawn, the nonlinear ones are named |
| DOE: Evaluate Design      | Correlation plot of the design's continuous factors |
| Fit Model, Fit Group (least squares) | Actual by Predicted scatter plot per response, from the prediction formula saved from that model |
| Fit Group of other analyses | Each analysis it holds, as above |

Anything else (other Fit Model personalities, the Profiler, control charts, DOE, ...) is named in a
warning when the file opens, and listed in the file summary with the reason. What a converted viewer
leaves out, such as a `By()` split, is listed under the analysis in the summary, and **Script** there
shows the script itself.

Where JMP keeps scripts in the file was worked out from a few files, so a table whose scripts aren't
found opens without them; its data is unaffected.

## Limitations

The `.jmp` format is proprietary and undocumented. This reader is a port of the independent
[JMPReader.jl](https://github.com/jaakkor2/JMPReader.jl). Files saved in an older JMP format are
rejected with a message asking you to re-save them in a current JMP version.

## Development

The design, the format as far as it is known, and the roadmap are in
[docs/JMP_PLUGIN.md](docs/JMP_PLUGIN.md).

```bash
grok test --host <stand>       # package tests: the parser in the Node pass, the platform parts in the browser
npm run corpus -- <folder>     # read every .jmp file in a folder and report (for files that can't be committed):
                               # columns read, saved analyses and how they convert, row-state flags seen
npm run large -- --mb 100      # write a synthetic 100 MB file and time reading it
grok test --host <stand> --benchmark --category "JMP: performance"   # the same in the browser
```

Behaviour is specified as Gherkin features in `bdd/features/` (see `libraries/bdd`):

```bash
node ../../libraries/bdd/bin/grok-bdd.js compile
DATAGROK_URL=<stand url> node ../../libraries/bdd/bin/grok-bdd.js run --reporter=list
```

## Test files

`files/` is published with the package, so its files also show in Browse under
Files › App Data › JmpReader.

| Folder | What it holds | Used by |
|---|---|---|
| `jmpreader/` | JMPReader.jl's 19 test files (MIT), unchanged | parser tests |
| `golden/` | what JMPReader.jl reads from each of them, as JSON; `tools/jmp2json.jl` regenerates it | reference tests, which compare our parse with it |
| `tests/` | files made here: a truncated file, CSV under a `.jmp` name, a column layout no JMP writes, synthetic saved analyses; one JMP 15 file shared by a user and anonymized | robustness and analysis tests |
| `public/` | 12 research tables published openly (CC0, CC BY 4.0), saved by JMP 10–18 | tests of real files, BDD features, the how-to guide |

Where each file comes from and its licence: [CREDITS.md](CREDITS.md) and
[docs/JMP_PLUGIN.md](docs/JMP_PLUGIN.md) §9.2.

See [CREDITS.md](CREDITS.md) for third-party attributions.
