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
- **From the Browse tree:** click or double-click a `.jmp` file to preview it: the JMP version
  that saved it, its size, when it was saved, its saved analyses, and the data. **Open** in the
  ribbon opens the table with its analyses.

If some columns can't be read, the rest of the table still opens. The preview lists the columns
that were left out under **Columns not imported**, and opening from disk names them in a warning.
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
kept as the `jmp.formula` tag (the computed values are imported) and **Spec Limits** as
`jmp.specLimits`. Every property is also kept as JSL text in `jmp.properties`.

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
| Graph Builder             | Scatter plot, line chart, bar chart, box plot, histogram or pie chart, by its first element; a trellis plot of it when the graph has panels |
| Tabulate                  | Pivot table                                |

Anything else (Fit Model, control charts, DOE, ...) is listed in the preview with the reason, and
named in a warning when the file opens. What a converted viewer leaves out, such as a `By()` split, is
listed under the analysis in the preview, and **Script** there shows the script itself.

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

`files/golden/` holds what JMPReader.jl reads from each test file; `tools/jmp2json.jl` regenerates
it.

See [CREDITS.md](CREDITS.md) for third-party attributions.
