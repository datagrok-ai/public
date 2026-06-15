---
title: "Data types"
sidebar_position: 6
---

When a script runs, Datagrok converts each input parameter to a native value in the
target language and converts the script's output back to a Datagrok value. The table
below shows the equivalent type used for each Datagrok parameter type, so you can
write the script body knowing exactly what you receive and what you must return.

## Type correspondence across languages

| Datagrok       | Python                | R                    | Octave                  | Julia                    | JavaScript            |
|----------------|-----------------------|----------------------|-------------------------|--------------------------|-----------------------|
| `int`          | `int`                 | `integer`            | `int64` (cast)          | `Int64`                  | `number`              |
| `double`       | `float`               | `numeric`            | `double`                | `Float64`                | `number`              |
| `bool`         | `bool`                | `logical`            | `logical`               | `Bool`                   | `boolean`             |
| `string`       | `str`                 | `character`          | `char`                  | `String`                 | `string`              |
| `datetime`     | `datetime.datetime`   | `POSIXct`            | `datenum` (double)      | `Dates.DateTime`         | `dayjs` instance      |
| `dataframe`    | `pandas.DataFrame`    | `data.frame`         | cell array              | `DataFrames.DataFrame`   | `DG.DataFrame`        |
| `column`       | `str` (column name)   | `character` (name)   | `char` (name)           | `Symbol`                 | `string` (name)       |
| `column_list`  | `list[str]`           | `character` vector   | cell of strings         | `Vector{Symbol}`         | `string[]`            |
| `list<string>` | `list[str]`           | `character` vector   | cell array              | `Vector{String}`         | `string[]`            |
| `map`          | `dict`                | named `list`         | `struct`                | `Dict`                   | object                |
| `file`         | path `str`            | path `character`     | path `char`             | path `String`            | `FileInfo` (`.data`)  |
| `blob`         | `bytes`               | `raw`                | `uint8` array           | `Vector{UInt8}`          | `Uint8Array`          |
| `graphics`     | `matplotlib.Figure`   | base `plot()`        | any plot                | `Plots.jl` + `display()` | `DG.Viewer` / element |

Notes:

- Python, R, Octave, and Julia scripts run server-side in the
  [Jupyter Kernel Gateway](../advanced-scripting/under-the-hood.mdx). DataFrames are
  serialised as CSV (or Parquet when enabled) over HTTP. Graphics outputs are
  captured from the kernel's display channel (PNG or SVG).
- JavaScript scripts run in the browser. `DG.DataFrame`, `dayjs`, and `FileInfo` are
  passed by reference without serialisation. A `graphics` output is any element or
  viewer returned by the script.
- `column` is passed as the column **name**, not the column object — look it up
  inside the script (for example, `df[xName]` in Python, `df[!, xName]` in Julia,
  or `df.col(xName)` in JavaScript).
- `list` requires an explicit element type, currently `list<string>` only. For
  numeric or mixed lists, pass a single-column `dataframe`.

## See also

- [Complex data types](complex-input-output.md) — dataframe, file, blob, and graphics usage examples
- [Function annotations](../../../datagrok/concepts/functions/func-params-annotation.md) — how to declare parameter types
