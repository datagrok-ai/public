---
title: "Run comparison"
keywords:
  - model run comparison
  - compare runs
  - comparison index column
  - comparison split column
---

The run comparison tool compares outputs of model runs side by side: pick a model, add
historical runs (or raw workspace tables) to the comparison set, choose what to compare, and
get a chart preview. The tool is a preliminary data selection with a quick preview — complex
analysis belongs in the Datagrok workspace, where any comparison can be exported as a snapshot
(**Results > Open in workspace**: the long-format data table plus the chart with its current
options).

## Matching

Scalars and dataframe columns are matched across runs by name: exact, then normalized
(case/underscore/dash-insensitive), then fuzzy. Values with mismatching units never match;
partially missing units produce a warning icon.

## Index and split columns

Column comparison requires an index (x axis) column per table. An optional split column
marks a categorical column whose values are separate series within one run.

* By default only **int, bigint, and string** columns are offered as indexes; float and
  datetime indexes are enabled with the **Float indexes** / **Datetime indexes** toggles.
* Split columns must be string-typed and different from the index.
* With **Merge same functions** (on by default), tables produced by the same function
  collapse into one picker row; a selection applies to all runs.
* Charts are built by concatenating raw rows of all runs into one long dataframe; the line
  chart splits by run and, when set, by the split column. There is no row alignment or
  tolerance logic — differing grids simply chart as separate lines.

## Default index/split annotations

A model can declare the default index and split columns on its dataframe output:

```javascript
//output: dataframe result {comparisonIndex: time; comparisonSplit: species}
```

When a run is added to the comparison:

* `comparisonIndex` pre-fills the index picker, and the column is offered even when its
  type is behind a toggle (an explicit annotation outranks the type gating);
* `comparisonSplit` pre-fills the split picker, provided the column is a valid split
  (string-typed, not the index);
* annotations never overwrite a selection the user has already made, and a missing column
  name is silently ignored;
* merged rows work naturally: all merged tables come from the same function, so they carry
  the same annotation and agree on the default.

Raw workspace tables have no annotations.

## Multiple values

When the selected column target has compatible siblings (same runs, tables, index and split
columns, and a line-chartable index), the **Multiple values** toggle (or Shift+click on a
compare row) enables selecting several values at once. Each value becomes a stacked
line-chart panel with runs (and split categories) as lines inside.
