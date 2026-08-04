---
title: "Run comparison"
keywords:
  - model run comparison
  - compare runs
  - comparison index column
  - comparison split column
  - time series comparison
  - scatterplot comparison
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

Columns additionally require their tables to be comparable: the index columns must
name-match, and the split columns must either both be unset or name-match. Beyond that,
matching stays per column — the tables don't need to share their full column sets.

Matched values whose data is the same in every run (equal scalar values, or equal
value, index, and split column contents) are hidden from the compare list by default.
Numbers count as equal when they differ by less than 0.1%. Turn off
**Hide equal values** to show them.

## Fine-tuning matches

Each compare row maps the same value across runs: one item per run. To chart several
series within one run, use a [split column](#index-and-split-columns) instead — that is
the native pattern for model results.

Automatic matching picks the best candidate per run, which is not always the desired one —
for example, when a run has several similar tables. Expand a compare row (the chevron on
column targets) to see every compatible item grouped by run, with checkboxes:

* items picked by the automatic matching are checked and marked **(auto)**;
* other compatible items are listed unchecked — checking one makes it the run's series
  and unchecks the run's previous pick (one item per run);
* unchecking an item excludes the run from the chart;
* raw workspace tables are checked in every target whose name matches exactly (up to case
  and separators), not just the best one; fuzzy name matches are listed unchecked;
* **Reset to auto** discards the manual toggles for that target.

Toggles survive index and split changes, and the coverage count next to each target
reflects the checked items.

## Index and split columns

Column comparison requires an index (x axis) column per table. An optional split column
marks a categorical column whose values are separate series within one run.

* Numeric (int, bigint, float), datetime, and string columns are offered as indexes.
* Split columns must be string-typed and different from the index.
* With **Merge same functions** (on by default), tables produced by the same function
  collapse into one picker row; a selection applies to all runs.
* Charts are built by concatenating raw rows of all runs into one long dataframe; the line
  chart splits by run and, when set, by the split column. There is no row alignment or
  tolerance logic — differing grids simply chart as separate lines.

## Row semantics

When a table's index column is numeric or datetime, a selector appears next to its
pickers saying what the table's rows mean. The default is **series**: rows form an
indexed sequence, and values chart as lines connected along the index. A value charts in
a non-default mode only when **every** participating table has that mode set. Until then
it charts as a plain series, and the runs whose tables lack the setting are flagged with
an amber chip under Results.

### Relative timeseries

The **relative timeseries** mode charts the table on an elapsed-time axis, which lets
runs with different time representations — a datetime column in a workflow dataframe,
plain float seconds in a CSV — share one line chart.

* Numeric indexes get a **units** selector (ms, s, min, h, days) saying what the numbers
  mean. It is prefilled from the column's units metadata when recognizable. Datetime
  indexes need no units.
* Runs are aligned automatically: each table's earliest time point is treated
  as 0, so runs overlay for shape comparison regardless of when they were recorded.
* The chart axis is a float column labeled with its unit, for example `time (s)`. The
  unit comes from the first participating numeric index; a comparison of only datetime
  indexes picks a readable unit from the data range automatically.

### Independent points

The **independent points** mode treats rows as unordered observations and charts them as
unconnected (index, value) points instead of a line. Row order and grid alignment stop
mattering, so this mode suits unordered data — optimization populations, sampled or
stochastic results, per-item tables. And since any numeric column can serve as the index,
it also gives value-vs-value views: pick one output quantity as the index and compare it
against another across runs.

* Points are colored by the split column when one is set, otherwise by run. When the
  split takes the color, runs are distinguished by marker shapes.
* With [multiple values](#multiple-values) selected, all values chart on one scatterplot
  sharing the y axis, colored by value name with runs as marker shapes. Whether the
  values are comparable on one scale (units, magnitude) is up to you.

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

When the selected column target has compatible siblings (values sharing at least one run
from the same tables, with a line-chartable index), the **Multiple values** toggle (or
Shift+click on a compare row) enables selecting several values at once. Each value becomes
a stacked line-chart panel with runs (and split categories) as lines inside. In the
[independent points mode](#independent-points), the values share one scatterplot instead,
colored by value name.

Editing matches never exits the mode: a run that is missing or re-sourced from another
table in one of the selected values shows as a gap in that panel, and the value's row is
marked **partial**. Reverting the edit restores the chart exactly.
