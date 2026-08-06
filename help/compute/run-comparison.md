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
analysis belongs in the Datagrok workspace, where any comparison can be
[exported](#exporting) as a snapshot: the data table plus the chart with its current
options.

## Matching

Scalars and dataframe columns are matched across runs by name: exact, then normalized
(case/underscore/dash-insensitive), then fuzzy. Values with mismatching units never match;
partially missing units produce a warning icon.

Columns additionally require their tables to be comparable: the index columns must
name-match, and the split columns must either both be unset or name-match. Beyond that,
matching stays per column — the tables don't need to share their full column sets.

Matched values whose data is the same in every run (equal scalar values, or equal
value, index, and split column contents) are hidden from the compare list by default.
Numbers count as equal when they differ by less than 0.1%, and datetime index
values only when they match exactly. Turn off **Hide equal values** to show them.

## Fine-tuning matches

Each compare row maps the same value across runs. Several items of one run can match
the row, but only one of them is active at a time: the active item is the run's charted
series. To chart several series within one run, use a
[split column](#index-and-split-columns) instead — that is the native pattern for
model results.

When an input and an output of the same step share a caption, the parameter name is
shown after it, for example `Data (result)`, so the two rows stay distinguishable.

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

## Default index annotations

A model can declare how its dataframe output should be indexed in a comparison with the
`comparison` option — a small JSON object on the output:

```javascript
//output: dataframe result {comparison: {"index": "time", "split": "species", "mode": "timeseries", "units": "s"}}
```

| Key | Value | Pre-fills |
|-----|-------|-----------|
| `index` | column name | the index (x axis) picker |
| `split` | column name | the split picker |
| `mode` | `series`, `timeseries`, or `points` | the [row semantics](#row-semantics) selector |
| `units` | `ms`, `s`, `min`, `h`, or `days` | the time units selector |

All keys are optional. When a run is added to the comparison:

* `index` pre-fills the index picker, and the column is offered even when its
  type is behind a toggle (an explicit annotation outranks the type gating);
* `split` pre-fills the split picker, provided the column is a valid split
  (string-typed, not the index);
* `mode` pre-fills the row semantics selector. `units` applies only together with
  `mode: timeseries` on a numeric index and outranks the column's units metadata;
* annotations never overwrite a selection the user has already made; a missing column
  name or an unknown mode or units value is silently ignored;
* merged rows work naturally: all merged tables come from the same function, so they carry
  the same annotation and agree on the default.

Raw workspace tables have no annotations.

## Multiple values

When the selected value has compatible siblings, the **Multiple values** toggle (or
Shift+click on a compare row) enables selecting several values at once. Unchecking the
last selected value turns the mode off, keeping that value as the regular selection. Every scalar is
compatible with every other scalar. Column values are compatible when they share at least
one run from the same tables and their index is numeric or datetime everywhere. Scalars
and columns never mix in one selection.

Several columns chart as stacked line-chart panels with runs (and split categories) as
lines inside. In the [independent points mode](#independent-points), the values share one
scatterplot instead, colored by value name. A comparison of several values needs a numeric
or datetime index — if an index moves to a string column, the chart is replaced with a
hint until the index is changed back or a single value is selected.

Several scalars chart on a radar chart: one axis per value, one polygon per run. A switch
above the chart flips it to a parallel-coordinates plot. With two values the tool always
uses the parallel-coordinates plot, because a two-axis radar degenerates to a line. The
radar viewer comes with the Charts package — without it, the parallel-coordinates plot is
used for any number of values.

Editing matches never exits the mode: a run that is missing or re-sourced from another
table in one of the selected values shows as a gap in that panel, and the value's row is
marked **partial**. Reverting the edit restores the chart exactly.

## Exporting

**Results > Export...** asks for a snapshot name and offers two ways out:

* **Open in workspace** adds the comparison to the workspace: the data table plus the
  chart with its current options. Tweak the chart before exporting — the tweaks carry
  over.
* **Save & share** opens the platform project dialog on the same snapshot, so the
  comparison can be saved as a project and shared without leaving the tool. On platforms
  without the project dialog, only **Open in workspace** is offered.
