---
feature: curves
target_layer: playwright
pyramid_layer: ui-smoke
coverage_type: regression
priority: p0
related_bugs:
  - GROK-13061
  - GROK-16166
  - GROK-17589
  - GROK-20662
  - GROK-20664
  - public#2100
  - public#2101
---

# Curve cells in the grid

The Curves package detects columns that hold fitted-curve data (semantic type `fit`) and draws a
chart in every cell. Clicking a curve cell shows three panes in the Context Panel: **Options**
(chart and series options, with the **Level** they are set at), **Chart** (the chart, enlarged) and
**Fit** (the fit statistics of each series). Double-clicking a cell opens the chart in the
**Edit chart** dialog.

## Setup

A clean session. Open `System:DemoFiles/curves.csv` (30 rows) from **Browse > Files > Demo**. It has
five curve columns: `multiple prefit` and `multiple styled series` (two series in a cell),
`styled proprtional with IC50` (one series), `log min/max candlesticks`, `custom fit function`.

## Scenarios

### 1. Curve columns are detected and drawn (GROK-13061)

1. Open `curves.csv`.
2. **Verify:** the table has 30 rows; the five curve columns are recognised as curves (semantic type
   `fit`); every cell of them shows a chart, not raw JSON text.
3. **Verify:** the console has no errors.

### 2. The Context Panel of a curve cell (GROK-20662, GROK-20664, public#2101)

1. Click a cell of `multiple styled series` (a cell with two series).
2. **Verify:** the Context Panel shows the **Options**, **Chart** and **Fit** panes.
3. Expand **Options**.
4. **Verify:** **Level** is `Column` and offers `Dataframe`, `Column`, `Cell`.
5. **Verify:** **Fit Function** offers `linear`, `sigmoid`, `log-linear`, `exponential`,
   `4pl-regression`, `4pl-dose-response` and has no empty entry (GROK-20662: clearing the fit
   function used to leave no way back).
6. **Verify:** **Droplines** offers `IC10`, `IC20`, `IC50`, `IC80`, `IC90`.
7. **Verify:** the **Statistics** list shows all its entries (`rSquared` … `pIC50`) without a
   scrollbar of its own (GROK-20664).
8. **Verify:** **Statistics mode** and **Aggregation** are shown: the cell holds several curves.
9. Click a cell of `styled proprtional with IC50` (one series).
10. **Verify:** **Fit Function** is shown; **Statistics mode** and **Aggregation** are not: one
    curve has nothing to summarise (public#2101).

### 3. The Edit chart dialog (public#2100)

1. Double-click a cell of `multiple prefit`.
2. **Verify:** the **Edit chart** dialog opens with the chart; it can be resized by its corner.
3. Click **Close**.
4. **Verify:** the dialog closes; the grid is unchanged; the console has no errors.

The rest of public#2100 (switching series and outliers inside the dialog) is not implemented — the
issue is open; do not test it.

### 4. Curve columns in other UI: filter panel and pinned column (GROK-17589, GROK-16166)

1. Open the **Filter Panel** (toolbar filter icon).
2. **Verify:** the panel opens; the console has no errors (GROK-17589).
3. Right-click the header of `multiple prefit` and choose **Pin > Pin Column**.
4. **Verify:** the column is pinned at the left and draws its charts (no empty cells); the console has
   no `Could not paint cell for pinned column` error (GROK-16166).

The rich curve field of the **Forms** viewer (GROK-19473, public#3534) is covered in
`Viewers/FormsViewer/forms.md`, Scenario 4.
