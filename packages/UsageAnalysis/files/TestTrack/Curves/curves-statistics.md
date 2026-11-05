---
feature: curves
target_layer: playwright
pyramid_layer: ui-smoke
coverage_type: regression
priority: p0
related_bugs:
  - GROK-12140
  - GROK-13408
  - GROK-16140
  - GROK-17492
  - GROK-17637
  - public#1751
  - public#2103
  - public#2310
---

# Fit statistics and statistic columns

The **Fit** pane of a curve cell lists the statistics of each series, named by the series' own fit
function (`R²`, `AUC`, `Top`, `Slope`, `IC50`, `Bottom`, `Y at IC50`, `Max Y`, `Min Y`, `pIC50` for a
sigmoid). The **+** next to a statistic adds it for the whole column as a calculated column, placed
right after the curve column. For a cell with several series, **Series** = `aggregated` and
**Aggregation type** summarise the series (an averaged IC50 is a geometric mean, GROK-17637).

## Setup

A clean session. Open `System:DemoFiles/curves.csv` from **Browse > Files > Demo**. The Context
Panel is open.

## Scenarios

### 1. A statistic of one series becomes a column (GROK-12140, public#1751, public#2103)

1. Click a cell of `styled proprtional with IC50`. Expand **Fit**.
2. **Verify:** the pane shows the series name and its statistics, among them **IC50** and **R²**.
3. Click **+** next to **IC50** (tooltip `Calculate ic50 for the whole column`).
4. **Verify:** a numeric column `styled proprtional with IC50 1 ic50` is added right after the curve
   column, filled for the rows whose curve fits; the console has no errors.

### 2. An aggregated statistic becomes a column (public#2310, GROK-16140)

1. Click a cell of `multiple styled series` (two series). Expand **Fit**.
2. **Verify:** **Series** offers `all` and `aggregated`.
3. Set **Series** to `aggregated` and **Aggregation type** to `avg`.
4. Click **+** next to **IC50** (tooltip `Calculate ic50 avg for the whole column`).
5. **Verify:** a numeric column `multiple styled series avg ic50` is added right after the curve
   column; the console has no errors.

### 3. A series without a name (GROK-17492)

1. Open a table with a curve column whose series have no `"name"` (for example, add the column
   `curvecol` of `curves-data-robustness.md`, Scenario 1).
2. Click a cell of it. In **Fit**, click **+** next to **R²**.
3. **Verify:** a statistic column `curvecol … rSquared` is added; the console has no errors.

### 4. Toggling an outlier refits the curve and its statistic column (GROK-13408, GROK-17637)

1. Add the `ic50` column of `styled proprtional with IC50` (Scenario 1).
2. Make the curve column wide and the rows tall enough for the points to be clickable.
3. Click a point of the curve in row 1.
4. **Verify (manual):** the point is drawn as an outlier, the curve is refitted, and the value of
   the `ic50` column in row 1 changes; clicking the point again restores both.

Manual: the points have no names a test can click.

### 5. Statistics drawn on the plot (public#2101)

1. Click a cell of `multiple styled series`. In **Options**, check `ic50` and `rSquared` under
   **Statistics**; set **Statistics mode** to `aggregated`, then `both`.
2. **Verify (manual):** the plot shows the statistics per series, then one summary line, then both;
   the font of the statistics does not change when the filter panel opens (GROK-20664).
