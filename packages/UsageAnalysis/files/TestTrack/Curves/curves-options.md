---
feature: curves
target_layer: playwright
pyramid_layer: ui-smoke
coverage_type: regression
priority: p1
related_bugs:
  - GROK-12239
  - GROK-12366
  - GROK-13420
  - GROK-14877
  - GROK-15703
  - GROK-17637
  - GROK-20663
  - public#2754
  - public#2924
---

# Curve options and the level they are set at

The **Options** pane of a curve cell edits series options (fit function, colors, markers, line
style, connect dots, fit line, confidence interval, droplines) and chart options (axis bounds, title,
axis names, log scales, legend, statistics). The **Level** input says where a change is stored:
`Column` (the default) reaches every cell of the column, `Cell` only the clicked cell, `Dataframe`
every curve column of the table. Column and Dataframe options are kept by layouts and projects.

## Setup

A clean session. Open `System:DemoFiles/curves.csv` from **Browse > Files > Demo**. The Context
Panel is open.

## Scenarios

### 1. A Column-level option reaches every cell of the column (GROK-12239)

1. Click a cell of `styled proprtional with IC50` in row 1. Expand **Options**; **Level** is
   `Column`.
2. Set **Fit Function** to `linear`.
3. Click the cell of the same column in row 5.
4. **Verify:** **Fit Function** reads `linear`.
5. Click a cell of `custom fit function`.
6. **Verify (manual):** **Fit Function** reads `Polynomial`: another column is not changed. If the
   Context Panel shows only **Actions** instead of Options / Chart / Fit, click the cell once more
   ([#4306](https://github.com/datagrok-ai/public/issues/4306): the panel loses the curve panes a
   moment after the click).
7. Go back to `styled proprtional with IC50` and set **Fit Function** back to `sigmoid`.

### 2. A Cell-level option stays in its cell (manual)

1. Click a cell of `styled proprtional with IC50` in row 1. Set **Level** to `Cell`.
2. Set **Title** to `Only this one`.
3. **Verify:** the cell draws the title when it is big enough (manual).
4. Click the cell in row 3 of the same column.
5. **Verify:** **Title** reads `Dose-Response curves`, the title the data names.
6. Click the cell in row 2 again.
7. **Verify:** **Title** reads `Only this one`.

### 3. Connect dots without and with the fit line (public#2924, public#2754)

1. Click a cell of `styled proprtional with IC50`. **Level** `Column`.
2. Uncheck **Show Fit Line**, check **Connect Dots**.
3. **Verify (manual):** every cell of the column draws straight segments between its points and no
   fitted curve; the console has no errors.
4. Check **Show Fit Line** again.
5. **Verify (manual):** both the connected points and the fitted curve are drawn.
6. Uncheck **Connect Dots**.

### 4. Markers, outliers, confidence interval (GROK-13420, GROK-14877, GROK-15703, GROK-12366)

1. Click a cell of `multiple styled series`. **Level** `Column`.
2. Set **Marker** to `square`, **Outlier Marker** to `diamond`, **Outlier Color** to red.
3. Check **Confidence Interval**.
4. **Verify (manual):** points are squares, outliers are red diamonds, each series has its own
   confidence band; the console has no errors.

### 5. Show Points: candlesticks need replicate points (GROK-20663)

1. Click a cell of `log min/max candlesticks` (several points share an x).
2. Set **Show Points** to `candlesticks`, then `both`.
3. **Verify (manual):** candlesticks are drawn for the shared x values.
4. On `styled proprtional with IC50` (one point per x) the same setting draws nothing new — that is
   by design (GROK-20663, Won't fix).

### 6. Column options are kept by a layout and a project (GROK-17637)

1. Click a cell of `styled proprtional with IC50`. **Level** `Column`. Set **Fit Function** to
   `log-linear`.
2. Save the layout, close the table, open `curves.csv` again and apply the layout.
3. Click a cell of `styled proprtional with IC50`.
4. **Verify:** **Fit Function** reads `log-linear`.
5. Save the project (**SAVE** on the ribbon), close all, open the project.
6. **Verify:** **Fit Function** of the column still reads `log-linear`; the console has no errors.
7. Delete the project and the layout.
