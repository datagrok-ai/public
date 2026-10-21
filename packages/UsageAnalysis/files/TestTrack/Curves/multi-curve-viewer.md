---
feature: curves
target_layer: playwright
pyramid_layer: ui-smoke
coverage_type: regression
priority: p0
related_bugs:
  - GROK-14069
  - GROK-19014
  - GROK-19025
  - GROK-19028
  - GROK-20112
  - GROK-20659
  - GROK-20661
  - public#1645
---

# Multi curve viewer

The **MultiCurveViewer** (Add Viewer: `MultiCurveViewer`) overlays curves from several rows and
curve columns on one plot: the selected rows, the current row and the row under the mouse. Its
**Curves** property picks the curve columns, **Legend** names the curves by a column. Its other
properties are the chart options of a curve cell (axis bounds, log scales, title, axis names, show
legend, statistics), plus **Merge Column Series**, **Merge Cell Series** and **Show Outliers**.

## Setup

A clean session.

## Scenarios

### 1. Added from the Add Viewer dialog (public#1645)

1. Open `System:DemoFiles/curves.csv`.
2. Click **Add viewer** on the toolbar, search `curve`.
3. **Verify:** the `MultiCurveViewer` card is enabled.
4. Click it.
5. **Verify:** the dialog closes and the view holds one MultiCurveViewer; the console has no errors.
6. Click its **Gear** icon.
7. **Verify:** **Curves** lists one of the five curve columns (`1 / 5`).

### 2. Not offered for a table without curves (GROK-19014, GROK-14069)

1. Open `System:DemoFiles/demog-1000.csv`.
2. Click **Add viewer**, search `curve`.
3. **Verify:** the `MultiCurveViewer` card is disabled. Press **Esc**: the dialog closes, nothing is
   added and the console has no errors.

### 3. Curve columns, selection and the current row

1. On `curves.csv` with a MultiCurveViewer, open **Curves** (`...`) and click **All**, **OK**.
2. **Verify:** **Curves** reads `5 / 5`; the console has no errors.
3. Select rows 1–3 in the grid.
4. **Verify (manual):** the viewer overlays the curves of the three rows from every curve column;
   the Y axis takes in all of them (GROK-20661 — none is cut at the top).
5. Uncheck **Show Selected Rows Curves**.
6. **Verify (manual):** only the current row's curves stay.

### 4. Merge Cell Series on and off leaves the grid as it was (GROK-19028) — manual

1. Click row 1. Check **Merge Cell Series**, then uncheck it.
2. **Verify:** the viewer shows the series as before; the grid cell of row 1 in `multiple prefit`
   still holds two series (its **Options** pane: **Merge Series** unchecked, **Fit** pane lists two
   series).
3. **Verify:** the console has no errors.

### 5. Use in Trellis (GROK-19025)

1. On the viewer, open the hamburger menu > **General** > **Use in Trellis**.
2. **Verify:** a Trellis plot is added with **Curves** as the inner viewer; its Context Panel has a
   **MultiCurveViewer** tab with the viewer's properties.

### 6. The legend and its tooltips (GROK-20112, GROK-20659)

Covered by the guide `Curves/bdd/features/guides/curves-legend.feature` on the Assay Curves demo:
the legend takes the free corner, a shortened name gives the whole of it on hover, `+N more` lists
the rest, **Show Legend** hides and shows it.

### 7. Closing the viewer (Curves 1.13.1)

1. Close the MultiCurveViewer (its **X** icon). Select some rows, click other rows.
2. **Verify:** the console has no errors (a closed viewer used to keep re-rendering).
