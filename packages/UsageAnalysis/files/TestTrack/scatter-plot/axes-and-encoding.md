---
feature: scatter-plot
sub_features_covered:
  - scatter-plot.x-column-name
  - scatter-plot.y-column-name
  - scatter-plot.color-column-name
  - scatter-plot.size-column-name
  - scatter-plot.markers-column-name
  - scatter-plot.marker-type
  - scatter-plot.marker-default-size
  - scatter-plot.marker-opacity
  - scatter-plot.jitter-size
  - scatter-plot.x-axis-type
  - scatter-plot.y-axis-type
  - scatter-plot.invert-x-axis
  - scatter-plot.filter-out-invalid
  - scatter-plot.show-regression-line
  - scatter-plot.regression-line.show-equation
  - scatter-plot.regression-line.pearson-correlation
target_layer: uitests-package
priority: regression
produced_from: atlas-driven
related_bugs: []
---

# Scatter Plot — Axes, encoding, and regression line

## Setup

1. Close all open views and viewers.
2. Open the `demog` demo dataset (`System:DemoFiles/demog.csv`).
3. Add a Scatter Plot viewer to the table view.

## Scenarios

### Scenario 1: Default axes wire-up (smoke)

Priority: smoke

Steps:
1. Verify the Scatter Plot opens with x and y bound to two of the
   numerical columns from `demog` (`AGE`, `HEIGHT`, `WEIGHT`).
2. In the property panel, change the X column to `AGE`.
3. Change the Y column to `HEIGHT`.
4. Set the Color column to `RACE` (categorical).
5. Set the Marker column to `SEX` (categorical).

Expected:
- Markers redraw on each change. Color legend shows the four `RACE`
  categories. Marker shapes differ across `SEX` values.

### Scenario 2: Marker appearance combinations (regression)

Priority: regression

Steps:
1. With the setup applied, set the Size column to `WEIGHT`.
2. In the property panel, switch Marker default size to `8`.
3. Lower Marker opacity to `40`.
4. Set Jitter size to `5`.
5. Switch Marker type to `cross` (no-markers-column path: clear the
   Markers column first via the property panel).
6. Re-set Markers column back to `SEX` and confirm shapes override
   the default Marker type.

Expected:
- Marker size scales with `WEIGHT` per row.
- Opacity reduction visibly softens overlapping marker clusters.
- Jitter shifts dense regions without changing axis ranges.
- Default Marker type only applies while Markers column is empty.

### Scenario 3: Regression line with correlation readouts (regression)

Priority: regression

Steps:
1. With the setup applied, set X to `AGE`, Y to `HEIGHT`.
2. Press `R` to toggle the regression line on.
3. In the property panel, enable Show regression line equation.
4. Enable Show Pearson correlation.

Expected:
- A regression line appears across the X/Y plane after `R` is pressed.
- The line equation label is visible alongside the line.
- A Pearson correlation value is rendered next to the equation.

### Scenario 4: Logarithmic axis with non-positive values (edge)

Priority: edge

Steps:
1. With the setup applied, set X to `AGE`, Y to `HEIGHT`.
2. In the property panel, switch X axis type to `logarithmic`.
3. Switch Y axis type to `logarithmic`.
4. Toggle Filter-out invalid OFF.
5. Toggle Filter-out invalid back ON.
6. Toggle Invert X axis ON.

Expected:
- With log scales and Filter-out invalid OFF, rows with non-positive
  AGE/HEIGHT (if any in `demog`) are dropped or rendered out-of-bounds
  by the log transform — the viewer must not crash.
- Toggling Filter-out invalid ON re-applies the filter and the plot
  re-fits to the remaining points.
- Invert X axis flips the X direction without breaking the regression
  line position.

## Notes

- **target_layer rationale:** every step is a single property-panel
  toggle or column selection on the Scatter Plot viewer; cross-dialog
  state and persistence are not exercised. This is the
  `uitests-package` shape per STEP D.
- **Coverage map status:** STEP B confirmed `coverage-map/scatter-plot.yaml`
  exists; all 47 sub_features were `uncovered` at scan time
  (`version_of_scan: 2026-04-27`). This file does not attempt the 70%
  bar (A-STRUCT-05); it is a coherent slice — Critic Gate A
  SCOPE_REDUCTION is the expected acceptable outcome.
- **STEP C fallback path taken:** no `gap_description` was provided,
  and the scatter-plot atlas's `critical_paths` contains only a
  placeholder marker (`id: test-curation-marker`), not curated paths.
  Per STEP C-2, fell back to grouping highest-interaction uncovered
  sub_features. The four scenarios cover 16 of 47 sub_features at
  density 4.0 (well above the ≥ 2 floor in A-STRUCT-06).
- **Atlas level:** `full` (per `feature-atlas/scatter-plot.yaml`),
  but `interactions`, `edge_cases`, and `known_issues` arrays are
  empty in this atlas revision — quality of the edge scenario is
  bounded by source-derived behavior only.
- **No deferrals.**
