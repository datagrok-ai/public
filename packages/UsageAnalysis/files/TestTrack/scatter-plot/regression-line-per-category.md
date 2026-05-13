---
feature: scatter-plot
sub_features_covered:
  - scatter-plot.regression-line.per-category
  - scatter-plot.show-regression-line
  - scatter-plot.color-column-name
  - scatter-plot.regression-line.show-equation
  - scatter-plot.regression-line.pearson-correlation
target_layer: uitests-package
priority: regression
produced_from: atlas-driven
related_bugs: []
---

# Scatter Plot — Regression line per category

## Setup

1. Close all open views and viewers.
2. Open the `demog` demo dataset (`System:DemoFiles/demog.csv`).
3. Add a Scatter Plot viewer.
4. In the property panel, set X to `AGE` and Y to `HEIGHT`.

## Scenarios

### Scenario 1: Per-category regression line driven by categorical Color (regression)

Priority: regression

Steps:
1. Press `R` to enable the regression line — verify a single line
   appears across the entire plot.
2. In the property panel, set Color column to `RACE` (categorical with
   four values in `demog`).
3. Enable the Per-category regression line property
   (`scatter-plot.regression-line.per-category`).
4. Enable Show regression line equation.

Expected:
- Step 3: a separate regression line is drawn for each `RACE`
  category, color-matched to the legend swatch for that category.
- The Per-category property only takes effect when the Color column is
  categorical (numerical Color must collapse back to a single line —
  see Scenario 2).
- Step 4: an equation is drawn next to each per-category line.

### Scenario 2: Per-category falls back to single line when Color is numerical (edge)

Priority: edge

Steps:
1. With the setup and `R` toggled on, enable the Per-category
   regression line property.
2. Set Color column to `AGE` (numerical).
3. Enable Show Pearson correlation.

Expected:
- The viewer renders a single regression line, not one per bin —
  per-category mode is meaningful only for categorical Color.
- Pearson correlation is shown alongside the single line.
- Switching Color column back to `RACE` (categorical) returns the
  viewer to per-category mode without re-toggling the property.

## Notes

- **gap_description:** `"regression line not tested per-category"` —
  STEP C-1 matched `scatter-plot.regression-line.per-category`
  directly. Parent (`scatter-plot.show-regression-line`) and the
  triggering interaction (`scatter-plot.color-column-name`, must be
  categorical) were pulled in to satisfy density rule A-STRUCT-06.
- **target_layer rationale:** every step is property-panel work on a
  single viewer with no cross-dialog state — `uitests-package` per
  STEP D.
- **Coverage map status:** at scan time `version_of_scan: 2026-04-27`,
  `scatter-plot.regression-line.per-category` had
  `coverage_verdict: uncovered`. After this scenario passes, the
  Coverage Mapper should move it to `partial` (still missing the
  categorical-with-many-categories visual stress case, deferred
  below).
- **Deferrals:**
  - Visual quality of overlapping per-category lines under heavy
    overlap — deferred, requires pixel-precision visual diff which
    `uitests-package` cannot drive (per A-MERIT-02 / Lattice Rule 13,
    real dependency cited).
- **Atlas notes:** the per-category sub_feature has no
  `interactions[]` entries and no `keyboard_shortcuts[]`. The
  property name in the property panel is sourced from
  `scatterplot_look.dart:L267` per atlas; Automator should resolve
  the exact label string at spec time.
