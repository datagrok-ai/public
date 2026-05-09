---
feature: legend
sub_features_covered:
  - legend.use-custom-color-coding
  - legend.item.color-picker
  - legend.allow-item-coloring
target_layer: playwright
coverage_type: edge
priority: p0
pyramid_layer: integration
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Viewers/Legend/color-consistency.md
migration_date: 2026-05-07
migration_report: color-consistency-migration-report.md
related_bugs:
  - GROK-17438
  - github-3132
  - GROK-17278
---

## Setup

1. Open SPGI (`System:DemoFiles/SPGI.csv`)
2. Add six viewers via the JS API (omit Scatter plot — it is covered by `scatterplot.md`): **Histogram**, **Line chart**, **Bar chart**, **Pie chart**, **Trellis plot**, **Box plot**
3. Set the categorical legend on each viewer to `Stereo Category`, using the viewer's own legend-source property:
   * Histogram, Line chart, Bar chart → **Split**
   * Pie chart → **Category**
   * Trellis plot → **X**
   * Box plot → **Category**

## Scenarios

### 1. Categorical color coding from grid propagates to every viewer's legend

1. In the grid, enable **Categorical color coding** for `Stereo Category`
2. Change at least two category colors in the grid (e.g. `R_ONE` → red, `S_UNKN` → green)
3. Verify every viewer's legend reflects the new colors (the column is the single source of truth for legend coloring)

### 2. Per-category color change via legend propagates back to grid and other viewers

1. On the Bar chart, hover a legend swatch and click the color picker icon
2. Pick a new color for one category and click **OK**
3. Verify the change propagates to the grid's column color coding
4. Verify the change propagates to every other viewer's legend (column is single source of truth, not a per-viewer state)

### 3. Layout persistence — customized palette survives reload

1. Save the layout
2. Re-apply the saved layout (allow ≥3 s settle for legend rebuild)
3. Verify every viewer still shows the customized palette after the round-trip

### 4. Project persistence — customized palette survives close/reopen

1. Save the project
2. Close the project, then reopen it
3. Verify the customized palette persists across the project round-trip on every viewer

## Notes

- Specialty: 6-viewer color propagation; grid as single source of truth for legend coloring. Scatter plot deliberately omitted (covered by `scatterplot.md` Section 5).
- Delegates standard legend UI flows (visibility / position / save-dialog widgets) to `visibility-and-positioning.md`.
- `pyramid_layer: integration` per chain Rule 4 — multi-viewer co-existence verifying bidirectional sync grid↔viewers↔legend.
- `priority: p0` — golden-path color sync. The "single source of truth" invariant is foundational; if it breaks, every per-viewer color customization workflow degrades.
- Three related bugs:
  - GROK-17438 (legend disappears after color change on shared-legend viewers) — strong cross-cutting match.
  - github-3132 (sequential color changes reset previous) — directly tested by Scenario 2 (one color change is the positive case; bug-focused spec covers sequential).
  - GROK-17278 (linechart colors not saved to layout/project) — Scenarios 3 and 4 verify the positive baseline.
- Helpers: `loginToDatagrok`, `softStep`, `specTestOptions`, `stepErrors` from `spec-login`. No new helpers proposed.

---
{
  "order": 4,
  "datasets": ["System:DemoFiles/SPGI.csv"]
}
