---
feature: charts
target_layer: playwright
coverage_type: smoke
priority: p0
realizes_atlas: [charts.cp.open-viewer-with-required-columns, charts.cp.configure-via-property-panel]
realizes: [charts.chord, charts.globe, charts.group-analysis, charts.multiplot, charts.radar, charts.sankey, charts.sunburst, charts.surface-plot, charts.tree, charts.word-cloud]
realized_as:
  - radar-spec.ts
pyramid_layer: ui-smoke
ui_coverage_responsibility:
  - add-viewer-radar
  - viewer-property-panel-gear
ui_coverage_delegated_to: null
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Charts/radar.md
migration_date: 2026-05-07
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
related_bugs:
  - GROK-18085
---

# Radar viewer (Charts package)

Smoke scenario for the Charts package Radar viewer. Verifies the
"Add Viewer → Radar" entry point on two distinct datasets, and the
property panel (Gear icon → Context Panel) for four representative
properties: switching the bound table, selection check-boxes, the
chosen Values count, and style/color changes. This is the primary
Add-Viewer + property-panel smoke test for the Charts section;
sibling scenarios `sunburst.md` and `tree.md` cover their own
specialty UI surfaces.

Note: GROK-18085 (Radar table-rebind breaking on project save/reopen)
is related to this viewer but is not exercised by this scenario — see
`radar-save-reopen-bug.md` for that reproduction.

## Setup

A clean Datagrok session is the only shared setup. Each step opens
its own dataset; no fixture chaining required (`depends_on: []` per
chain rev 1).

## Scenarios

### Add Radar viewer to two datasets and exercise the property panel

1. Open `System:DemoFiles/geo/earthquakes.csv` (e.g. via
   **Browse** > **Files**, or **File** > **Open**).
2. On the Menu Ribbon, click **Add viewer** and select **Radar**.
3. **Verify:** the Radar viewer opens for `earthquakes.csv` with
   no errors.
4. Open `System:DemoFiles/demog.csv`.
5. On the Menu Ribbon, click **Add viewer** and select **Radar**.
6. **Verify:** the Radar viewer opens for `demog.csv` with no errors.
7. On the Radar viewer (on `demog.csv`), click the **Gear** icon.
8. **Verify:** the **Context Panel** opens with the Radar viewer's
   properties.
9. **Switching tables:** in the Context Panel, change the bound
   `table` property between `earthquakes` and `demog`.
   **Verify:** the Radar viewer rebinds and re-renders against the
   new table without errors.
10. **Check-boxes in selection:** select two or more rows in the
    bound table (Ctrl+Click in the grid, or
    `df.selection.set(...)` equivalent UI gesture).
    **Verify:** the selected-row lines are reflected on the Radar
    viewer (current / mouseover row lines).
11. **Increasing and decreasing the amount of chosen Values:** in
    the Context Panel, change the columns chosen as Values (the
    color column / values bag).
    **Verify:** the Radar viewer redraws with more / fewer axes
    matching the chosen Values.
12. **Style (color) changes:** in the Context Panel, change a
    Style-category color (e.g. `currentRowColor`,
    `mouseOverRowColor`, `lineColor`, or `backgroundMinColor` /
    `backgroundMaxColor`).
    **Verify:** the Radar viewer redraws with the new color.
13. **Check all properties (broad sweep):** scroll through the
    Context Panel and confirm every Radar property is visible
    and editable. Spot-check that toggling representative
    properties does not throw console errors. (See Notes — this
    step preserves the original scenario's broad "Check all
    properties" instruction without enumerating each radar.*
    property; the four bullets above are the explicit Main things
    the original calls out as MUST be reflected on the viewer.)

## Notes

- Step 13 ("check all properties") is a representative sweep, not an
  exhaustive one — it doesn't enumerate every `radar.*` property
  (title, min/max percentile, show-current-row, show-tooltip,
  color-column, color-palette, show-min-max, legend-visibility)
  individually. Steps 3 and 9-12 already verify the specific
  properties the manual scenario calls out; the sweep covers the
  rest at a shallower level.

---
{
  "order": 28,
  "datasets": ["System:DemoFiles/demog.csv"]
}
