---
feature: legend
realizes_atlas: []
realizes: []
realized_as:
  - line-chart-spec.ts
target_layer: playwright
coverage_type: regression
priority: p1
pyramid_layer: integration
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Viewers/Legend/line-chart.md
migration_date: 2026-05-07
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities:
  - scenario-4-step-3-in-plot-column-selector-vs-ycolumnnames-js-api-path
  - scenario-4-step-4-the-corresponding-legend-block-must-update
scope_reductions: []
related_bugs:
  - GROK-17222
  - GROK-17278
  - GROK-19083
---

# Line chart — legend behavior with Split and Multi Axis

Verifies Line chart-specific legend behavior: setting **Split** produces a categorical legend, enabling **Multi Axis** splits that into a per-subplot legend block, and charting multiple Y columns gives each line its own legend block that updates correctly when a Y column is replaced. Also checks that all of this survives layout and project reloads.

## Setup

1. Open SPGI (`System:DemoFiles/SPGI.csv`)
2. Add a **Line chart**

## Scenarios

### 1. `Split = Series` produces a categorical legend

1. On the Line chart, set **Split** = `Series`
2. Verify the legend lists every distinct value of `Series` (typically 4–10 categories)

(Visual check «each category gets a distinct palette color, no two adjacent share» moved to `legend-ui.md` §1, 2026-05-08 — pixel-level palette ordering is undecidable at the JS API surface.)

### 2. Multi Axis — per-subplot legend driven by `Split`

1. On the **Context Panel > Data**, enable **Multi Axis** (`lc.props.multiAxis = true`)
2. Verify each subplot has its own legend driven by `Split` (per-subplot legend block, not a single shared legend)

(Visual check «each Y line gets its own subplot» moved to `legend-ui.md` §2, 2026-05-08 — subplot SVG layout requires DOM measurement that the Line chart does not expose.)

### 3. Layout and project persistence — Split + Multi Axis + per-line legends [coverage_type: edge]

1. Save the layout
2. Re-apply the saved layout (allow ≥3 s settle)
3. Verify `Split`, `Multi Axis`, and per-line legends all survive the layout round-trip
4. Save the project, close it, and reopen it
5. Verify `Split`, `Multi Axis`, and per-line legends all survive the project round-trip

### 4. Multi-Y configuration — both lines render their own legend block, Y replacement updates the affected legend [coverage_type: edge]

1. Configure two Y columns: `lc.props.yColumnNames = ['Average Mass', 'TPSA']`
2. Verify both lines render and each gets its own legend block driven by `Split`
3. Replace one Y column via the in-plot column selector (or set `lc.props.yColumnNames = ['Average Mass', 'NIBR logP']` if the in-plot selector is hover-only)
4. Save the layout, then re-apply it
5. Verify the new Y column choice and its updated legend block persist after the layout round-trip
6. Save the project, close it, and reopen it
7. Verify the new Y column choice and its updated legend block persist across the project round-trip

(Visual check «the corresponding legend block updates to reflect the new Y column's values» moved to `legend-ui.md` §3, 2026-05-08 — per-line legend blocks are not individually labeled in the DOM, single `[name="legend"]` element for the whole chart.)

## Notes

- Standard legend UI flows (visibility, position, color picker, save-dialog behavior) are covered once in `visibility-and-positioning.md` and not repeated here.
- Related bugs:
  - GROK-17222 ("Line chart: legend is not consistent with filtering") — this scenario doesn't include a filter step; the filtering-specific repro is covered separately by `legend-grok-17222-spec.ts`.
  - GROK-17278 ("Linechart: changed colors are not saved to the layout and project") — this scenario doesn't change colors via the legend; that repro is covered separately by `legend-grok-17278-spec.ts`.
  - GROK-19083 (markers <-> legend sync) — weak match: the Line chart doesn't use markers.
- The visual, pixel-level checks (distinct palette colors, each Y line getting its own subplot, per-Y legend block updates) were moved to the companion manual-QA file `legend-ui.md` (sections 1-3) — they can't be verified through the DOM. This file keeps the checks that can be verified through the JS API / UI (legend item count, property equality).

---
{
  "order": 6,
  "datasets": ["System:DemoFiles/SPGI.csv"]
}
