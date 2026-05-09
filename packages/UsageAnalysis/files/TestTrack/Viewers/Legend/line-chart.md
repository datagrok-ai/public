---
feature: legend
sub_features_covered:
  - legend.column
  - legend.refresh.on-data-change
target_layer: playwright
coverage_type: regression
priority: p1
pyramid_layer: integration
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Viewers/Legend/line-chart.md
migration_date: 2026-05-07
migration_report: line-chart-migration-report.md
related_bugs:
  - GROK-17222
  - GROK-17278
  - GROK-19083
---

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
3. Replace one Y column via the in-plot column selector (or set `lc.props.yColumnNames = ['Average Mass', 'NIBR logP']` if the in-plot selector is hover-only — see Unresolved ambiguities)
4. Save the layout, then re-apply it
5. Verify the new Y column choice and its updated legend block persist after the layout round-trip
6. Save the project, close it, and reopen it
7. Verify the new Y column choice and its updated legend block persist across the project round-trip

(Visual check «the corresponding legend block updates to reflect the new Y column's values» moved to `legend-ui.md` §3, 2026-05-08 — per-line legend blocks are not individually labeled in the DOM, single `[name="legend"]` element for the whole chart.)

## Notes

- Specialty: Line chart-specific legend behaviors — `Split` produces categorical legend, Multi Axis splits one legend into per-subplot blocks, multi-Y configurations support per-Y legend blocks, Y-column replacement updates the affected block.
- Delegates standard legend UI flows (visibility / position / color picker / save-dialog widgets) to `visibility-and-positioning.md`.
- `pyramid_layer: integration` per chain Rule 4 — multi-subsystem within Line chart (Split + Multi Axis + multi-Y configuration + persistence).
- Three related bugs:
  - GROK-17222 ("Line chart: legend is not consistent with filtering") — direct title match. This scenario does not include a filter step (the bug repro adds a filter); covered by the bug-focused spec `legend-grok-17222-spec.ts` (proposed by chain).
  - GROK-17278 ("Linechart: changed colors are not saved to the layout and project") — direct title match. This scenario does not change colors via legend (the bug repro changes a category color); covered by the bug-focused spec `legend-grok-17278-spec.ts`.
  - GROK-19083 (markers ↔ legend sync) — column intersection only; line chart does not use markers. Weak match.
- Helpers: `loginToDatagrok`, `softStep`, `specTestOptions`, `stepErrors` from `spec-login`. No new helpers proposed.
- Visual / pixel-level checks split into companion `legend-ui.md` (`target_layer: ui-only`, manual QA) on 2026-05-08: §1 «distinct palette colors», §2 «each Y line own subplot», §3 «per-Y legend block updates». Spec body retains JS-API proxy invariants (legend item count, prop equality) for these flows.

---
{
  "order": 6,
  "datasets": ["System:DemoFiles/SPGI.csv"]
}
