---
feature: charts
target_layer: playwright
coverage_type: regression
priority: p1
realizes_atlas: [charts.cp.open-viewer-with-required-columns, charts.cp.configure-via-property-panel]
realizes: [charts.tree]
realized_as:
  - tree-spec.ts
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Charts/tree.md
migration_date: 2026-05-07
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities:
  - setup-step-set-the-hierarchy-to-control-sex-and-race
  - filtered-count-semantics
scope_reductions: []
related_bugs: [github-3221, github-3245]
---

# Tree viewer — Filter Panel coordination

Checks that the Tree viewer, configured with a CONTROL → SEX → RACE
hierarchy, stays in sync with the Filter Panel: applying and clearing
a categorical filter updates the visible row count correctly and the
viewer keeps rendering without errors.

## Setup

1. Open `System:DemoFiles/demog.csv` (5850 rows).
2. On the Menu Ribbon, click **Add viewer** and select **Tree**. Tree viewer attaches to the demog table view.
3. Set the Tree viewer's **Hierarchy** to columns **CONTROL**, **SEX**, and **RACE** (in that order).

## Scenarios

### Filter Panel coordination — Tree viewer ↔ Filter Panel

Verifies that the Filter Panel toggles drive the dataframe filter
bitset and the Tree viewer reflects the filtered row set without
errors. The original scenario interleaved Tree branch Shift+Click
(canvas) with Filter Panel toggles to test the
selection-∧-filter intersection contract — the canvas
Shift+Click portion is **moved to** `charts-ui.md` (ui-only
manual scenario) since canvas hit-testing has no automatable
equivalent. This playwright-layer scenario covers the Filter
Panel control flow + filter cardinality coordination only.

1. Open the **Filter Panel** (via Toolbox `Filters` section header
   or ribbon Filter icon).
2. In the Filter Panel, set **CONTROL = true** (apply categorical
   filter via the categorical control or `df.filter` bitset since
   categorical checkboxes are canvas-rendered per `filters.md`).
3. Verify that the filter cardinality is **39** (rows with
   CONTROL=true in demog.csv).
4. In the Filter Panel, clear the **CONTROL = true** filter.
5. Verify that the filter cardinality returns to **5850** (full
   dataset).

> **Tree branch Shift+Click multi-selection** (the original
> Steps 1, 4 of the scenario) is moved to `charts-ui.md` —
> manual-only ui-only scenario. The selection-∧-filter
> cardinality contract (originally asserted as count=0 at
> CONTROL=true ∧ initial selection, count=2 after extension,
> count=176 after filter clear) is preserved in the manual
> scenario for human QA verification.

## Notes

- Tree branches are canvas-rendered (ECharts); the canvas Shift+Click
  multi-select gesture is moved to `charts-ui.md` (ui-only manual
  scenario). The spec previously used a `df.selection` bitset proxy
  as an ambiguous substitute; that proxy has been removed, and the
  gesture is now documented in `charts-ui.md` for human QA
  verification. This scenario covers Filter Panel control flow +
  filter cardinality coordination only.
- The original scenario doesn't specify the UI path for setting the
  Hierarchy (column-bag drag vs. property-panel hierarchy editor) —
  this ambiguity is preserved as-is in the wording.
- This scenario covers only the base Tree viewer attached to a
  hierarchy, plus Filter Panel coordination — it doesn't exercise
  most of the Tree viewer's ~16 configurable properties (layout,
  orient, initialTreeDepth, symbol/symbolSize, fontSize,
  moleculeSize, labelRotate, showCounts, color-palette,
  showMouseOverLine, size/color-aggregation, onClick, includeNulls).
  The github-3221 8-enhancement bundle and github-3245 Row Source ×
  On Click state machine are covered separately — see
  `tree-enhancements-bundle-bug.md` and
  `tree-rowsource-onclick-state-bug.md`.
