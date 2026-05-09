---
feature: charts
sub_features_covered: [charts.tree]
target_layer: playwright
coverage_type: regression
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Charts/tree.md
migration_date: 2026-05-07
migration_report: tree-migration-report.md
related_bugs: [github-3221, github-3245]
---

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

- Dataset: `System:DemoFiles/demog.csv` (declared in scenario footer, `order: 30`).
- Tree branches are canvas-rendered (ECharts); the canvas Shift+Click
  multi-select gesture is **moved to** `charts-ui.md` (ui-only manual
  scenario) per `charts-remediate-2026-05-09` user directive.
  Previously the spec used a `df.selection` bitset proxy with
  `[AMBIGUOUS]` warning; now the proxy is removed and the gesture is
  documented in `charts-ui.md` for human QA verification. This
  scenario covers Filter Panel control flow + filter cardinality
  coordination only.
- The original scenario does not specify the UI path for setting the Hierarchy (column-bag drag vs. property-panel hierarchy editor). Migrator preserved the original wording; see migration report's `Unresolved ambiguities`.
- Atlas declares 16 sub-features under `charts.tree` (layout, orient, initialTreeDepth, symbol, symbolSize, fontSize, moleculeSize, labelRotate, showCounts, color-palette, showMouseOverLine, size-aggregation, color-aggregation, on-click, includeNulls). The original scenario exercises **none** of those properties — this scenario covers only the base `charts.tree` viewer attached to a hierarchy plus the collaborative-filtering coordination with the Filter Panel. Property-coverage gaps (notably the github-3221 8-enhancement bundle and github-3245 Row Source × On Click state machine) are surfaced for follow-up bug-focused scenarios; see migration report.
