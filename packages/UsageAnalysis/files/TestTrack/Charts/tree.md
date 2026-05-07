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

### Collaborative filtering — Tree viewer ↔ Filter Panel

Verifies that branch selection in the Tree viewer and filter toggles in the Filter Panel coordinate consistently: the visible row set is always the intersection of the Tree-driven selection bitset and the Filter-Panel-driven filter bitset, with the row counts asserted at each transition.

1. In the Tree viewer, hold **Shift** and click to multi-select the following three branches:
   - **All → false → F → Asian**
   - **All → false → F → Black**
   - **All → false → M → Asian**
2. In the **Filter Panel**, set **CONTROL = true**.
3. Verify that the filtered row count (selection ∧ filter) is **0**.
4. In the Tree viewer, extend the selection by **Shift + Click** on the branch **All → true → F → Black**.
5. Verify that the filtered row count (selection ∧ filter) is **2**.
6. In the **Filter Panel**, clear the **CONTROL = true** filter.
7. Verify that the filtered row count is **176**.

## Notes

- Dataset: `System:DemoFiles/demog.csv` (declared in scenario footer, `order: 30`).
- Tree branches are canvas-rendered (ECharts); the existing `tree-spec.ts` documents canvas-coord Shift+Click as not reliably synthesizable and uses a programmatic `df.selection` bitset fallback. This migration preserves the original UI-step wording verbatim and leaves the synthesis decision to Automator (current spec authoritative behavior is the bitset fallback with `[AMBIGUOUS]` warning).
- The original scenario does not specify the UI path for setting the Hierarchy (column-bag drag vs. property-panel hierarchy editor). Migrator preserved the original wording; see migration report's `Unresolved ambiguities`.
- Atlas declares 16 sub-features under `charts.tree` (layout, orient, initialTreeDepth, symbol, symbolSize, fontSize, moleculeSize, labelRotate, showCounts, color-palette, showMouseOverLine, size-aggregation, color-aggregation, on-click, includeNulls). The original scenario exercises **none** of those properties — this scenario covers only the base `charts.tree` viewer attached to a hierarchy plus the collaborative-filtering coordination with the Filter Panel. Property-coverage gaps (notably the github-3221 8-enhancement bundle and github-3245 Row Source × On Click state machine) are surfaced for follow-up bug-focused scenarios; see migration report.
