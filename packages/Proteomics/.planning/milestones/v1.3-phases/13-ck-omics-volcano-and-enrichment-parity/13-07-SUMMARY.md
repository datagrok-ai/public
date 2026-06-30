---
phase: 13-ck-omics-volcano-and-enrichment-parity
plan: 07
status: complete
gap_closure: true
requirements: [R2]
files_modified:
  - src/viewers/enrichment-viewers.ts
  - src/tests/enrichment-visualization.ts
tests:
  category: "Enrichment Visualization"
  before: 8
  after: 10
  added: ["openEnrichmentVisualization co-docks volcano on enrichment view (directional)", "openEnrichmentVisualization co-docks volcano on enrichment view (single-direction)"]
---

# 13-07 — Co-dock linked volcano on enrichment view

## What was built

`openEnrichmentVisualization` now co-docks a `proteinDf`-bound `createVolcanoPlot`
viewer LEFT of the dot/bar grid on the enrichment `TableView` in BOTH layout
branches (directional-split and single-direction fallback). The cross-DataFrame
selection wiring is untouched — `wireEnrichmentToVolcano` still mutates
`proteinDf.selection`, and the co-docked volcano renders that shared `BitSet`
without any new subscriptions.

## UX option chosen

Option B from the UAT diagnosis: keep the dot/bar charts on the enrichment view
but **co-dock a viewer of `proteinDf`** (specifically a linked volcano) on the
same view. Rationale:

- The directional charts already deserve their existing real estate; moving
  them to the protein view would have crowded the volcano the user opened
  through Visualize | Volcano.
- The protein `TableView` is untouched, so users coming from Show All
  Visualizations or the Visualize | Volcano menu see no regression.
- `wireEnrichmentToVolcano` keeps its single-responsibility
  `proteinDf.selection`-mutation contract — no new subscription endpoint, no
  side-effect coupling between layouts.

## createVolcanoPlot idempotency

`createVolcanoPlot` is a `create*` factory with no shell side-effects. Calling
it twice on the same `proteinDf` is safe because `volcano.ts` reuses
`NEG_LOG_COL` and `DIRECTION_COL` in place (Pitfall 5 — never remove+re-add),
and `applyThresholdLines` filters formula lines by prefix before adding new
ones. So re-opening the enrichment view on the same `proteinDf` does not stack
formula lines or duplicate columns.

## Guard

`createVolcanoPlot` throws when `adj.p-value` is missing
(`pickMetricColumn`). The `enrichmentCharts` fallback path in `package.ts`
passes the same `df` as both args when no DE-complete candidate exists, so
this case is real. Both branches wrap the co-dock in `try/catch` that warns
and continues — the enrichment view still renders without the linked volcano
in that case.

## Signature change

`openEnrichmentVisualization` now returns the `DG.TableView` it created or
reused. The return value is purely additive — existing callers (`package.ts`
`enrichmentAnalysis` / `enrichmentCharts`) ignore it. The new return value is
what makes the layout-regression tests reliable: it sidesteps the
`v.dataFrame === enrichDf` identity check on a fresh JS wrapper through
`grok.shell.tableViews` (the platform-internal `toJs` wrapper cache is not
guaranteed to round-trip JS identity for freshly-constructed DataFrames,
which made the first test pass fail — see the test-file inline comment).

## Tests added

- **Test B (directional-split)** — builds an enrichDf with Up/Up/Down/Down rows
  plus an `Intersection` column; calls `openEnrichmentVisualization`; asserts
  at least one viewer on the returned TV is bound to `proteinDf` and its type
  is `Scatter plot`; sets `enrichDf.currentRowIdx = 0` and asserts
  `proteinDf.selection.trueCount >= 1` (the data-layer wire still fires
  through the new dock).
- **Test C (single-direction fallback)** — same shape without a `Direction`
  column, exercising the second branch.

The existing 8 tests in the category (5 `createTopNEnrichmentDf` + 3
`wireEnrichmentToVolcano`) are unchanged — fixture extension on
`makeMockProteinDf` added an `adj.p-value` `Float32Array` column (bulk-init
via `DG.Column.fromFloat32Array` per memory `feedback_dg_column_bulk_init`)
but the existing tests don't read that column.

## Verification

- `npm run build` clean (no TypeScript errors).
- `grok test --category "Enrichment Visualization"` — **10/10 pass**.
- `grok test --category "Enrichment"` — **12/12 pass** (no regression to
  the data-layer cross-link contract in `src/tests/enrichment.ts`).

## Follow-up

13-UAT.md Test 2 should be re-run on a real Spectronaut Candidates DE result
to confirm Gap 1 closes from the user's perspective. The automated tests
cover the layout contract (dock present, viewer type, wire fires); the
visual confirmation that "highlight lands next to the dot/bar charts" is
human-UAT only.
