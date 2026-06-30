---
phase: 13-ck-omics-volcano-and-enrichment-parity
plan: 09
status: complete
gap_closure: true
requirements: [R2]
files_modified:
  - src/viewers/enrichment-viewers.ts
  - src/tests/enrichment-visualization.ts
tests:
  category: "Enrichment Visualization"
  before: 10
  after: 13
  added:
    - "openEnrichmentVisualization docks dot/bar on protein view (directional)"
    - "openEnrichmentVisualization docks dot/bar on protein view (single-direction)"
    - "openEnrichmentVisualization switches focus to protein TableView"
---

# 13-09 — Option A: dock enrichment charts on the protein view

## What was built

`openEnrichmentVisualization` now resolves the protein TableView by
Dart-handle identity scan over `grok.shell.tableViews`, computes
`chartHost = proteinTv ?? tv`, and docks the Up/Down dot+bar viewers (or
the single-direction fallback dot+bar pair) on `chartHost`. After all
docking completes, `grok.shell.v = proteinTv` switches focus back to the
view holding the user's existing volcano.

The 13-07 `createVolcanoPlot(proteinDf, {title: 'Volcano (linked)'})`
co-dock is removed from both layout branches, along with the
`import {createVolcanoPlot} from './volcano';` it depended on. The
smart-filter banner stays on `tv` (the enrichment view) — it describes
the merged enrichment DataFrame, not the protein DataFrame.

The data-layer cross-link (`wireEnrichmentToVolcano` mutating
`proteinDf.selection`) is unchanged; only the dock target and post-dock
focus change.

## UX option chosen

**Option A** — dock the dot/bar viewers on the protein TableView; the
user's pre-existing volcano (created by Proteomics | Visualize | Volcano
Plot or the DE auto-open) is the cross-link target.

The post-13-07 UAT round (Test 2, severity major) confirmed Option B (a
co-docked "Volcano (linked)" at ratio 0.4 on the enrichment view) was
insufficient:
- A duplicate volcano on a separate view did not match the user's mental
  model of "their" volcano (the one on the table they imported).
- The user still landed on the enrichment view after the dialog closed,
  so the dot/bar charts were still physically associated with it.
- The verbatim complaint — "enrichment charts are associated with the
  enrichment table and not with the table that has the volcano plot" —
  is a direct request for the charts to live on the protein view.

## Fallback semantics

When `proteinTv` cannot be resolved (e.g. `enrichmentCharts` fallback
where `df` is passed as both args, or the banner-only tests that don't
call `addTableView(proteinDf)`), `chartHost` evaluates to `tv` and the
function falls back to the pre-13-09 dock target. The `if (proteinTv)`
guard on the focus switch ensures it never throws when the protein view
is absent. The two banner tests (`enrichmentBannerRendersWhenTagSet`,
`enrichmentBannerAbsentWhenTagUnset`) intentionally do not pre-stage a
protein view, so they now exercise this fallback path as a free
regression net on top of their original banner-presence assertions.

## Test updates

- **Test B (directional)** — opens proteinTv via `grok.shell.addTableView`
  before the call, snapshots `countAllViewers(proteinTv)`, then asserts
  the count grew by ≥4 (Up dot + Up bar + Down dot + Down bar) and that
  `enrichTv` carries 0 proteinDf-bound viewers (the 13-07 co-dock is
  gone). Cross-link assertion preserved.
- **Test C (single-direction)** — same shape, asserts ≥+2 viewers added
  to the protein view.
- **Test D (focus, new)** — asserts
  `(grok.shell.v?.dataFrame as any)?.dart === (proteinDf as any).dart`
  after `openEnrichmentVisualization` returns.
- Helper `countAllViewers` added (because `tv.viewers` is `Iterable`,
  not array — no `.length`); unused `firstViewerBoundTo` removed.

## Verification

- `npm run build` succeeds with no TypeScript errors.
- `grok test --category "Enrichment Visualization"` — 13/13 pass
  (5 createTopNEnrichmentDf + 3 wireEnrichmentToVolcano + 3 layout/focus
  + 2 banner).

## Human-UAT follow-up

13-UAT Test 2 should be re-run on a real BP DMD/WT Spectronaut Candidates
DE result. Expected: after Enrichment runs, the user lands on the
protein view; the Up/Down dot+bar viewers are docked to the right of
their existing volcano on that view; clicking term rows on either the
dot+bar (now on the protein view) or the merged enrichment grid (on the
enrichment tab) highlights matching proteins on the original volcano in
the focused view. This closes the verbatim "enrichment charts associated
with the enrichment table" complaint.
