---
status: diagnosed
trigger: "enrichment-cross-link-not-firing: enrichment dot/bar charts feel associated with the enrichment table and not with the table that has the volcano plot — Phase-9 cross-DataFrame link either not firing live or visually imperceptible"
created: 2026-05-22T00:00:00Z
updated: 2026-05-28T00:00:00Z
---

## Current Focus

hypothesis: CONFIRMED — perception failure due to layout (not a wiring bug). The cross-link DOES fire and DOES set proteinDf.selection correctly, but the visible effect (orange-highlighted points on volcano) lives on the protein TableView, while the user is looking at the enrichment TableView where the dot/bar charts dock. The protein view is backgrounded by addTableView(enrichmentDf) right before openEnrichmentVisualization.
test: Traced DataFrame identity through runEnrichmentPipeline → showEnrichmentDialog → openEnrichmentVisualization; verified what wireEnrichmentToVolcano mutates; confirmed dock-target view; verified Spectronaut Candidates parser sets SEMTYPE.GENE_SYMBOL.
expecting: Diagnosis confirmed. Reporting to planner.
next_action: Report ROOT CAUSE FOUND to caller; planner designs fix (likely: co-locate volcano on enrichment view OR dock enrichment charts onto protein view OR add visible feedback such as a tooltip/banner on the enrichment side).

## Symptoms

expected: After running Proteomics | Analyze | Enrichment on Spectronaut Candidates, selecting a term row in the merged enrichment DataFrame should highlight that term's proteins on the volcano of the protein DataFrame via enrichDf.onCurrentRowChanged subscription
actual: User reports charts feel isolated to the enrichment table; no perceptible cross-link to the protein table's volcano
errors: None reported. Unit test "cross-link regression" reportedly passes.
reproduction: Import Spectronaut Candidates → open volcano on protein DF → Proteomics | Analyze | Enrichment (WikiPathways on) → click term row in merged enrichment table → expect highlight on protein DF volcano
started: 2026-05-21 during Phase 13 UAT (first real-data exercise of 13-02 directional split + cross-link interaction)

## Eliminated

- hypothesis: append() returns a new DataFrame that orphans the wireEnrichmentToVolcano subscription
  evidence: src/analysis/enrichment.ts L356-366 reassigns `enrichmentDf = enrichmentDf.append(...)` BEFORE returning. The returned object is the post-merge DataFrame, and showEnrichmentDialog L458-459 passes the SAME `result.enrichmentDf` to both `addTableView` and `openEnrichmentVisualization`. Identity is preserved end-to-end. The cross-link wire IS bound to the right object.
  timestamp: 2026-05-28

- hypothesis: wireEnrichmentToVolcano sets only a tag or an invisible bitset (filter), so no rendering effect
  evidence: src/viewers/enrichment-viewers.ts L122-132 sets `proteinDf.selection` (DG BitSet) and calls `fireChanged()`. Datagrok scatter plots render selected rows as highlighted markers. The mutation is rendered — but on the protein DF's volcano viewer, which is not in the currently-focused view.
  timestamp: 2026-05-28

- hypothesis: Gene-symbol lookup fails on Spectronaut Candidates because the column isn't semantically tagged
  evidence: src/parsers/spectronaut-candidates-parser.ts L239 sets `geneCol.semType = SEMTYPE.GENE_SYMBOL`. findColumn(SEMTYPE.GENE_SYMBOL,...) at wireEnrichmentToVolcano L96 resolves it. Intersection column emits comma-space-joined gene names; the handler splits on `,` and trims (L119). Names match.
  timestamp: 2026-05-28

- hypothesis: Per-viewer subscription mismatch — user must click in dot/bar instead of in the table
  evidence: Partially mitigated, not the root cause. openEnrichmentVisualization L199-203 wires THREE subscriptions: merged enrichDf, upTop, downTop. Clicking in either the merged enrichment grid (in the enrichment TableView) or in the dot/bar charts (built on upTop/downTop) fires onCurrentRowChanged on the correct DataFrame and triggers the handler. The handler runs. So clicking is not the issue.
  timestamp: 2026-05-28

## Evidence

- timestamp: 2026-05-28
  checked: js-api DataFrame.append signature
  found: src/dataframe/data-frame.ts L412 — `append(t2, inPlace=false, columnsToAppend=null)` returns `new DataFrame(api.grok_DataFrame_Append(...))`. Default is NOT in-place; a fresh wrapper around a new dart handle is returned.
  implication: Any subscription bound to the pre-append DataFrame would be orphaned. The code knows this — it reassigns `enrichmentDf` to the append return value before any wiring is done. No orphan.

- timestamp: 2026-05-28
  checked: src/analysis/enrichment.ts L353-366 (Step 4 merge) and L450-459 (dialog onOK)
  found: `runEnrichmentPipeline` returns `{enrichmentDf, ...}` where `enrichmentDf` IS the merged frame (line 358 reassigns in loop, then line 359-365 re-apply name/tag/color). `showEnrichmentDialog` then does `grok.shell.addTableView(result.enrichmentDf); openEnrichmentVisualization(result.enrichmentDf, df);` — same object identity through to wiring.
  implication: Hypothesis 1 (orphaned subscription via append) is REFUTED. The cross-link object identity is correct.

- timestamp: 2026-05-28
  checked: src/viewers/enrichment-viewers.ts L162-217 (openEnrichmentVisualization)
  found: The function resolves `tv = grok.shell.tableViews where dataFrame===enrichDf` (already added in showEnrichmentDialog). All dot/bar charts are docked on `tv` (the enrichment table view) via `tv.dockManager.dock(...)`. The protein DataFrame and its volcano are NOT touched — they live on a separate TableView.
  implication: The four directional viewers and the enrichment grid are all on the ENRICHMENT view. The volcano is on the PROTEIN view. The user sees enrichment view after `addTableView(enrichmentDf)` switches focus.

- timestamp: 2026-05-28
  checked: src/viewers/enrichment-viewers.ts L92-134 (wireEnrichmentToVolcano)
  found: The subscription handler sets `proteinDf.selection` (clears all, sets matching rows, fireChanged). It does NOT touch the enrichment DF, does NOT add a viewer, does NOT show a notification or status text. The ONLY observable effect lives on `proteinDf` and its currently-bound viewers — i.e. the volcano on the protein TableView.
  implication: With the user focused on the enrichment view, the selection change happens off-screen. No co-located feedback exists in the enrichment view to indicate the cross-link fired.

- timestamp: 2026-05-28
  checked: src/package.ts L405-432 (top-menu handlers) and src/analysis/enrichment.ts L450-459 (dialog onOK)
  found: `enrichmentAnalysis` calls `showEnrichmentDialog(df)`. The dialog's onOK calls `grok.shell.addTableView(result.enrichmentDf)` — Datagrok's standard behavior is to switch focus to the newly added view. Then `openEnrichmentVisualization` docks viewers onto that new view. The protein TableView (with the volcano) is backgrounded.
  implication: The flow guarantees the user lands on the enrichment view post-pipeline. They cannot see the volcano without manually switching views — and they have no signal that switching would reveal the cross-link effect.

- timestamp: 2026-05-28
  checked: src/tests/enrichment.ts L195-213 ("wireEnrichmentToVolcano selects matching proteins")
  found: The test constructs synthetic `enrichDf` and `proteinDf` via `DataFrame.fromColumns`, calls `wireEnrichmentToVolcano`, sets `currentRowIdx = 0`, and asserts `proteinDf.selection.trueCount === 2`. The test does NOT exercise a merged (post-append) enrichDf, does NOT exercise the dock layout, and does NOT verify any user-visible rendering. It only verifies the selection-bitset contract.
  implication: The test correctly passes because the wire-up contract IS correctly implemented at the data layer. The unit test cannot detect the perception failure, which is a UI/UX layout issue, not a data-layer bug.

- timestamp: 2026-05-28
  checked: Spectronaut Candidates → DE-complete tag → enrichmentCharts handler (alternate entry, src/package.ts L412-432)
  found: This handler is `Proteomics | Visualize | Enrichment Charts...` — it's a re-open path for an already-existing enrichment DF. It re-runs `openEnrichmentVisualization(df, candidates[0])` with the user's current view being the enrichment DF. Same outcome: viewers dock on the enrichment view, the volcano remains on the protein view.
  implication: Even the re-open path inherits the same layout issue — the diagnosis applies to both Phase-9 (single direction) and Phase 13-02 (split) flows.

## Resolution

root_cause: The cross-DataFrame link IS firing correctly and IS mutating `proteinDf.selection` as designed. The user-perceived failure is a layout/visibility issue, not a wiring bug. (a) `showEnrichmentDialog` calls `grok.shell.addTableView(result.enrichmentDf)` (src/analysis/enrichment.ts L458), which switches focus to a new TableView dedicated to the enrichment DataFrame. (b) `openEnrichmentVisualization` (src/viewers/enrichment-viewers.ts L162-217) docks all four directional dot/bar charts — and the enrichment grid — onto THAT enrichment TableView. (c) `wireEnrichmentToVolcano` (L92-134) mutates only `proteinDf.selection` (with `fireChanged()`); the only visible effect of a row-click is highlighted markers on the volcano scatter plot — which lives on the SEPARATE protein TableView that the user is no longer looking at. (d) There is no co-located feedback on the enrichment view (no status text, no co-docked volcano, no indicator) to tell the user that something just happened on the protein view. Combined, this produces the verbatim symptom: "enrichment charts are associated with the enrichment table and not with the table that has the volcano plot." The user cannot perceive the cross-link because its effect renders off-screen, on a view that was backgrounded by the act of opening enrichment.

fix:
verification:
files_changed: []
