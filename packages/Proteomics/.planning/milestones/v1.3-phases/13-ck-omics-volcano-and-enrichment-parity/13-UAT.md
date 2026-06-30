---
status: diagnosed
phase: 13-ck-omics-volcano-and-enrichment-parity
source: [13-01-SUMMARY.md, 13-02-SUMMARY.md, 13-03-SUMMARY.md, 13-04-SUMMARY.md, 13-05-SUMMARY.md, 13-06-SUMMARY.md]
started: 2026-05-19T13:39:00Z
updated: 2026-05-28T13:50:00Z
---

## Current Test

[testing complete]

## Tests

### 1. Candidates volcano visual parity (sign-normalization)
expected: Import the BP DMD/WT Spectronaut Candidates file, open Proteomics | Visualize | Volcano Plot. Point pattern/orientation matches ~/Downloads/ck/DMD_vs_WT/volcano_plots/ with NO manual Comparison flip needed (13-03 per-row sign normalization). Candidates is pre-computed DE — Proteomics | Analyze | DE correctly shows the Candidates guard message (expected behavior, not an issue).
result: pass

### 2. Directional enrichment layout + WikiPathways source
expected: Proteomics | Analyze | Enrichment with WikiPathways checked (default on) docks Up and Down dot+bar charts side-by-side (4 viewers); merged table has a Direction column; WP terms appear; selecting a term still highlights the linked proteins on the volcano (Phase-9 cross-link preserved).
result: issue
reported: "enrichment charts are associated with the enrichment table and not with the table that has to volcano plot associated with it"
severity: major

### 3. Live Volcano Options metric/color toggle stays synchronized
expected: Proteomics | Visualize | Volcano Options… → switching significance metric adj.p-value ↔ p-value moves the Y axis, the up/down/NS coloring, and the horizontal threshold line together (no stale disagreement). Switching Color by significance ↔ location recolors by the locked 11-colour subcellular palette. The p-value choice is absent when the input has no p-value column.
result: issue
reported: "changing options kicks off some calculations that never return" — later amended: "previous test - eventually returned" (screenshot showed locked 11-colour subcellular palette, Y-axis negLog10P, threshold line at −log10(0.05) = 1.301, direction column populated — i.e. spec functionally met but first-time Color-by-location fetch is slow enough to feel like a hang; no progress UI)
severity: minor

### 4. Multi-contrast Candidates Filter viewer
expected: Importing a Candidates file with >1 distinct Comparison auto-docks a native Filters viewer scoped to the Comparison column (docked RIGHT). A single-comparison file docks no such filter.
result: pass
note: First synthesized fixture had identical log2FC in both blocks (filter propagation ambiguous); re-synthesized with +1.5 offset on block B's AVG Log2 Ratio so filter propagation is visually distinguishable. Filter docks scoped to Comparison column (built-in Flags row-state section is intrinsic to DG.Viewer.filters and unrelated to column scoping).

### 5. UniProt subcellular-location fetch + cross-session cache
expected: Coloring the volcano by Subcellular Location fetches real categories via the UniProt stream (fetchProxy); classifications look biologically sane for BP DMD/WT; a second session loads the cached map quickly (userDataStorage, __schema_v) without re-fetching.
result: pass

### 6. Report-path DE defaults to declared contrast (D-09)
expected: Import the matching Spectronaut **Report** (PG report — per-sample intensities, via Proteomics | Import | Spectronaut Report), run the Annotate → (Normalize/Impute) → Proteomics | Analyze | Differential Expression pipeline. The DE dialog's Comparison dropdown defaults to the declared contrast (numerator = group1 = Spectronaut's declared Numerator), so the resulting volcano has the correct orientation WITHOUT manually flipping the dropdown.
result: pass

## Summary

total: 6
passed: 4
issues: 2
pending: 0
skipped: 0
blocked: 0

## Gaps

- truth: "Running directional enrichment docks Up/Down dot+bar charts and selecting a term highlights the linked proteins on the volcano (Phase-9 cross-DataFrame link)"
  status: failed
  reason: "User reported: enrichment charts are associated with the enrichment table and not with the table that has to volcano plot associated with it"
  severity: major
  test: 2
  root_cause: "The cross-DF link IS firing correctly — wireEnrichmentToVolcano mutates proteinDf.selection (clearAll → set matching rows → fireChanged) on every click. The bug is layout/visibility: showEnrichmentDialog calls grok.shell.addTableView(enrichmentDf) which switches focus to a new TableView, and openEnrichmentVisualization then docks all 4 directional viewers + the merged grid onto that enrichment TableView. The volcano remains on the now-backgrounded protein TableView (never touched), so the highlight renders on a view the user is not looking at. Three subscriptions are wired (merged + upTop + downTop) so click targets are correct; the existing test asserts only the data-layer contract (selection.trueCount === 2) on synthetic frames and cannot catch the post-addTableView focus-switch."
  artifacts:
    - path: "src/viewers/enrichment-viewers.ts"
      issue: "openEnrichmentVisualization (L162-217) docks all 4 directional viewers + grid onto enrichmentDf's TableView; wireEnrichmentToVolcano (L92-134) updates only proteinDf.selection with no co-located feedback on the enrichment view"
    - path: "src/analysis/enrichment.ts"
      issue: "showEnrichmentDialog onOK calls grok.shell.addTableView(result.enrichmentDf) (L458-459), which switches focus away from the protein view before openEnrichmentVisualization runs"
    - path: "src/package.ts"
      issue: "enrichmentAnalysis (L405-410) and enrichmentCharts (L412-432) both terminate in openEnrichmentVisualization(enrichDf, proteinDf) — both inherit the same layout pattern"
    - path: "src/tests/enrichment.ts"
      issue: "L195-213 cross-link test asserts data-layer contract only on synthetic frames; never exercises post-addTableView dock layout or user-visible rendering"
  missing:
    - "UX layout decision (planner picks one): (A) dock the Up/Down dot+bar viewers onto the *protein* TableView (alongside the volcano) instead of the enrichment view; (B) keep them on the enrichment view but also co-dock a viewer of proteinDf (volcano or a selection-filtered gene grid) onto the enrichment view; (C) keep current layout and add explicit co-located feedback (status line / toast / auto-bring-protein-view-to-front on row-click)"
    - "Integration/layout test that asserts the cross-link is visible to the user post-addTableView (not just the synthetic-frame selection.trueCount check)"
  debug_session: .planning/debug/enrichment-cross-link-not-firing.md

- truth: "Toggling Volcano Options metric or color-by recomputes the volcano with timely user feedback — Y axis, up/down/NS coloring, threshold lines, and colorColumnName all update together, and the user is not left wondering whether the app has hung"
  status: failed
  reason: "User initially reported: 'changing options kicks off some calculations that never return'; later amended to 'previous test - eventually returned' after a screenshot confirmed the spec is functionally met (locked 11-colour subcellular palette, Y bound to negLog10P, threshold line at 1.301, direction column populated). Root concern is UX: first-time Color-by-location triggers ~80 chunked UniProt fetchProxy calls (13-04) with no progress UI, which feels like a hang."
  severity: minor
  test: 3
  root_cause: "Three compounding conditions: (R1 PRIMARY) No progress signal during fetch — volcanoOptions creates DG.TaskBarProgressIndicator but pi.update(...) is NEVER called in the call chain; getSubcellularLocations accepts no progress callback; no grok.shell.info, no chunk counter, no percent — static indeterminate label is the only feedback. (R2 PRIMARY) Both UniProt fetch passes in getSubcellularLocations are strictly sequential `for…of chunk { await fetchProxy }` with no Promise.all and no concurrency — for an 8k-accession Spectronaut Candidates file at ACC_CHUNK=100 that's 80 sequential UniProt stream round-trips plus the GENE_CHUNK=20 Pass-2 fallback, so 80–300+ s of opaque waiting. (R3 SECONDARY) ensureLocationColumn never short-circuits when the column is already populated — every toggle re-walks Primary Protein ID, re-calls getSubcellularLocations, re-bulk-inits the column, even when the accession set hasn't changed. Also: userDataStorage.put is atomic-at-end of both passes — interrupted sessions lose every fetched accession."
  artifacts:
    - path: "src/package.ts"
      issue: "volcanoOptions L317-324 creates pi = DG.TaskBarProgressIndicator.create('Updating volcano...') and closes in finally but never calls pi.update(); dialog dismisses immediately on OK with no replacement spinner, leaving the volcano frozen in old state for the full fetch duration"
    - path: "src/analysis/subcellular-location.ts"
      issue: "Pass 1 (L227-247) and Pass 2 (L260-286) use `for (const group of chunk(...)) { await grok.dapi.fetchProxy(url) }` — strictly sequential, no Promise.all, no concurrency cap, no progress-callback parameter; userDataStorage.put cache write at L296-303 is atomic-at-end so interrupted sessions lose all work"
    - path: "src/viewers/volcano.ts"
      issue: "ensureLocationColumn L87-116 unconditionally re-runs accession scan + getSubcellularLocations + bulk-init on every invocation; column-reuse on L111 preserves scatter binding but does NOT short-circuit the underlying work"
  missing:
    - "Thread a `(done, total, phase)` progress callback through getSubcellularLocations → ensureLocationColumn → volcanoOptions so the existing pi can call pi.update(percent, 'Fetching subcellular locations: N/80')"
    - "Pre-OK toast (e.g. grok.shell.info('Fetching subcellular locations for first-time use; subsequent toggles will be cached')) so users know what they signed up for"
    - "Bounded Promise.all concurrency in getSubcellularLocations fetch loops (e.g. 4–8 concurrent chunks) to cut wall-clock 4–8× without overwhelming rest.uniprot.org"
    - "Incremental userDataStorage cache writes (every N chunks) so an interrupted session retains partial progress"
    - "Optional: accession-set hash as a column tag on Subcellular Location so a repeat toggle with the same accession set short-circuits the entire pipeline (handles R3)"
  debug_session: .planning/debug/volcano-options-first-fetch-no-progress.md
