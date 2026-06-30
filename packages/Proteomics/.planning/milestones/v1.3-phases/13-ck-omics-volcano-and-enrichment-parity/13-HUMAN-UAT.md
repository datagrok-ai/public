---
status: complete
phase: 13-ck-omics-volcano-and-enrichment-parity
source: [13-VERIFICATION.md]
started: 2026-05-18T01:26:58Z
updated: 2026-06-04T00:00:00Z
rounds:
  - {date: 2026-05-28, outcome: "Test 1 issue (major) → Phase 14 follow-ups; Tests 2/3/4 issues → 13-07/13-08; Test 5 pass"}
  - {date: 2026-06-04, outcome: "Tests 2/3/4 re-run live against 13-09/13-10 + R4 fix → all pass; two new shipping bugs caught and fixed inline"}
---

## Current Test

[testing complete]

## Tests

### 1. Volcano visual parity vs the client CK-omics figure
expected: Importing the BP DMD/WT Spectronaut Candidates file and opening Proteomics | Visualize | Volcano Plot produces a figure whose point pattern/orientation matches `~/Downloads/ck/DMD_vs_WT/volcano_plots/` — with NO manual Comparison flip needed (the sign-normalization + declared-contrast default handle it).
result: issue
reported: "distribution matches. title does not match. axis titles could be more friendly. protein names for certain points are shown by default. up vs down is rephrased to 'enriched in WT' and 'enriched in DMD'. legend is not shown by default. search for proteins not obvious."
severity: major
gaps:
  - title: shows generic "Volcano" instead of contrast-aware "Volcano Plot: DMD vs WT"
  - axis_labels: shows raw column names "negLog10P"/"log2FC" instead of human "-Log10(Q-value)"/"Log2 Fold Change"
  - point_labels: top significant proteins not labeled by default (CK-omics labels ~15-20 top hits)
  - direction_semantics: legend reads generic "direction" with red/blue; CK-omics uses biology-aware "Enriched in WT" (magenta) / "Enriched in DMD" (cyan) with category counts
  - legend_visibility: no default legend; CK-omics shows a side legend with per-category counts (n=339/7864/126)
  - search_affordance: no visible protein search/highlight input; CK-omics has a "Search/Highlight Proteins" textbox below the plot
  - parity_pass: distribution shape and orientation DO match (no manual Comparison flip needed)

### 2. Directional enrichment layout + WikiPathways source
expected: Running Proteomics | Analyze | Enrichment with WikiPathways checked (default on) docks Up and Down dot+bar charts side-by-side (4 viewers), the merged table has a Direction column, WP terms appear, and selecting a term still highlights the linked proteins on the volcano (Phase-9 cross-link).
result: pass
note: "Screenshot confirms all 4 viewers docked (Up/Down × dot/bar), Direction column present, Source=WP for all rows (235 enrichment results), orange-highlighted protein on linked volcano confirms cross-link."

### 3. Live Volcano Options metric/color toggle stays synchronized
expected: Proteomics | Visualize | Volcano Options… → switching significance metric adj.p-value ↔ p-value moves the Y axis, the up/down/NS coloring, and the horizontal threshold line together (no stale disagreement); switching Color by significance ↔ location recolors by the locked 11-colour subcellular palette. p-value choice is absent when the input has no p-value column.
result: pass
adjacent_issues:
  - "Volcano Options dialog does not preload current values on second invocation (always shows defaults)."
  - "Switching Color → Location is slow; needs a progress indicator and/or faster pipeline."

### 4. Multi-contrast Candidates Filter viewer
expected: Importing a Candidates file with >1 distinct Comparison auto-docks a native Filters viewer scoped to the Comparison column (RIGHT); a single-comparison file docks no such filter.
result: pass
adjacent_issues:
  - "Filters viewer shows TWO columns instead of one: 'Flags' (Valid/significant) + 'Comparison (group1/group2)' joined by OR. Confirmed it's a single Filters viewer (not two). 'Flags' isn't created by our parsers — it's a Spectronaut Candidates input column, and the Filters viewer auto-includes it despite our `{columnNames: [cmpCol.name]}` scoping."

### 5. UniProt subcellular-location fetch + cross-session cache
expected: Coloring the volcano by Subcellular Location fetches real categories via the UniProt stream (fetchProxy), classifications look biologically sane for BP DMD/WT, and a second session loads the cached map quickly (userDataStorage, __schema_v) without re-fetching.
result: pass

## Summary

total: 5
passed: 4
issues: 1
pending: 0
skipped: 0
blocked: 0
adjacent_issues: 3
updated: 2026-05-28

## Gaps

- truth: "Volcano produces a figure whose point pattern/orientation matches CK-omics — with NO manual Comparison flip needed"
  status: partial
  reason: "Distribution and orientation match (no flip needed). But visual parity gaps: missing contrast-aware title, raw column names on axes instead of human labels, no default point labels for top hits, generic 'direction' coloring instead of 'Enriched in WT'/'Enriched in DMD' with magenta/cyan, no default legend with category counts, no protein search/highlight input."
  severity: major
  test: 1
  artifacts: []
  missing:
    - contrast-aware title in createVolcanoPlot (use proteomics.de_method + Numerator/Denominator from groups tag to render "Volcano Plot: <Num> vs <Den>")
    - axis label rewriting (negLog10P → "-Log10(Q-value)" or "-Log10(p-value)" depending on metric; log2FC → "Log2 Fold Change")
    - default label-top-N-by-Q overlay (CK-omics labels ~15-20 most significant)
    - direction category labels driven from group names (Enriched in <left-group>/<right-group>) with locked magenta/cyan palette + counts
    - legend visible by default with per-category n-counts
    - inline protein search/highlight input (or surface existing Volcano Options search affordance prominently)

- truth: "Volcano Options dialog reflects current viewer state when re-opened"
  status: failed
  reason: "User reported: re-opening Volcano Options shows defaults, not the values currently applied to the volcano (no state read-back from the live viewer/df tags)."
  severity: minor
  test: 3
  artifacts: []
  missing:
    - read current metric/color settings from viewer instance (or proteomics.* tags) and prefill ui.input.* fields in showVolcanoOptionsDialog
    - decide source of truth: viewer.getOptions() vs. df.tags vs. inputs' last-set values

- truth: "Multi-contrast Candidates dock a Filters viewer scoped to ONLY the Comparison column"
  status: failed
  reason: "User screenshot confirms a single Filters viewer is docked but it includes TWO filter columns: 'Flags' (Valid/significant from the Spectronaut Candidates input) AND 'Comparison (group1/group2)' joined by OR. Our intent (src/package.ts:128) was Comparison-only via {columnNames: [cmpCol.name]}, but the platform Filters viewer auto-appends the Flags column anyway."
  severity: minor
  test: 4
  artifacts:
    - "src/package.ts:117 dockComparisonFilterIfMultiContrast — passes {columnNames: [cmpCol.name]} but result includes Flags too"
    - "screenshot: Flags (16658/4930) | Comparison (DMD/WT, Treated/Control)"
  missing:
    - reproduce in isolation: build DG.Viewer.filters(df, {columnNames: ['Comparison (group1/group2)']}) on a candidates df and inspect getOptions() to see what columnNames the viewer actually adopts
    - if columnNames truly is ignored / merged with defaults: set the option post-create via viewer.setOptions({columnNames: [...]}), or build via DG.Viewer.fromType('Filters', df, {...}) with explicit `filters` array
    - if the viewer simply auto-includes low-cardinality string columns, look for a 'showHistogram'/'autoAddColumns'/'excludeColumns' style opt-out in the Filters viewer options

- truth: "Color → Location switch completes quickly and shows progress"
  status: failed
  reason: "User reported: switching Color by Location is slow with no visible progress indicator. Likely re-fetching/recomputing the UniProt subcellular map for visible rows."
  severity: minor
  test: 3
  artifacts: []
  missing:
    - wrap the location-color path in a progress indicator (DG.TaskBarProgressIndicator or ui.setUpdateIndicator on the viewer host)
    - profile + cache hot path: confirm we hit the cross-session cache (userDataStorage __schema_v) instead of refetching; batch+memoize per-protein lookups; render incrementally so first paint isn't blocked by full-table classification

## Round 3 Re-run — 2026-06-04 (post-13-09 / 13-10)

Driven live in Chrome against localhost Datagrok with Playwright. Fixtures:
`~/Downloads/ck/2026-03-25_BP_EMT_DMD_Candidates.tsv` (8329 rows, single
contrast) and `~/Downloads/ck/2026-05-21_BP_DMD_WT_two_comparisons.tsv`
(16658 rows, two contrasts).

### R2 — Enrichment chart dock target (Test 2 re-run, post-13-09)
result: pass (live)
observed:
  - openEnrichmentVisualization docked 4 dot/bar viewers on the protein TableView (the one with the volcano), not on a new enrichment view
  - countViewersBoundTo(enrichTv, proteinDf) == 0 — 13-07 'Volcano (linked)' co-dock fully gone
  - grok.shell.v ended on the protein view after the dialog closed
  - Cross-link confirmed end-to-end: setting enrichDf.currentRowIdx=0 (Intersection="SLC34A2, ADGRF5") highlighted rows 3193 + 5977 on proteinDf, "2 selected rows" surfaced in the context panel
shipping_bugs_found_and_fixed_inline:
  - "Initial 13-09 implementation switched focus to proteinTv AFTER docking. addTableView(enrichDf) had backgrounded the protein view, so dock-manager splitters were 0-dimensional, the Dart scatter plot's smart-labels feature crashed with 'drawImage on canvas of width or height 0', and every docked viewer rendered at w=0/h=0. Fix: move grok.shell.v = proteinTv to BEFORE the dock storm (fix(13-09): switch focus to protein view BEFORE docking)."
artifacts:
  - uat-13-04 / uat-13-07 / uat-13-08 (.playwright-mcp/) screenshots

### R1 / R6 — In-volcano busy overlay (Test 3 re-run, post-13-10)
result: pass (live)
observed:
  - userDataStorage('proteomics-subcell-loc') cleared to force first-time fetch
  - Overlay attached at t≈531ms with "Classifying subcellular locations…"
  - Phase rotated to "Fetching subcellular locations 1/84 (1%)" at t≈1288ms; ticked per chunk (2/84 at t=1741, 3/84 at t=2798)
  - Both fetch-acc (84 chunks) and fetch-gene (40 chunks) phases surfaced live; total wall time ≈ 88s
  - Overlay anchored on the volcano viewer (≈113×95 px, centered)
  - Volcano re-rendered colored by locked 11-category subcellular palette after fetch completion
warm_cache_no_strobe:
  - Sampled DOM at every requestAnimationFrame for 2s during warm-cache toggle (location → significance → location): 219 frames, 0 frames with overlay present. Recompute returns within a microtask on the LOCATION_HASH_TAG short-circuit so attach + detach happens entirely between paint cycles. Better than the planned <1-frame flash.
artifacts:
  - uat-13-09 (.playwright-mcp/) — overlay mid-fetch screenshot
  - uat-13-10 (.playwright-mcp/) — final colored volcano

### R4 — Multi-contrast Filters auto-dock (Test 4 re-run)
result: pass (live, after inline fix)
observed:
  - Two-contrast Candidates import auto-docked a Filters viewer (right, ratio 0.3)
  - After fix: three categorical filters render — Comparison (DMD/WT: 8329, Treated/Control: 8329), Display Name (gene-symbol checklist with row counts), Source ID (LOC IDs + search box for the 16654-row long-tail)
shipping_bugs_found_and_fixed_inline:
  - "dockComparisonFilterIfMultiContrast passed `filters: filterSpecs` (G4 typed shape from 14-03). The Datagrok options serializer dropped that key entirely — getOptions().look returned only showBoolCombinedFilter, with NO filters or columnNames, so the viewer docked empty in production. Static test (multi-comparison docks a Filters viewer) was satisfied because it only asserts viewer presence, not column wiring. Fix: swap typed `filters: filterSpecs` for `columnNames: string[]` (the stable shape that round-trips). The G4 Flags-exclusion verification gate became redundant and was removed — columnNames doesn't auto-include Flags (fix(13): Filters viewer config uses columnNames, not typed filters array)."
artifacts:
  - uat-13-11 (.playwright-mcp/) — broken state before fix (filters viewer collapsed/empty)
  - uat-13-12 (.playwright-mcp/) — post-fix verification with columnNames (rendered correctly)
  - uat-13-13 (.playwright-mcp/) — re-import via the fixed handler

## Updated Summary (post-Round 3)

total_rounds: 3
tests_now_passing: 5 of 5 (Tests 2/3/4 confirmed live; Tests 1/5 unchanged from Round 2)
new_shipping_bugs_caught_in_round_3: 2 (13-09 focus order, R4 Filters columnNames vs. typed filters)
follow_up_phase: Phase 14 absorbed the Test-1 gaps (title, axis labels, top-N labels, magenta/cyan direction, legend, search affordance)
