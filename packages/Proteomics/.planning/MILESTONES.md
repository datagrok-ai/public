# Milestones

## v1.3 Richer Dialogs & UX (Shipped: 2026-06-05)

**Phases completed:** 5 phases (10–14) + 2 single-commit hotfixes (999.1, 999.4)
**Plans:** 25 (Phases 10-11: 6, Phase 12: 4, Phase 13: 10 incl. 4 gap-closure, Phase 14: 5, hotfixes: 2)
**Timeline:** 2026-03-07 → 2026-06-05 (initial Phases 10-11 cut shipped 2026-03-10; extended through May-June with Phases 12/13/14 + two hotfixes responding to Cytokinetics engagement findings)
**Audit:** `milestones/v1.3-MILESTONE-AUDIT.md` — verdict PASSED (17/17 requirements, 14/14 human-UAT tests across 3 phases, integration verified end-to-end)
**Known deferred items at close:** 12 (see STATE.md Deferred Items — all diagnosed/perception-already-addressed/backlog enhancements, not blockers)

**Key accomplishments:**

*Phases 10-11 — original v1.3 cut (March):*
1. Spectronaut long-format TSV parser with auto-detection, pivot, CON__/REV__ filtering, and auto-group annotation
2. Quantile normalization (client-side), VSN normalization (R script with quantile fallback), kNN imputation with progress indicator
3. Multi-method normalization dialog with reactive box plot preview and Spectronaut pre-normalized warning
4. Multi-method imputation dialog with conditional parameters and valid-values protein filter
5. DE dialog with comparison direction picker, t-test method, and dynamic FC hint text
6. Descriptive viewer titles on volcano/PCA/heatmap and filename-based DataFrame naming

*Phase 12 — Spectronaut Input Coverage (May):*
7. `parseSpectronautStream` streams precursor-level Spectronaut TSV via `Blob.stream()` + `TextDecoderStream`, single-pass aggregating to the same wide protein×sample DataFrame with duckdb-parity filters/aggregates, bounded memory, bytes-read progress bar + explicit macrotask yield — wired into `importSpectronaut` behind a header sniff while the proven PG-level `file.text()` path stays byte-identical
8. Committed duckdb aggregation oracle + synthetic precursor fixture + verbatim JSON sidecar + extended `makeLongFormatTsv`; 7 new `grok test` cases lock the streaming path against the fixture/golden/sidecar (per-cell numeric equivalence within 1e-3 to both the text path and the duckdb golden)
9. `handleFields` returns a discriminated `'kept' | 'malformed' | 'filtered'` outcome so the streaming import surfaces "malformed line(s)" ONLY for genuinely unparseable rows

*Phase 13 — CK-omics Volcano + Enrichment Parity (June):*
10. Added `SEMTYPE.SUBCELLULAR_LOCATION` + mirrored detector; shared `src/analysis/subcellular-location.ts` with verbatim CK-omics 11-category classifier + locked hex palette, cached chunked fetchProxy stream fetch with reviewed-by-gene fallback
11. `runEnrichmentPipeline` runs one g:Profiler query per direction (up/down) over a shared all-detected background, merges into one Direction-tagged DataFrame, renders Up/Down dot+bar side-by-side, adds WikiPathways as a default-on source — Phase-9 cross-link untouched
12. Candidates rows whose declared comparison is the exact reverse of canonical orientation are flipped per-row (log2FC negated, AVG Group Quantity + Condition swapped, Comparison relabeled) in a pure parser; multi-contrast files dock a native Comparison Filters viewer
13. Volcano refactor: metric-parameterized Y/direction helpers, locked-palette Subcellular Location color column, `recomputeVolcano` that moves Y + up/down/NS + threshold lines + color together — wired to a `Volcano Options...` dialog with adj.p-value/significance defaults
14. Wave-5 gap closures: 13-09 dock enrichment dot/bar on protein TableView (Option A), 13-10 in-volcano busy overlay for the first subcellular-location fetch

*Phase 14 — CK-omics Analyst-Experience Enhancements (June):*
15. Gene-label resolver foundation: 4 new SEMTYPEs + detectors, `src/utils/gene-label-resolver.ts` verbatim port from CK-omics, wired into all 5 parsers
16. Volcano polish: D-04 magenta/cyan/gray + group-name legend, D-03 top-15 labels, D-06 live counter overlay, title + axis labels, dialog state preload, metric-aware progress wording
17. Unified Filters viewer scoping (no Flags leak), Display Name + Source ID free-text search to `df.selection` union with top-N
18. UniProt panel per-group magenta/cyan bars + Group-Mean Correlation viewer (Numerator Mean / Denominator Mean + inline Pearson/Spearman + y=x diagonal)
19. Smart pathway filter port — verbatim CK-omics `apply_smart_pathway_filtering`: enrichment dialog checkbox, pipeline insertion before `buildEnrichmentDf`, banner above grid

*Hotfixes (June):*
20. 999.1: Annotate Experiment dialog pre-fills from existing `proteomics.groups` (closes Phase 12 UAT obs-A data-loss risk — auto-grouped Spectronaut DMD/WT would otherwise be silently overwritten by an empty OK)
21. 999.4: Spectronaut streaming↔text malformed-counter parity (closes 12-REVIEW.md WR-01/WR-03 — empty-protein + both-null-casts now `'filtered'` silent matching the text path; only structurally truncated rows surface as "malformed")

---

## v1.2 Biological Interpretation (Shipped: 2026-03-07)

**Phases completed:** 2 phases, 4 plans, 12 commits
**Lines added:** 938 TypeScript/JS
**Timeline:** 1 day (2026-03-07)
**Git range:** 046404a1ee..9c3d587880

**Key accomplishments:**

1. g:Profiler API integration for gene ID mapping (UniProt to gene symbols) with 9 organism support
2. GO/KEGG/Reactome overrepresentation analysis with proteomics-specific background set and FDR correction
3. Enrichment dot plot and bar chart visualization with top-N filtering
4. Cross-DataFrame selection wiring: select enrichment term -> highlight member proteins on volcano plot
5. 15 unit tests (7 enrichment + 8 enrichment visualization)

---

## v1.1 Parser & QC (Shipped: 2026-03-07)

**Phases completed:** 2 phases, 5 plans, 16 commits
**Lines added:** 1,725 TypeScript/JS
**Timeline:** 2 days (2026-03-06 to 2026-03-07)
**Git range:** d2371aedb1..7bcdd9c764

**Key accomplishments:**

1. Generic CSV/TSV matrix import with column mapping dialog, auto-suggestion, log2 detection, and live preview
2. QC dashboard with 7 linked viewers (MA plot, trend line, CV, correlation, missing heatmap, missing bar, box plot) in tiled layout
3. Shared parser utilities extracted for reuse across MaxQuant and generic parsers
4. 18 unit tests (11 generic parser + 7 QC computation)
5. BigInt column support and DG.Viewer factory pattern for auxiliary DataFrames

---

## v1.0 MVP (Shipped: 2026-03-05)

**Phases completed:** 5 phases, 12 plans, 36 commits
**Lines of code:** 2,530 TypeScript/JS
**Timeline:** 6 days (2026-02-28 to 2026-03-05)
**Git range:** f89fd20..6a0ec4e

**Key accomplishments:**

1. MaxQuant proteinGroups.txt parser with auto-detection, contaminant filtering, log2 transform, and semantic types
2. Full analysis pipeline: group annotation, median normalization, MinProb imputation, Welch's t-test + BH FDR
3. Interactive visualizations: volcano plot, expression heatmap with clustering, PCA colored by group
4. UniProt info panel with auto-trigger on protein ID columns
5. DEqMS/limma R scripts with three-level fallback chain (DEqMS -> limma -> client-side t-test)
6. Phase 5 hardening: detector regex fixes, heatmap filter isolation, progress indicators

### Known Gaps

- Heatmap dendrogram not rendering (UAT Test 5) — tree rendering code commented out, dFunc.apply passes wrong argument types
- IMPORT-06 (Demo in Datagrok Demo app) — deferred to v2 backlog

---
