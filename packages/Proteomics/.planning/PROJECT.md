# Proteomics Package for Datagrok

## What This Is

A Datagrok extension package for mass spectrometry-based proteomics data analysis. Scientists can import proteomics data from MaxQuant, Spectronaut, or any CSV/TSV matrix, run a complete analysis pipeline with multi-method normalization (median/quantile/VSN), imputation (MinProb/kNN/zero/mean/median), and differential expression (DEqMS/limma/t-test) through polished dialogs with visual feedback, assess data quality through a linked QC dashboard, perform gene ontology and pathway enrichment analysis, and visualize results with interactive volcano plots, enrichment charts, heatmaps, and PCA -- all within Datagrok.

## Core Value

Scientists can import proteomics data and go from raw protein quantification to differential expression results with biological interpretation -- all within Datagrok, with no tool-switching or file exports.

## Requirements

### Validated

- ✓ MaxQuant proteinGroups.txt import with auto-detection, filtering, log2 transform, semantic types -- v1.0
- ✓ Bundled demo dataset -- v1.0
- ✓ Experiment annotation dialog with group persistence -- v1.0
- ✓ Median centering normalization (client-side TypeScript) -- v1.0
- ✓ MinProb imputation with configurable parameters -- v1.0
- ✓ Welch's t-test + BH FDR differential expression (client-side TypeScript) -- v1.0
- ✓ DEqMS/limma R scripts with three-level fallback chain -- v1.0
- ✓ DE result columns with correct semantic types -- v1.0
- ✓ Volcano plot with threshold lines and direction coloring -- v1.0
- ✓ Expression heatmap with hierarchical clustering -- v1.0
- ✓ PCA plot colored by experimental group -- v1.0
- ✓ UniProt info panel with auto-trigger on ProteinId columns -- v1.0
- ✓ Generic CSV/TSV matrix import with column mapping dialog and auto-suggestion -- v1.1
- ✓ QC dashboard with 7 linked viewers (MA, CV, correlation, missing values, intensity distributions) -- v1.1
- ✓ Gene ID mapping (UniProt to gene symbols via g:Profiler, 9 organisms) -- v1.2
- ✓ GO/KEGG/Reactome overrepresentation analysis with proteomics background -- v1.2
- ✓ Enrichment visualization (dot plot, bar chart, interactive table) -- v1.2
- ✓ Cross-DataFrame volcano integration (select enrichment term -> highlight proteins) -- v1.2
- ✓ Spectronaut long-format TSV import with auto-detection, pivot, filtering, auto-groups -- v1.3
- ✓ Quantile normalization (client-side TypeScript) -- v1.3
- ✓ VSN normalization (R script with quantile fallback) -- v1.3
- ✓ kNN imputation with progress indicator -- v1.3
- ✓ Zero, mean, median imputation methods -- v1.3
- ✓ Multi-method normalization dialog with before/after distribution plots and Spectronaut pre-normalized warning -- v1.3
- ✓ Multi-method imputation dialog with conditional parameters and valid-values filter -- v1.3
- ✓ DE dialog with comparison direction picker and method-specific parameters -- v1.3
- ✓ Descriptive viewer titles on volcano, PCA, and heatmap -- v1.3
- ✓ DataFrame retains imported filename as name -- v1.3
- ✓ Gene-label resolution at parse time (ENSRNOG/ENSMUSG/ENSG/ENSDARG/MGP/LOC/RGD/AABR → Display Name + Source ID with `*`/`†` provenance) -- v1.3 Phase 14
- ✓ Live filter-aware counter overlay on the volcano (Visible Proteins recompute on filter/selection/property change) -- v1.3 Phase 14
- ✓ UniProt panel per-group quantities (SVG bar chart, mean ± SD, magenta/cyan per D-04) -- v1.3 Phase 14
- ✓ Group-Mean Correlation viewer (Numerator Mean vs Denominator Mean, inline Pearson r + Spearman ρ, y=x diagonal) -- v1.3 Phase 14
- ✓ Smart hierarchical pathway filtering (generic-parent drop, per-source cap, default-on with dialog opt-out) -- v1.3 Phase 14
- ✓ Filters viewer G4 root-cause fix + D-05 unified protein search (typed filter spec, Display Name / Source ID free-text, df.selection highlight-not-hide) -- v1.3 Phase 14
- ✓ Volcano polish: D-04 magenta/cyan/gray + group-name legend, synthesized title + axis labels, default top-15 labels, metric-aware progress wording -- v1.3 Phase 14

### Active

**Milestone v1.4: Cross-Team Review** -- Cytokinetics-driven cross-team workflow milestone.
Requirements TBD (defined by REQUIREMENTS.md after research phase).

## Current Milestone: v1.4 Cross-Team Review

**Goal:** Enable proteomics experts to publish frozen analyses for biologist consumers, run longitudinal screening campaigns with weekly QC trust, and compare compound impact across runs -- all within Datagrok, no tool-switching.

**Target features:**
- Read-only analysis publishing as a Datagrok Project (target-keyed, trimmed columns -- protein ID / gene / log2FC / p / adj.p / sig, optional enrichment -- no raw intensities), with reviewer-group permissions
- Sample-level SPC tracking across weekly assays (intensity median, missingness %, control-correlation, protein count) with Shewhart/Nelson drift rules and run-level pass/flag tagging
- Screening-campaign comparison -- first-class run/compound data model + side-by-side volcano comparison viewer with cross-run highlighting + Chem package integration for structures

**Key context:**
- Cytokinetics-driven; next check-in is the milestone close
- Tight scope: defers perf debt (limma cold-start, heatmap progress), March-2026 dialog polish cluster, parser expansion (FragPipe MaxLFQ, DIA-NN), CK-omics cluster C, and 999.2 banner-wording fix
- Phase numbering continues from v1.3 (starts at Phase 15)

### Future

- Additional parsers: FragPipe, DIA-NN, Proteome Discoverer
- Expanded import sources: local drive, Datagrok folders, existing DataFrames
- UniProt batch annotation
- GSEA (rank-based gene set enrichment, beyond ORA)
- Functional annotation clustering
- STRING PPI network viewer

### Out of Scope

- Raw instrument file processing (.raw, .mzML) -- search engines handle this
- De novo peptide sequencing -- specialized tools exist
- Custom volcano JsViewer -- Datagrok's scatter plot handles this natively with formula lines
- Single-cell proteomics workflows -- future
- Spatial proteomics visualization -- future
- Multi-omics integration -- future
- Instrument control or LIMS -- out of scope for analysis package
- Pathway diagram visualization (colored KEGG maps) -- licensing issues
- GO DAG visualization -- large DAG rendering slow; rarely useful in practice
- Spectronaut peptide-level quantification -- doubles parser complexity for marginal gain
- Tail-Robust Quantile Normalization -- niche improvement over standard quantile
- Batch effect correction (ComBat) -- different statistical procedure requiring batch metadata
- Multi-group DE (ANOVA-style) -- fundamentally different statistical tests

## Context

Shipped v1.3 closed 2026-06-05 with ~12,985 LOC (7,607 src + 4,558 tests + 152 detectors.js + 121 R scripts + 547 tools/SQL).
Tech stack: TypeScript, webpack, Datagrok API, R (limma/DEqMS/VSN via Bioconductor), g:Profiler REST API, UniProt REST API.
Package scaffold at `packages/Proteomics/` with full build pipeline.
v1.3 milestone spanned 2026-02-28 → 2026-06-05; ~270 commits across the package; Phases 10-11 in March, Phases 12-14 + hotfixes May–June driven by Cytokinetics engagement.

Known issues / carry-forward to next milestone:
- Heatmap dendrogram not rendering (tree code commented out, dFunc.apply argument mismatch) -- carried from v1.0
- limma-de-slow / showheatmap-hangs -- performance debt from 2026-03 (`.planning/debug/`, status: diagnosed)
- Phase 12 + 14 `VALIDATION.md` drafts (`nyquist_compliant: false`) -- documentation hygiene
- 20 pending enhancement todos in `.planning/todos/pending/` (mostly March-2026 dialog/title polish ideas; backlog)
- Backlog phases 999.2 (banner wording — cosmetic) + 999.3 (publish read-only analysis for biologist review — strategic v1.5 spine)

## Constraints

- **Tech stack**: Must follow Datagrok package conventions (TypeScript, webpack, grok CLI tools, function metadata comments)
- **Platform**: All analysis functions run either client-side (TypeScript) or server-side (R scripts via Datagrok compute infrastructure)
- **Webpack externals**: Must not bundle datagrok-api, rxjs, cash-dom, dayjs, openchemlib, wu
- **Viewers**: Prefer leveraging existing Datagrok viewers (scatter plot, heatmap, PCA) over building custom JsViewers where possible
- **Demo data**: Use publicly available MaxQuant datasets rather than proprietary data

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Use Datagrok scatter plot for volcano instead of custom JsViewer | Platform already provides interactive scatter with formula lines, linked selection, and hover labels | ✓ Good -- works well, no custom viewer needed |
| R for differential expression (not TypeScript) | limma/DEqMS are the gold standard; reimplementing moderated t-statistics in TS would be error-prone | ✓ Good -- three-level fallback chain handles R unavailability gracefully |
| Client-side TS for normalization and imputation | Median centering and MinProb are simple enough to run in-browser with instant feedback | ✓ Good -- instant feedback, no server dependency |
| MaxQuant first, other parsers later | MaxQuant is most cited search engine; proves the pattern for other parsers | ✓ Good -- pattern established, generic and Spectronaut parsers now extend it |
| UniProt info panel (not batch annotation) | Lighter weight, immediate value on click, avoids slow bulk API calls on import | ✓ Good -- responsive UX |
| DataFrame clone for heatmap filter isolation | Heatmap needs top-N filter without affecting shared DataFrame used by volcano | ✓ Good -- fixed cross-viewer filter leak |
| Synthetic demo data over real public dataset | Controlled edge case coverage for tests | ✓ Good -- predictable test assertions |
| PCA in separate table view | Sample-level data has different row count than protein-level DataFrame | ✓ Good -- clean separation |
| Shared parser utilities (v1.1) | Extract reusable functions from MaxQuant parser for generic parser | ✓ Good -- clean separation of format-specific vs generic logic |
| DG.Viewer factory for auxiliary DataFrames (v1.1) | Avoid grok.shell.addTable() for DataFrames used only as viewer data sources | ✓ Good -- fixed circular JSON error in QC dashboard |
| Moving-average for loess approximation (v1.1) | Simple sliding window sufficient for QC bias detection | ✓ Good -- fast, no external dependency |
| ScatterPlot as dot plot (v1.2) | Reuse existing viewer via sizeColumnName/colorColumnName instead of custom rendering | ✓ Good -- zero custom viewer code |
| g:Profiler for enrichment (v1.2) | Well-maintained REST API, supports custom background, multiple organisms | ✓ Good -- no local database needed |
| Cross-DF selection via onCurrentRowChanged (v1.2) | Subscribe to enrichment table row changes, parse Intersection column, set selection on protein DF | ✓ Good -- clean integration without tight coupling |
| Two-phase v1.3 structure: algorithms first, dialogs second (v1.3) | Algorithms (Phase 10) provide foundation; dialogs (Phase 11) consume them | ✓ Good -- clean dependency chain |
| Zero new npm dependencies for v1.3 (v1.3) | Quantile norm and kNN in pure TypeScript, VSN via R script | ✓ Good -- no bundle size increase |
| Spectronaut pivot via Map<protein, Map<sample, ibaq>> (v1.3) | First-encountered-value deduplication handles PG.IBAQ being constant per protein+sample | ✓ Good -- simple and correct |
| Unconditional preNormalized tag for Spectronaut (v1.3) | Spectronaut always normalizes regardless of export format | ✓ Good -- downstream warning always fires |
| Viewer title passthrough pattern (v1.3) | Viewer functions accept optional title, callers pass analysis context | ✓ Good -- backward-compatible, flexible |
| Parser no-name policy (v1.3) | Parsers return unnamed DataFrames; import handlers set df.name from filename | ✓ Good -- single responsibility |
| Header-sniff routing for Spectronaut precursor vs PG-level (v1.3 Phase 12) | First-line column-name check picks streaming-vs-text path; PG-level remains byte-identical | ✓ Good -- additive, no regression |
| Streaming Spectronaut: `Blob.stream()` + `TextDecoderStream` + duckdb-parity bounded-aggregation Map (v1.3 Phase 12) | Avoids V8 512MB string ceiling; single-pass aggregation; byte-bounded memory | ✓ Good -- proven on 2.6GB reference file |
| `LineOutcome` discriminated classification (v1.3 Phase 12 + 999.4 hotfix) | Splits structurally-malformed from by-design-filtered drops; streaming↔text user-message parity | ✓ Good -- closes false corruption signal |
| `seedAnnotationDialogInputs` reads from `proteomics.groups` (v1.3 999.1 hotfix) | Dialog pre-fills from existing auto-grouped state; stale column refs drop silently | ✓ Good -- closes silent data-loss risk |
| Subcellular-location classifier port verbatim from CK-omics (v1.3 Phase 13) | Preserves client's exact 11-category taxonomy + reviewed-priority + GO-CC fallback | ✓ Good -- semantic parity with CK-omics deliverables |
| Gene-label resolver: ENSRNOG/ENSMUSG/ENSG/ENSDARG/MGP/LOC/RGD/AABR detection + Ensembl resolution (v1.3 Phase 14) | At-parse-time enrichment of predicted-only annotations; `*`/`†` provenance marks | ✓ Good -- maps un-readable IDs to human-readable labels |
| Filters viewer config via `columnNames` array, not typed filters spec (v1.3 Phase 13 round-3) | Typed filter spec was stripped by platform serializer; columnNames is the stable contract | ✓ Good -- viewer docks with intended columns |
| Capture-restore search wiring: `df.filter` restore + `df.selection` write (v1.3 Phase 14) | Highlight-not-hide for search matches keeps NS cloud visible | ✓ Good -- preserves volcano semantics under search |
| Two-PR / hotfix close pattern post-phase-shipping (v1.3 999.x) | Single-commit hotfixes (not full phase ceremony) for narrow UAT/review follow-ups | ✓ Good -- low overhead, atomic |

## Evolution

This document evolves at phase transitions and milestone boundaries.

**After each phase transition** (via `/gsd-transition`):
1. Requirements invalidated? -- Move to Out of Scope with reason
2. Requirements validated? -- Move to Validated with phase reference
3. New requirements emerged? -- Add to Active
4. Decisions to log? -- Add to Key Decisions
5. "What This Is" still accurate? -- Update if drifted

**After each milestone** (via `/gsd:complete-milestone`):
1. Full review of all sections
2. Core Value check -- still the right priority?
3. Audit Out of Scope -- reasons still valid?
4. Update Context with current state

---
*Last updated: 2026-06-06 -- v1.4 Cross-Team Review milestone started*
