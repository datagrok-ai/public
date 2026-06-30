# Project Retrospective

*A living document updated after each milestone. Lessons feed forward into future planning.*

## Milestone: v1.0 — MVP

**Shipped:** 2026-03-05
**Phases:** 5 | **Plans:** 12 | **Commits:** 36

### What Was Built
- MaxQuant proteinGroups.txt parser with auto-detection, filtering, log2 transform, semantic types
- Full analysis pipeline: group annotation, median normalization, MinProb imputation, Welch's t-test + BH FDR
- Interactive visualizations: volcano plot, expression heatmap with clustering, PCA colored by group
- UniProt info panel with auto-trigger on protein ID columns
- DEqMS/limma R scripts with three-level fallback chain
- Phase 5 hardening: detector regex fixes, heatmap filter isolation, progress indicators

### What Worked
- Leveraging Datagrok's built-in viewers (scatter plot for volcano, heatmap, PCA) instead of custom JsViewers saved significant effort
- Three-level DE fallback (DEqMS -> limma -> client-side t-test) gave resilience without complexity
- DataFrame clone pattern cleanly solved heatmap filter isolation from volcano
- Synthetic demo data enabled deterministic tests with controlled edge cases
- Milestone audit before Phase 5 identified real integration risks that were then fixed

### What Was Inefficient
- Dendrogram rendering attempted but not fully working — TreeHelper/DendrogramService API proved harder than expected
- Some ROADMAP.md progress tracking fell out of sync (checkboxes not updated as plans completed)
- UAT revealed gene name fallback gap that could have been caught by tests earlier
- Audit found integration risks (R script name case mismatch, detector regex) that could have been caught by integration tests

### Patterns Established
- Factory function pattern for viewer creation (createVolcanoPlot, createExpressionHeatmap, createPcaPlot)
- DataFrame metadata tags for pipeline state tracking (proteomics.normalized, proteomics.imputed, proteomics.de_complete)
- Column detection utility with semantic type -> name hint -> null fallback chain
- Function metadata comment pattern (replacing decorator polyfill from scaffold)
- Z-score normalization via temporary columns to preserve original data

### Key Lessons
1. Always verify R function call names match what grok api generates (PascalCase) — silent fallbacks hide mismatches
2. Detector regex must be tested against all column name formats upfront (space-delimited, parenthesized, prefixed)
3. When viewers share a DataFrame, filter isolation via clone is cleaner than trying to coordinate filter state
4. Milestone audits are valuable — they caught 3 high/medium integration risks before shipping

### Cost Observations
- Total execution time: ~34 minutes across 12 plans
- Average plan execution: 2.8 minutes
- Phase 5 (hardening) was fastest at 1.7 min/plan — fixing known issues is faster than building new features
- Timeline: 6 calendar days

---

## Milestone: v1.1 — Parser & QC

**Shipped:** 2026-03-07
**Phases:** 2 | **Plans:** 5 | **Commits:** 16

### What Was Built
- Generic CSV/TSV matrix import with column mapping dialog, auto-suggestion, log2 detection, and live preview
- QC dashboard with 7 linked viewers in tiled dock layout (MA plot, trend line, CV, correlation, missing heatmap, missing bar, box plot)
- Shared parser utilities extracted from MaxQuant parser for reuse
- 18 unit tests (11 generic parser + 7 QC computation)

### What Worked
- Extracting shared utilities first (Plan 06-01) before building the generic parser made 06-02 straightforward
- DG.Viewer factory methods solved the auxiliary DataFrame viewer creation problem cleanly
- Moving-average approximation for loess trend was sufficient — no need for a full statistical library
- Verification passes caught real bugs (BigInt column support, auto-suggest gaps) before they became UAT failures

### What Was Inefficient
- BigInt column type not considered initially — required a verification-phase fix across all shared utility functions
- Phase 7 needed a gap closure plan (07-03) to fix viewer creation pattern discovered during UAT
- Dock layout API (DockManager.dock()) required trial-and-error to discover that refNode must be a DockNode, not a Viewer

### Patterns Established
- File import pattern: openFile -> parse -> show mapping dialog -> assign semTypes -> addTableView
- Auxiliary DataFrame pattern: create DF locally, use DG.Viewer.factory(df, opts), dock directly — never addTable
- Column cleanup on re-run: remove and recreate computed columns to prevent accumulation
- ensureFreshFloat helper for idempotent column operations

### Key Lessons
1. Always test with BigInt column types — Datagrok BIG_INT columns return BigInt values that break Math.* functions
2. DG.Viewer factory methods are the correct way to create viewers for auxiliary DataFrames — avoid grok.shell.addTable()
3. DockManager.dock() refNode parameter requires DockNode (return value from previous dock()), not Viewer objects
4. Showing all numeric columns as candidates (not just keyword matches) gives users more flexibility

### Cost Observations
- Total execution time: ~22 minutes across 5 plans (excluding human verification pauses)
- Average plan execution: 4.4 minutes
- Phase 6 slower (16h wall clock) due to human verification pause between plans
- Timeline: 2 calendar days

---

## Milestone: v1.2 — Biological Interpretation

**Shipped:** 2026-03-07
**Phases:** 2 | **Plans:** 4 | **Commits:** 12

### What Was Built
- g:Profiler API integration for gene ID mapping (UniProt to gene symbols, 9 organisms)
- GO/KEGG/Reactome overrepresentation analysis with proteomics-specific background set
- Enrichment dot plot and bar chart with top-N filtering
- Cross-DataFrame volcano integration (select enrichment term -> highlight member proteins)
- 15 unit tests (7 enrichment + 8 enrichment visualization)

### What Worked
- Reusing ScatterPlotViewer as dot plot via sizeColumnName/colorColumnName — zero custom viewer code
- g:Profiler REST API provided everything needed (mapping, enrichment, custom background) in one service
- onCurrentRowChanged subscription for cross-DataFrame selection was clean and required no tight coupling
- Module-level subscription cleanup pattern prevented memory leaks on repeated dashboard opens

### What Was Inefficient
- Nothing significant — fastest milestone yet (7 min execution time across 4 plans)
- g:GOSt nested response structure (data.result[0].result) required a fallback handler but was easily handled

### Patterns Established
- External REST API client with AbortController timeout pattern
- Cross-DataFrame selection wiring via onCurrentRowChanged + Intersection column parsing
- Tag-based table discovery: grok.shell.tables.find with tag check for cross-table references
- Module-level subscription array with cleanup-on-reopen for memory safety

### Key Lessons
1. Reusing existing Datagrok viewers with configuration beats building custom viewers — ScatterPlot as dot plot worked perfectly
2. g:Profiler's adjusted p-values are already FDR-corrected — both P-value and FDR columns get the same value
3. Tag-based DataFrame discovery (grok.shell.tables.find) is the right pattern for cross-table references

### Cost Observations
- Total execution time: ~7 minutes across 4 plans
- Average plan execution: 1.75 minutes (fastest milestone)
- Benefit of established patterns: plans executed with zero deviations
- Timeline: 1 calendar day

---

## Milestone: v1.3 — Richer Dialogs & UX

**Shipped:** 2026-03-10
**Phases:** 2 | **Plans:** 6 | **Commits:** 20

### What Was Built
- Spectronaut long-format TSV parser with auto-detection, pivot, CON__/REV__ filtering, and auto-group annotation
- Quantile normalization (client-side), VSN normalization (R script with fallback), kNN imputation with progress indicator
- Multi-method normalization dialog with reactive box plot preview and pre-normalized data warning
- Multi-method imputation dialog with conditional parameter visibility and valid-values protein filter
- DE dialog with comparison direction picker, t-test as explicit method, and dynamic FC hint text
- Descriptive viewer titles and filename-based DataFrame naming across all import/viewer paths

### What Worked
- Two-phase structure (algorithms first → dialogs second) gave clean dependency chain with no circular imports
- Zero new npm dependencies — quantile norm and kNN in pure TypeScript, VSN via existing R infrastructure
- Established patterns from v1.0-v1.2 (viewer factory, parser utilities, column detection) made v1.3 plans execute with zero deviations
- Milestone audit (16/16 requirements, 5/5 E2E flows) provided confidence before archival
- Gap closure via 11-04 (preNormalized tag fix) caught by UAT before milestone close

### What Was Inefficient
- Phase 11 needed a 4th plan (11-04) for a 1-line tag placement fix — could have been caught by testing both log2 and raw-intensity Spectronaut paths in 10-01
- ROADMAP.md progress table got out of sync (Phase 10/11 rows had inconsistent milestone/status columns)
- 6 unused imports in package.ts accumulated as functions were added to modules but called indirectly via dialogs

### Patterns Established
- Reactive dialog preview: clone df, apply method on clone, unpivot, render box plot, destroy on input change
- Inline warning banner pattern: conditional display via getTag check with non-blocking UX
- Container-div visibility toggle: wrap related inputs in ui.div, toggle display on method change
- Comparison direction parsing: split selected string on ' vs ' to determine numerator/denominator
- Parser no-name policy: parsers return unnamed DataFrames, import handlers set df.name from filename
- Viewer title passthrough: viewer functions accept optional title parameter, callers provide context

### Key Lessons
1. Spectronaut always normalizes its data regardless of export format — unconditional tagging is correct
2. Test both log2 and raw-intensity paths through all parsers to catch conditional-branch gaps early
3. Reactive previews should clone with only needed columns (df.clone(null, selected)) for performance
4. ANY-group threshold for valid-values filter prevents over-aggressive protein removal

### Cost Observations
- Total execution time: ~16 minutes across 6 plans
- Average plan execution: 2.7 minutes
- Fastest plans: 11-01 and 11-02 at 2 min (dialog expansion on top of existing algorithms)
- Timeline: 2 calendar days

---

## Cross-Milestone Trends

### Process Evolution

| Milestone | Commits | Phases | Key Change |
|-----------|---------|--------|------------|
| v1.0 | 36 | 5 | Initial milestone — established patterns |
| v1.1 | 16 | 2 | Shared utility extraction, DG.Viewer factory pattern |
| v1.2 | 12 | 2 | External API integration, cross-DF wiring |
| v1.3 | 20 | 2 | Multi-method dialogs, Spectronaut parser, reactive previews |

### Cumulative Quality

| Milestone | Tests | Known Gaps | Tech Debt Items |
|-----------|-------|------------|-----------------|
| v1.0 | 15 (unit) + 6 (UAT) | 1 (dendrogram) | 4 (low severity) |
| v1.1 | 33 (unit) + 10 (UAT) | 1 (dendrogram, carried) | 0 |
| v1.2 | 48 (unit) + 10 (UAT) | 1 (dendrogram, carried) | 0 |
| v1.3 | 62+ (unit) + 16 (UAT) | 1 (dendrogram, carried) | 1 (unused imports, cosmetic) |

### Top Lessons (Verified Across Milestones)

1. Milestone audits catch integration risks that unit tests miss — always audit before shipping
2. Extract shared utilities early — doing it before building the consumer keeps both clean
3. DG.Viewer factory methods, not addTable+addViewer, for auxiliary DataFrames
4. Reuse existing Datagrok viewers with configuration before building custom ones — confirmed v1.0 (volcano), v1.2 (dot plot)
5. Established patterns compound: execution speed stable at 2-3 min/plan from v1.2 onward thanks to accumulated patterns
6. Test both branches of conditional logic (e.g., log2 vs raw-intensity paths) — confirmed by v1.3 gap closure (11-04)
