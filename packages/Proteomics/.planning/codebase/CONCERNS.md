# Codebase Concerns

**Analysis Date:** 2026-05-11

---

## Performance

### R Session Cold-Start on Every DE Invocation

- **Issue:** The Datagrok scripting model spawns a fresh R process per `grok.functions.call()`. Loading base R + the limma package (BiocGenerics, Biobase, limma — ~3–5 MB of R bytecode) takes several seconds on every invocation.
- **Files:** `scripts/limma_de.R:18`, `scripts/deqms_de.R:24`, `src/analysis/differential-expression.ts:135`
- **Impact:** Every DE run with the "limma" or "DEqMS" method pays a 3–10 second cold-start penalty before any computation begins. This is the dominant latency in the DE workflow.
- **Fix approach:** Add an `environments/` directory with a Conda/Bioconductor environment spec declaring `bioconductor-limma` (and `bioconductor-deqms` for that path) so the platform can pre-warm a cached R session. Also consider a lightweight R health-check call during `initProteomics()` at `src/package.ts:56` to warm the R worker on package init. (Partial: `scripts/limma_de.R` already declares `#environment:` as of commit `489bb53ec9`.)

### DEqMS Sequential Double-R-Call Fallback

- **Issue:** When "DEqMS" is selected and `runDeqmsDE` fails (DEqMS not installed), the catch block at `src/analysis/differential-expression.ts:351` immediately calls `runLimmaDE`, incurring a second full R cold-start before reaching the client-side fallback.
- **Files:** `src/analysis/differential-expression.ts:340-370`
- **Impact:** Up to 2× the R cold-start penalty (potentially 6–20 seconds) for a method failure the user never knew happened.
- **Fix approach:** Probe R/limma availability once before the dialog is shown (or on first call), cache the result, and skip the sequential fallback chain when R is known unavailable.

### kNN Imputation O(N²) Per Imputed Row

- **Issue:** `imputeKnn` computes Euclidean distances from each row-with-missing to every other row, iterating over all `nRows * nRows` pairs in a nested loop in pure JS.
- **Files:** `src/analysis/imputation.ts:131-152`
- **Impact:** For a dataset with 3000 proteins (common in DDA proteomics), each imputed row requires 3000 distance computations. With many missing rows this becomes O(N²) with no batching or early exit. The progress indicator (`pi.update`) is called per row but does not yield, so long runs block the event loop.
- **Fix approach:** Batch the progress indicator yields with `await new Promise(resolve => setTimeout(resolve, 0))` every 50 rows to keep the UI responsive. Consider using Web Workers or switching to `@datagrok-libraries/ml` nearest-neighbor utilities if available.

### Dendrogram Package Cold-Start in Heatmap

- **Issue:** `createExpressionHeatmap` calls `DG.Func.find({package: 'Dendrogram', name: 'hierarchicalClustering'})` at `src/viewers/heatmap.ts:138-145`. On the first call in a session, the Dendrogram package may not be loaded, requiring the platform to initialize it — adding 2–10 seconds of latency.
- **Files:** `src/viewers/heatmap.ts:136-150`
- **Impact:** Heatmap clustering stalls silently on first invocation per session. The progress indicator created inside `createExpressionHeatmap` at line 136 is in addition to the one in `package.ts:152`, so there are two nested progress indicators — minor redundancy.
- **Fix approach:** No timeout is set on the `dFunc.apply()` call; if the Dendrogram function is found but the package is unresponsive, the call hangs without bound. Add a `Promise.race` with a timeout and fall back to significance-based sort as the existing catch block does.

### Sequential Worker Dispatches in Distance Matrix

- **Issue:** The now-commented-out Dendrogram tree-helper path iterated `calcDistanceMatrix` once per sample column sequentially (12+ Worker round-trips for 12 samples). The current `dFunc.apply` approach at `src/viewers/heatmap.ts:138-145` replaces this with a single function call; however the Worker dispatch overhead for the Dendrogram package's internal column loop remains.
- **Files:** `src/viewers/heatmap.ts:138-168` (large commented-out block on lines 151-169 — dead code that should be removed)
- **Impact:** Dead code at lines 151-169 clutters the file and documents an abandoned approach, creating confusion about what the actual algorithm is.
- **Fix approach:** Delete the commented-out block at lines 151-169.

### VSN Normalization R Cold-Start

- **Issue:** `vsnNormalize` calls `grok.functions.call('Proteomics:VsnNormalize', ...)` at `src/analysis/normalization.ts:134`. `scripts/vsn_normalize.R` has a hard `library(vsn)` with no fallback probing.
- **Files:** `src/analysis/normalization.ts:103-153`, `scripts/vsn_normalize.R:10`
- **Impact:** VSN normalization pays full R cold-start for every invocation. The JS catch block correctly falls back to quantile normalization, but no user-visible progress is shown during the R attempt — the dialog appears to hang for several seconds before the fallback warning appears.
- **Fix approach:** Show a progress indicator before the `grok.functions.call` in `vsnNormalize`. The `showNormalizationDialog` onOK handler at line 229 calls `await vsnNormalize(df, selected)` with no surrounding progress indicator.

---

## UX Debt

### QC Dashboard — Progress Indication Not Phase-Aware

- **Issue:** `openQcDashboard` creates a `DG.TaskBarProgressIndicator` at `src/viewers/qc-dashboard.ts:22` but never calls `pi.update(...)` between the compute steps (MA, loess trend, CV, viewer creation). The thin task-bar indicator fires once and stays static.
- **Files:** `src/viewers/qc-dashboard.ts:21-145`
- **Impact:** From user perspective, clicking "QC Dashboard" produces several seconds of UI silence. Identified as a demo-blocking issue at Cytokinetics 2026-05-04. Filed as todo `2026-05-04-add-progress-indication-to-qc-dashboard.md`.
- **Fix approach:** Call `pi.update(pct, message)` between each named compute step (MA, loess, CV, viewer creation). Target: user sees progress move within 1 second of menu click.

### QC Dashboard — Layout Navigability

- **Issue:** The docked grid layout in `openQcDashboard` (`src/viewers/qc-dashboard.ts:133-142`) tiles 7 viewers in a single flat level using `DG.DOCK_TYPE.RIGHT` and `DOWN`. For a typical 1920×1080 screen this produces very narrow panels and makes the heatmap-like missing-value grid unreadable.
- **Files:** `src/viewers/qc-dashboard.ts:133-142`
- **Impact:** Filed as todo `2026-05-04-reorganize-qc-dashboard-layout-for-navigability.md`. Biologics reviewers demoed to on 2026-05-04 could not identify the missing-values heatmap as such.
- **Fix approach:** Redesign the dock layout into two rows with proportional sizing. Consider collapsible/tabbed arrangement for secondary plots.

### Heatmap — Commented-Out Code from Prior Debug

- **Issue:** `src/viewers/heatmap.ts:151-169` contains a 19-line block of commented-out code from a refactored Dendrogram integration approach. The block imports `getTreeHelper`, `getDendrogramService`, and `DistanceMetric` which are no longer used.
- **Files:** `src/viewers/heatmap.ts:1-7` (unused imports), `src/viewers/heatmap.ts:151-169` (dead block)
- **Impact:** Dead imports increase bundle size slightly and confuse future maintainers. The `getDendrogramService` import from `@datagrok-libraries/bio/src/trees/dendrogram` at line 6 and `getTreeHelper`/`DistanceMetric` at lines 6-7 are imported but never called in the current code path.
- **Fix approach:** Remove lines 5-7 imports and delete the commented block at lines 151-169.

### Analytics-Level Help Missing

- **Issue:** No dialog in the package (normalization, imputation, DE, enrichment) has embedded help text, links, or tooltips explaining why each parameter matters statistically. Only method-level tooltips exist (e.g. `methodInput.setTooltip(...)` in `src/analysis/differential-expression.ts:278`).
- **Files:** `src/analysis/differential-expression.ts`, `src/analysis/normalization.ts`, `src/analysis/imputation.ts`, `src/analysis/enrichment.ts`
- **Impact:** Filed as todo `2026-05-04-add-explanatory-help-to-all-analytics-dialogs-and-viewers.md`. First-time users cannot determine what "downshift 1.8, width 0.3" means in the MinProb imputation dialog without external references.
- **Fix approach:** Add a `ui.link("?", "https://...")` help link or inline `ui.divText(...)` explanation beneath each parameter section in dialogs.

### DE Dialog — No Re-Run Capability

- **Issue:** `showDEDialog` at `src/analysis/differential-expression.ts:242` returns early if `proteomics.de_complete === 'true'` with a warning. There is no way to re-run DE with different parameters without manually clearing the tag.
- **Files:** `src/analysis/differential-expression.ts:242-245`, `src/package.ts:116-126`
- **Impact:** Filed as todo `2026-03-05-reset-de-analysis-to-allow-re-running-with-different-parameters.md`. Users who choose wrong thresholds on first run are blocked.
- **Fix approach:** Replace the early-return with a confirmation dialog ("DE already performed. Re-run with new parameters? This will overwrite existing results."). Remove existing DE-result columns before re-running.

---

## Data Model Gaps

### Two-Group-Only Experiment Model

- **Issue:** `GroupAssignment` in `src/analysis/experiment-setup.ts:8-11` is hard-coded to exactly two groups (`group1`, `group2`). The annotation dialog, DE, QC, heatmap, and PCA all assume this exact structure.
- **Files:** `src/analysis/experiment-setup.ts:8-11`, `src/analysis/differential-expression.ts:237-244`, `src/viewers/heatmap.ts:38-40`, `src/viewers/qc-computations.ts:35-38`
- **Impact:** Cannot handle time-course experiments, dose-response series, multi-condition screens, or any study design with more than two groups. Spectronaut auto-group only fires when exactly 2 conditions are detected (`src/parsers/spectronaut-parser.ts:129`).
- **Fix approach:** Extend `GroupAssignment` to support N groups, or introduce a separate `ComparisonPair` abstraction for DE that references two groups from an N-group design.

### No Multi-Experiment / Longitudinal Data Model

- **Issue:** The package has no concept of "screening run," "experiment batch," "week," or "campaign." Every analysis is a standalone single-DataFrame operation with no linkage to other analyses.
- **Files:** `src/analysis/experiment-setup.ts`, `src/package.ts`
- **Impact:** Blocks the three pending use cases filed 2026-05-11: cross-week compound comparison (`2026-05-11-compare-most-impactful-compounds-across-weekly-screening-runs.md`), SPC/longitudinal QC (`2026-05-11-track-control-sample-performance-across-weekly-assays-with-spc.md`), and cross-team sharing filed by target (`2026-05-11-share-analysis-read-only-with-biologics-team-filed-by-target.md`). All three require a data model that links multiple separate analyses.
- **Fix approach:** Introduce a "screening campaign" or "study" entity as a Datagrok project tag or a separate metadata table that groups individual analysis DataFrames under a common identifier.

### No Compound Entity

- **Issue:** There is no representation of a chemical compound in the package. Sample columns are named by condition+replicate strings only; no SMILES, InChI, or compound ID is attached.
- **Files:** `src/analysis/experiment-setup.ts`, `src/parsers/maxquant-parser.ts`, `src/parsers/spectronaut-parser.ts`
- **Impact:** Blocks the compound-comparison use case and any integration with the Datagrok Chem package (structure rendering, SAR, MCS). Filed as part of `2026-05-11-compare-most-impactful-compounds-across-weekly-screening-runs.md`.
- **Fix approach:** Add an optional "compound" annotation step in `showAnnotationDialog`, accepting SMILES or compound ID, stored as DataFrame/column tags.

---

## Feature Gaps

### Label-Free DDA/DIA Only — No TMT/iTRAQ/PTM

- **Issue:** Both parsers (`maxquant-parser.ts`, `spectronaut-parser.ts`) focus on LFQ/iBAQ intensities for label-free data. There is no support for tandem mass tag (TMT) or iTRAQ quantification (reporter ion intensities), phosphoproteomics (PTM site localization), or cross-linking MS data.
- **Files:** `src/parsers/maxquant-parser.ts:7` (only `lfq intensity`, `ibaq`, `reporter intensity`, `intensity` prefixes), `src/parsers/spectronaut-parser.ts:10` (PG.IBAQ only)
- **Impact:** The package cannot be used for isobaric-labeling workflows common in pharma proteomics. `reporter intensity` is listed in `INTENSITY_PREFIXES` at `src/parsers/maxquant-parser.ts:7` but there is no label deconvolution logic or TMT channel mapping.
- **Fix approach:** Add a TMT/iTRAQ parser path that reads reporter intensity columns, maps channel → sample label, and normalizes across channels.

### DE Method Coverage — No ANOVA, LME, or Paired Tests

- **Issue:** DE supports: limma moderated t-test, DEqMS, client-side Welch's t-test. No ANOVA (for N>2 groups), linear mixed effects, or paired-sample tests are available.
- **Files:** `src/analysis/differential-expression.ts:273-275` (method choice list)
- **Impact:** Time-course experiments, repeated-measures designs, and experiments with batch covariates cannot be analyzed correctly.
- **Fix approach:** Add at least one server-side paired/batch-corrected option (e.g. limma with blocking factor) and document the limitation in dialog tooltips.

### Enrichment — No GSEA or Gene Set Scoring

- **Issue:** Enrichment is over-representation analysis (ORA) via g:GOSt only. No gene set enrichment analysis (GSEA), CAMERA, or fgsea for ranked-list approaches. No support for custom gene sets.
- **Files:** `src/analysis/enrichment.ts:92-131`
- **Impact:** ORA is threshold-dependent and discards rank information. Filed implicitly in the analytics-help todo; researchers often expect GSEA.
- **Fix approach:** Add a "GSEA" method option that passes the full ranked protein list to a server-side R script using fgsea.

### Volcano Plot — No Auto-Zoom

- **Issue:** `createVolcanoPlot` does not adjust axis extents after adding threshold lines. If a dataset has extreme outliers, the threshold lines at `|log2FC| = 1.0` and `-log10(0.05) ≈ 1.3` may appear in the bottom-left corner of a large scatter.
- **Files:** `src/viewers/volcano.ts:71-121`
- **Impact:** Filed as todo `2026-03-05-auto-zoom-volcano-plot-to-show-all-points-and-threshold-lines.md`. First impression of the volcano is often a blank-looking plot until the user zooms manually.
- **Fix approach:** Set `sp.props.xMin`, `xMax`, `yMin`, `yMax` based on data range with 10% padding, ensuring threshold lines are always visible.

---

## External Coupling

### Raw `fetch()` to UniProt Without Proxy

- **Issue:** `fetchUniProtData` in `src/panels/uniprot-panel.ts:74-83` calls `fetch('https://rest.uniprot.org/...')` directly from the browser, bypassing the Datagrok proxy.
- **Files:** `src/panels/uniprot-panel.ts:75-83`
- **Impact:** Will fail with CORS errors in Datagrok environments that block direct browser-to-external-API calls. The CLAUDE.md rule explicitly states: use `grok.dapi.fetchProxy()` for external URLs to avoid CORS. Fails silently — the panel shows "Unable to fetch UniProt data for X" with no indication of why.
- **Fix approach:** Replace `fetch(url)` with `grok.dapi.fetchProxy(url)` in `fetchUniProtData`.

### Raw `fetch()` to g:Profiler Without Proxy

- **Issue:** `fetchWithTimeout` in `src/analysis/enrichment.ts:54-71` calls `fetch(url, {...options, signal: ...})` directly for both the `gConvert` and `gGOSt` API calls.
- **Files:** `src/analysis/enrichment.ts:54-71`, `src/analysis/enrichment.ts:77`, `src/analysis/enrichment.ts:99`
- **Impact:** Same CORS risk as the UniProt panel. g:Profiler calls will silently fail in proxy-only environments. The error message "g:Profiler API returned status X" or "request timed out" gives the user no path to resolution.
- **Fix approach:** Replace `fetch(url, options)` with `grok.dapi.fetchProxy(url, options)` in `fetchWithTimeout`. Note: `fetchProxy` may not support `AbortController`/`signal` — the timeout logic may need restructuring.

### UniProt Panel — No Result Caching

- **Issue:** `uniprotPanel` at `src/panels/uniprot-panel.ts:187-198` fetches UniProt data on every selection change with no caching. Clicking through a list of proteins fires a new HTTP request per protein.
- **Files:** `src/panels/uniprot-panel.ts:73-84`, `src/panels/uniprot-panel.ts:187-198`
- **Impact:** High network traffic and latency when browsing DE results. Filed as todo `2026-03-03-cache-uniprot-protein-content-to-avoid-repeated-api-calls.md`.
- **Fix approach:** Add a `Map<string, UniProtEntry>` module-level cache in `uniprot-panel.ts`. Cache keyed by accession; evict after N entries or a TTL.

### Dendrogram Cross-Package Dependency

- **Issue:** `src/viewers/heatmap.ts` depends on the Dendrogram package being installed and responsive via `DG.Func.find({package: 'Dendrogram', ...})`. The package has no version constraint or install check.
- **Files:** `src/viewers/heatmap.ts:138-145`
- **Impact:** If Dendrogram is not installed, the catch block fires silently and falls back to significance-based sort — correct behavior, but the user sees no indication that clustering was skipped. If Dendrogram is installed but a different version, the `inputs.length !== 4` guard at line 139 is the only version check.
- **Fix approach:** Surface a user-visible note (e.g. `grok.shell.info(...)`) when the Dendrogram fallback fires, so the user understands why the heatmap is not clustered.

---

## Code Quality

### Auto-Generated Files Must Not Be Edited

- **Files:** `src/package.g.ts` (110 lines), `src/package-api.ts` (107 lines)
- **Issue:** These files are auto-generated by `grok api`. Both files have a header comment indicating they are generated. Any manual edits to these files will be overwritten on the next `grok api` run.
- **Fix approach:** Do not edit these files directly. Place any extensions in `src/package.ts` using the `@grok.decorators.*` pattern.

### Decorator Polyfill in package.ts

- **Issue:** `src/package.ts:26-49` contains a hand-rolled polyfill that injects missing `grok.decorators.*` methods at runtime. The polyfill is guarded by `if (!grok.decorators)` but the inner `forEach` also checks each individual decorator — meaning any platform that has a partial `grok.decorators` object gets silently patched, potentially masking version-incompatibility bugs.
- **Files:** `src/package.ts:24-51`
- **Impact:** Runtime behavior depends on which Datagrok platform version is installed. If the platform adds a decorator with different semantics, the polyfill will shadow it.
- **Fix approach:** Pin the minimum required `datagrok-api` version in `package.json` that natively supports all required decorators, and remove the polyfill.

### QC Dashboard Progress Indicator Closed in `finally` Around Synchronous Code

- **Issue:** `openQcDashboard` at `src/viewers/qc-dashboard.ts:21-145` wraps all computation in `try/finally { pi.close(); }`. All the computation (MA, loess, CV, viewer creation, docking) is synchronous. The `pi.close()` fires immediately after the synchronous block completes — before the browser has a chance to paint the progress bar. The user likely never sees the indicator.
- **Files:** `src/viewers/qc-dashboard.ts:22`, `src/viewers/qc-dashboard.ts:143-145`
- **Impact:** The progress indicator is effectively invisible for the current synchronous implementation. Confirmed in the todo filed 2026-05-04.
- **Fix approach:** Wrap compute steps in `await new Promise(resolve => setTimeout(resolve, 0))` yields between stages so the event loop can paint the indicator before each heavy step.

### `showNormalizationDialog` Preview Leaks DG.Viewer Instances

- **Issue:** `updatePreview()` in `src/analysis/normalization.ts:197-218` creates a new `DG.Viewer.boxPlot(...)` each time the method or column selection changes, without calling `dispose()` on the previous viewer.
- **Files:** `src/analysis/normalization.ts:197-218`
- **Impact:** Each change to the method or column picker leaks a DG.Viewer. On a dialog that stays open while the user experiments, this can accumulate many undisposed viewers.
- **Fix approach:** Track the current preview viewer in a `let currentPreview: DG.Viewer | null = null` variable; call `currentPreview?.detach()` or `.dispose()` before creating a new one.

---

## Test Coverage Gaps

### R-Dependent Paths Not Tested

- **What's not tested:** `runLimmaDE`, `runDeqmsDE`, `vsnNormalize` (R call paths). All three fall back gracefully in test environments (R unavailable), but the actual R output parsing (column renaming from `p.value` → `p-value`, `adj.p.value` → `adj.p-value` at `src/analysis/differential-expression.ts:148-151`) is not covered by any test.
- **Files:** `src/analysis/differential-expression.ts:126-172`, `src/analysis/normalization.ts:103-153`, `src/tests/analysis.ts`
- **Risk:** Column name mismatches between the R output and the JS copy loop would silently produce all-null DE result columns, breaking downstream visualization with no error.
- **Priority:** High — the R output schema is fragile and not version-pinned.

### Heatmap and Viewer Rendering Not Tested

- **What's not tested:** `createExpressionHeatmap`, `createVolcanoPlot`, `createPcaPlot`, `openQcDashboard`. No test file covers viewer creation or Dendrogram integration.
- **Files:** `src/viewers/heatmap.ts`, `src/viewers/volcano.ts`, `src/viewers/pca-plot.ts`, `src/viewers/qc-dashboard.ts`
- **Risk:** Grid configuration, heatmap column visibility logic, Dendrogram fallback path, and formula-line addition are all untested. Regressions in viewer setup would only surface in manual testing.
- **Priority:** Medium — add at minimum a test that `createExpressionHeatmap` returns a `DG.Grid` with the correct column visibility for a known dataset.

### UniProt Panel and Enrichment API Calls Not Tested

- **What's not tested:** `fetchUniProtData`, `gConvert`, `gGOSt`, `runEnrichmentPipeline`. These are network-dependent and not mocked in any test.
- **Files:** `src/panels/uniprot-panel.ts:73-84`, `src/analysis/enrichment.ts:73-132`, `src/tests/enrichment.ts`
- **Risk:** Response-shape changes in the UniProt or g:Profiler APIs break panels/enrichment silently. The `buildEnrichmentDf` and `countSignificantProteins` functions are tested in isolation, but the full pipeline including network calls is not.
- **Priority:** Medium — add tests with mocked `fetch` responses (or `grok.dapi.fetchProxy`) that verify the JSON parsing and DataFrame construction paths.

### PCA Computation Not Tested

- **What's not tested:** `computePCA` in `src/analysis/pca.ts`. No test file for PCA computation.
- **Files:** `src/analysis/pca.ts`
- **Risk:** The custom Jacobi eigendecomposition (`pca.ts:33-119`) is a non-trivial algorithm implemented from scratch. Off-by-one errors in the covariance matrix transpose or eigenvalue sorting would silently produce wrong PC scores with no warning.
- **Priority:** High — the Jacobi implementation has no external validation. Add a test that verifies PC1 explains >50% of variance for a dataset with a known dominant axis, and that PC scores match a reference implementation (even a simple 2×2 case).

### DE Dialog Fallback Chain Not Integration-Tested

- **What's not tested:** The DEqMS → limma → client-side t-test fallback chain in `showDEDialog` (`src/analysis/differential-expression.ts:340-390`). Unit tests cover `runDifferentialExpression` (client-side) only.
- **Files:** `src/analysis/differential-expression.ts:340-390`
- **Risk:** The fallback chain's intermediate `pi.update()` calls, tag setting (`proteomics.de_method`), and `onComplete` callback invocation are all untested. A regression in error propagation could silently run the wrong method without telling the user.
- **Priority:** Medium.

---

## Scaling Limits

### UniProt Panel Fires on Every Row Change

- **Issue:** The `uniprotPanel` is a Datagrok context panel that fires on every `currentRow` change. No debounce or rate limiting is applied.
- **Files:** `src/panels/uniprot-panel.ts:187-198`
- **Impact:** Scrolling through a 1000-protein grid at keyboard-repeat speed fires 1000 fetch requests in rapid succession. Without caching or debounce, this saturates the network and the UniProt rate-limit.
- **Scaling path:** Add a 300ms debounce before triggering the fetch, combined with the caching fix above.

### kNN Imputation Blocking on Large Datasets

- **Issue:** `imputeKnn` in `src/analysis/imputation.ts:58-179` runs entirely on the JS main thread with no async yields. For 3000 proteins × 24 samples with 30% missing, this is approximately 900 rows × 2999 distance computations = ~2.7M inner-loop iterations in a blocking loop.
- **Files:** `src/analysis/imputation.ts:112-174`
- **Impact:** Large datasets (common in DDA proteomics with deep fractionation) will freeze the browser tab for several seconds. No hard limit is enforced.
- **Scaling path:** Move kNN to a Web Worker, or suggest a faster imputation method (MinProb, Mean) for datasets > 500 proteins × 12 samples, with a warning in the dialog.

---

*Concerns audit: 2026-05-11*
