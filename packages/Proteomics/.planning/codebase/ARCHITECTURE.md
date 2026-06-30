# Architecture

**Analysis Date:** 2026-05-11

## System Overview

```text
┌─────────────────────────────────────────────────────────────────────┐
│                  Datagrok Platform Top Menu                         │
│   "Proteomics | Import | ..." "| Analyze | ..." "| Visualize | ..." │
│   Registered via @grok.decorators.func({'top-menu': '...'})         │
│   Entry: src/package.ts (PackageFunctions class)                    │
└──────┬──────────────────┬────────────────────┬──────────────────────┘
       │                  │                    │
       ▼                  ▼                    ▼
┌────────────┐  ┌──────────────────┐  ┌───────────────────┐
│  parsers/  │  │    analysis/     │  │     viewers/      │
│ maxquant   │  │ experiment-setup │  │ volcano.ts        │
│ spectronaut│  │ normalization    │  │ heatmap.ts        │
│ generic    │  │ imputation       │  │ pca-plot.ts       │
│ shared-utils│ │ differential-    │  │ qc-dashboard.ts   │
└────┬───────┘  │  expression      │  │ enrichment-viewers│
     │          │ enrichment       │  └────────┬──────────┘
     │          │ pca.ts           │           │
     ▼          └────────┬─────────┘           │
┌──────────────────────────────────────────────┘
│         DG.DataFrame (in-memory, semantic-tagged)
│   Tags: proteomics.source, proteomics.groups (JSON),
│         proteomics.de_complete, proteomics.normalized,
│         proteomics.imputed, proteomics.enrichment,
│         proteomics.preNormalized, proteomics.de_method
└──────────────────────────────────────────────────────
         │               │
         ▼               ▼
┌──────────────┐  ┌────────────────────────────┐
│  panels/     │  │  scripts/ (R)              │
│ uniprot-panel│  │ limma_de.R                 │
│ (DG.Widget)  │  │ deqms_de.R                 │
└──────────────┘  │ vsn_normalize.R            │
                  │ Invoked via                │
                  │ grok.functions.call(...)   │
                  └────────────────────────────┘
```

## Component Responsibilities

| Component | Responsibility | Key Files |
|-----------|----------------|-----------|
| Entry point | Platform function registration; menu wiring; precondition gating | `src/package.ts` |
| Auto-generated wrappers | JSDoc-style metadata wrappers for platform function registry | `src/package.g.ts`, `src/package-api.ts` |
| Parsers | Raw file text → filtered, semantic-tagged `DG.DataFrame`; log2 transform | `src/parsers/maxquant-parser.ts`, `src/parsers/spectronaut-parser.ts`, `src/parsers/generic-parser.ts` |
| Shared parser utilities | log2 transform, log2 detection, delimiter detection, auto-suggest | `src/parsers/shared-utils.ts` |
| Experiment setup | Group annotation dialog; `proteomics.groups` tag (JSON-encoded) | `src/analysis/experiment-setup.ts` |
| Normalization | Median, quantile (client-side); VSN (R via `grok.functions.call`) | `src/analysis/normalization.ts` |
| Imputation | MinProb, kNN, Zero, Mean, Median — all client-side | `src/analysis/imputation.ts` |
| Differential expression | Welch's t-test (client-side); limma / DEqMS (R); cascading fallback chain | `src/analysis/differential-expression.ts` |
| PCA math | Jacobi eigendecomposition; sample-level DataFrame creation | `src/analysis/pca.ts` |
| Enrichment | g:Profiler gConvert + gGOSt API calls; result DataFrame construction | `src/analysis/enrichment.ts` |
| QC computations | MA, loess trend, CV, missingness matrix, unpivot helpers | `src/viewers/qc-computations.ts` |
| Viewers | Volcano, heatmap, PCA plot, QC dashboard, enrichment visualizations | `src/viewers/` |
| UniProt panel | Context panel widget fetching UniProt REST API | `src/panels/uniprot-panel.ts` |
| Column detection | Semantic-type-first column lookup with name-hint fallback | `src/utils/column-detection.ts` |
| Type constants | SEMTYPE constants; all semantic type strings live here | `src/utils/proteomics-types.ts` |
| Detectors | Platform-loaded `semTypeDetector` functions for auto-tagging columns | `detectors.js` |

## Pattern Overview

**Overall:** Functional pipeline — stateless pure functions operating on a shared `DG.DataFrame`, with workflow state carried entirely as DataFrame tags.

**Key Characteristics:**
- No class instances beyond the platform-required `PackageFunctions` and `_package` singleton
- State accumulates as `df.setTag(key, value)` calls — no external state store
- Viewers are created on demand, not maintained as long-lived objects
- Cross-DataFrame linking achieved via platform selection events (`onCurrentRowChanged`) and RxJS subscriptions
- R scripts are optional; every server-side call has a client-side fallback

## Layers

**Parsers Layer:**
- Purpose: Convert raw file text into a typed, filtered `DG.DataFrame` ready for analysis
- Location: `src/parsers/`
- Contains: `parseMaxQuantText`, `parseSpectronautText`, `showGenericImportDialog`, and shared utilities
- Depends on: `DG.DataFrame.fromCsv`, `src/utils/proteomics-types.ts`, `src/parsers/shared-utils.ts`
- Used by: `src/package.ts` menu handlers

**Analysis Layer:**
- Purpose: Mutate the DataFrame in-place (add columns, set tags) to advance the workflow state
- Location: `src/analysis/`
- Contains: `experiment-setup.ts`, `normalization.ts`, `imputation.ts`, `differential-expression.ts`, `pca.ts`, `enrichment.ts`
- Depends on: `@datagrok-libraries/statistics` (t-test, FDR); `grok.functions.call` for R scripts
- Used by: `src/package.ts` menu handlers and viewer constructors

**Viewers Layer:**
- Purpose: Create `DG.Viewer` instances that read from an already-analysed DataFrame
- Location: `src/viewers/`
- Contains: `volcano.ts`, `heatmap.ts`, `pca-plot.ts`, `qc-dashboard.ts`, `qc-computations.ts`, `enrichment-viewers.ts`
- Depends on: Analysis layer (reads tags and columns); Dendrogram package (optional, for hierarchical clustering)
- Used by: `src/package.ts` menu handlers

**Panels Layer:**
- Purpose: Context panel widget providing per-row protein metadata from external API
- Location: `src/panels/`
- Contains: `uniprot-panel.ts`
- Depends on: UniProt REST API (`fetch` directly, not via `grok.dapi.fetchProxy`)
- Used by: Platform via `@grok.decorators.panel` registration

**Utils Layer:**
- Purpose: Shared constants and cross-cutting helpers
- Location: `src/utils/`
- Contains: `proteomics-types.ts` (SEMTYPE constants), `column-detection.ts` (findColumn, findProteomicsColumns)
- Depends on: `DG` types only
- Used by: All other layers

## Data Flow

### Primary Analysis Path

1. User opens file via top-menu handler (`src/package.ts`)
2. Parser reads raw text → `parseMaxQuantText(text)` or `parseSpectronautText(text)` → returns filtered, log2-transformed `DG.DataFrame` (`src/parsers/`)
3. Parser sets `proteomics.source` tag and optionally `proteomics.preNormalized`, `proteomics.groups`
4. `grok.shell.addTableView(df)` — DataFrame enters platform session
5. User opens Annotate dialog → `showAnnotationDialog(df)` → sets `proteomics.groups` tag as JSON (`src/analysis/experiment-setup.ts`)
6. User opens Normalize → `showNormalizationDialog(df)` → mutates log2 columns in-place, sets `proteomics.normalized = 'true'` (`src/analysis/normalization.ts`)
7. User opens Impute → `showImputationDialog(df)` → fills nulls in log2 columns, sets `proteomics.imputed = 'true'` (`src/analysis/imputation.ts`)
8. User runs Differential Expression → `showDEDialog(df)` → adds `log2FC`, `p-value`, `adj.p-value`, `significant` columns, sets `proteomics.de_complete = 'true'` (`src/analysis/differential-expression.ts`)
9. Volcano and Heatmap check `proteomics.de_complete` tag before rendering (`src/viewers/volcano.ts`, `src/viewers/heatmap.ts`)
10. Enrichment checks `proteomics.de_complete`, calls g:Profiler API, produces separate enrichment `DG.DataFrame` tagged `proteomics.enrichment = 'true'` (`src/analysis/enrichment.ts`)

### R Script Invocation

1. `runLimmaDE()` builds a clean expression `DG.DataFrame` with simple column names (`s1`, `s2`, ...) to avoid encoding issues
2. Calls `grok.functions.call('Proteomics:LimmaDE', {exprDf, nGroup1, fcThreshold, pThreshold})`
3. R script (`scripts/limma_de.R`, registered as `#name: limmaDE`) runs in server R environment with `bioconductor-limma`
4. Result `DG.DataFrame` (columns: `log2FC`, `p.value`, `adj.p.value`, `significant`) is copied back to the original `df`
5. On any failure, falls back to client-side Welch's t-test (`src/analysis/differential-expression.ts:runDifferentialExpression`)

### PCA Path (Sample-Level DataFrame)

1. `computePCA(df, intensityCols, groupAssignment)` builds a transposed matrix (samples × proteins), centers, computes covariance, Jacobi-decomposes (`src/analysis/pca.ts`)
2. Returns a **new** `pcaDf` with rows = samples (not proteins) containing PC1, PC2, Group, SampleName columns
3. `grok.shell.addTableView(pcaDf)` creates a **separate** table view — PCA never shares a table view with the protein DataFrame
4. `createPcaPlot(df, ...)` returns `{viewer, pcaDf}` pair; caller docks viewer into the pcaDf table view (`src/viewers/pca-plot.ts`)

### Cross-DataFrame Selection (Enrichment → Volcano)

1. `wireEnrichmentToVolcano(enrichDf, proteinDf)` subscribes to `enrichDf.onCurrentRowChanged` (`src/viewers/enrichment-viewers.ts`)
2. When a term row is selected, parses `Intersection` column → set of gene symbols
3. Looks up matching rows in `proteinDf` by gene symbol column
4. Sets `proteinDf.selection` bits → volcano/heatmap highlights linked proteins
5. Subscriptions stored in module-level `activeSubscriptions[]` and unsubscribed on re-open

**State Management:**
- Workflow state carried in `DG.DataFrame` tags (string key-value pairs); no separate state object
- Group assignments serialized as JSON in `proteomics.groups` tag
- Tags accumulate monotonically; no step can be undone without recreating the DataFrame

## Key Abstractions

**GroupAssignment Interface:**
- Purpose: Describes which intensity columns belong to each of two experimental groups
- Definition: `{group1: {name: string; columns: string[]}; group2: {name: string; columns: string[]}}`
- Persisted as: `df.getTag('proteomics.groups')` (JSON string)
- Used by: normalization, imputation, DE, PCA, heatmap, QC dashboard, enrichment
- File: `src/analysis/experiment-setup.ts`

**SEMTYPE Constants:**
- Purpose: Type-safe semantic type strings shared across parsers, analysis, and detectors
- Values: `PROTEIN_ID = 'Proteomics-ProteinId'`, `GENE_SYMBOL = 'Proteomics-GeneSymbol'`, `LOG2FC = 'Proteomics-Log2FC'`, `P_VALUE = 'Proteomics-PValue'`, `INTENSITY = 'Proteomics-Intensity'`
- File: `src/utils/proteomics-types.ts`

**findColumn / findProteomicsColumns:**
- Purpose: Column lookup strategy: semantic type first, name-hint heuristic second
- Pattern: `findColumn(df, SEMTYPE.GENE_SYMBOL, ['gene name', 'gene symbol'])`
- File: `src/utils/column-detection.ts`

## Entry Points

**PackageFunctions class:**
- Location: `src/package.ts`
- Triggers: Platform top-menu clicks; `@grok.decorators.func`, `@grok.decorators.init`, `@grok.decorators.panel`
- Responsibilities: Input validation (`grok.shell.tv?.dataFrame`), precondition tag checks, delegation to domain modules, viewer docking

**detectors.js:**
- Location: `detectors.js` (root, plain JS — loaded separately by platform before the bundle)
- Triggers: Platform column semantic type detection on DataFrame load
- Responsibilities: Auto-assign `Proteomics-ProteinId`, `Proteomics-GeneSymbol`, `Proteomics-Log2FC`, `Proteomics-PValue`, `Proteomics-Intensity` semantic types by column name heuristic

**R Scripts:**
- Location: `scripts/limma_de.R`, `scripts/deqms_de.R`, `scripts/vsn_normalize.R`
- Triggers: `grok.functions.call('Proteomics:LimmaDE', ...)`, `grok.functions.call('Proteomics:DeqmsDE', ...)`, `grok.functions.call('Proteomics:VsnNormalize', ...)`
- Environment: Conda environment with `bioconductor-limma` (declared via `#environment:` header in R scripts)

## Architectural Constraints

- **Threading:** Single-threaded browser event loop; all client-side computation (kNN imputation, Jacobi PCA, quantile normalization) blocks the UI thread. `DG.TaskBarProgressIndicator` provides visual feedback but does not off-load work.
- **Global state:** Module-level `activeSubscriptions: rxjs.Subscription[]` in `src/viewers/enrichment-viewers.ts` — must be unsubscribed before re-opening enrichment visualization to avoid duplicate handlers.
- **DataFrame mutability:** All analysis steps mutate the original `df` in-place. This means steps cannot be undone without re-importing. The heatmap intentionally clones a subset DataFrame (`df.clone(filter)`) to isolate its sort/filter state from the protein-level view.
- **PCA creates a new DataFrame:** PCA produces a sample-level DataFrame with a different row count than the protein DataFrame. This must be opened in a separate `TableView` — never added as a viewer to the protein table view.
- **Circular imports:** `src/analysis/normalization.ts` imports `unpivotIntensities` from `src/viewers/qc-computations.ts` (analysis → viewers direction); this is the one cross-layer dependency that runs against the natural data flow direction.
- **Dendrogram dependency:** Heatmap hierarchical clustering uses `DG.Func.find({package: 'Dendrogram', name: 'hierarchicalClustering'})` — a runtime lookup for an optional cross-package dependency. Falls back to significance-based sort if Dendrogram is unavailable.

## Anti-Patterns

### Using `fetch()` Directly for External APIs

**What happens:** `src/panels/uniprot-panel.ts` and `src/analysis/enrichment.ts` call `fetch()` directly for UniProt REST and g:Profiler APIs.
**Why it's wrong:** The Datagrok platform provides `grok.dapi.fetchProxy()` for external HTTP calls to avoid CORS issues in the browser environment; raw `fetch()` may work in development but fail in production deployments.
**Do this instead:** Replace `fetch(url, ...)` with `await grok.dapi.fetchProxy(url, ...)` in `src/panels/uniprot-panel.ts:fetchUniProtData` and `src/analysis/enrichment.ts:fetchWithTimeout`.

### Mutation of the In-Place DataFrame Without Idempotency Guard

**What happens:** `showDEDialog` guards against re-running (`proteomics.de_complete` check), but normalization and imputation add columns without removing previous ones. Re-calling these on an already-processed DataFrame would duplicate columns.
**Why it's wrong:** Duplicate `log2FC`, `p-value`, `adj.p-value` columns would corrupt downstream viewers silently.
**Do this instead:** Use the `ensureFreshFloat` pattern from `src/viewers/qc-computations.ts` — remove the column if it exists before adding a new one.

## Error Handling

**Strategy:** Cascading try/catch with graceful degradation to client-side alternatives.

**Patterns:**
- R script calls (`runLimmaDE`, `runDeqmsDE`, `vsnNormalize`) all have `try/catch` blocks that fall back to client-side equivalents with `grok.shell.warning(...)` notification
- DEqMS → limma → t-test cascade in `showDEDialog` (`src/analysis/differential-expression.ts:314-391`)
- VSN → quantile fallback in `src/analysis/normalization.ts:vsnNormalize`
- Heatmap prerequisite failures and Dendrogram unavailability logged via `console.warn()`

## Cross-Cutting Concerns

**Logging:** `console.warn()` for recoverable errors (R unavailable, Dendrogram unavailable, PCA ellipse failure); `grok.shell.info()` for success messages; `grok.shell.warning()` for precondition violations; `grok.shell.error()` for hard failures.
**Validation:** Precondition checks at the entry point layer (`src/package.ts`) before delegating to domain modules; redundant checks in some domain functions (e.g., `createExpressionHeatmap` checks `proteomics.de_complete` again).
**Progress indicators:** `DG.TaskBarProgressIndicator.create(...)` / `.update(...)` / `.close()` in `finally` blocks for long-running operations (heatmap, kNN imputation, DE).

---

*Architecture analysis: 2026-05-11*
