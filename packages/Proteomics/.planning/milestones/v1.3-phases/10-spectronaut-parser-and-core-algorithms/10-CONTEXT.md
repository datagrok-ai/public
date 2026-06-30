# Phase 10: Spectronaut Parser and Core Algorithms - Context

**Gathered:** 2026-03-07
**Status:** Ready for planning

<domain>
## Phase Boundary

Spectronaut long-format TSV import with pivot-to-wide, plus quantile normalization, VSN normalization, kNN imputation, and zero/mean/median imputation as callable functions registered with Datagrok. Existing dialogs remain unchanged -- Phase 11 expands them with method selectors and advanced UX.

</domain>

<decisions>
## Implementation Decisions

### Spectronaut pivot strategy
- Use PG.IBAQ as the quantity column for protein-level intensities (no aggregation needed)
- Name pivoted sample columns as `R.Condition + "_" + R.Replicate` (e.g., "HYE mix A_1")
- Store R.FileName as a column tag on each pivoted intensity column for traceability
- Q-value filtering: configurable threshold parameter with default 0.01; non-numeric values ('Profiled', 'NaN') treated as passing
- Contaminant/decoy filtering: filter by PG.ProteinGroups prefix (CON__, REV__), matching MaxQuant parser pattern

### Pre-normalized data handling
- Detect pre-normalization using existing `detectLog2Status()` heuristic from shared-utils.ts
- If data detected as pre-normalized, skip log2 transform and use `copyAsLog2Columns()` instead
- Tag as `proteomics.preNormalized = 'true'` (new tag, separate from `proteomics.normalized`)
- Retain both raw and log2 intensity columns in output DataFrame (VSN needs raw intensities)

### kNN imputation behavior
- Configurable k parameter with default k=10
- Row-wise (protein-wise) imputation: find k most similar proteins based on intensity pattern across samples
- Euclidean distance on shared non-missing values
- When fewer than k neighbors available, use available neighbors; fall back to column mean if zero neighbors
- Progress indicator via `grok.shell.startProgress()` with percentage updates per protein row

### Functions and Datagrok registration
- Export pure TypeScript functions: `quantileNormalize()`, `vsnNormalize()`, `imputeKnn()`, `imputeZero()`, `imputeMean()`, `imputeMedian()`
- Register all functions via `//name:` metadata so they appear in Datagrok's function browser
- Existing normalization/imputation dialogs remain unchanged -- Phase 11 adds method selectors
- VSN auto-falls back to quantile normalization if R is unavailable (matches limma/DEqMS fallback pattern)

### Claude's Discretion
- Exact distance weighting scheme for kNN (uniform vs inverse-distance)
- VSN R script parameter choices (lts.quantile, etc.)
- Spectronaut protein ID column extraction and primary column logic

</decisions>

<specifics>
## Specific Ideas

- Follow MaxQuant parser pattern exactly: filter -> clone -> assign semtypes -> transform intensities
- Reuse shared-utils.ts functions: `log2TransformColumns()`, `copyAsLog2Columns()`, `detectLog2Status()`, `addPrimaryColumnIfNeeded()`
- VSN R script should follow limma_de.R pattern: `grok.functions.call()` with DataFrame arguments, try/catch with fallback

</specifics>

<code_context>
## Existing Code Insights

### Reusable Assets
- `shared-utils.ts`: `log2TransformColumns()`, `copyAsLog2Columns()`, `detectLog2Status()`, `detectDelimiter()`, `addPrimaryColumnIfNeeded()`, `autoSuggestProteinIdColumn()`
- `maxquant-parser.ts`: Complete parser pattern (filter markers, clone, semtypes, transform) to mirror
- `normalization.ts`: `medianNormalize()` -- pattern for new normalization functions (in-place, `getRawData()`, `fireValuesChanged()`, tag setting)
- `imputation.ts`: `imputeMinProb()` -- pattern for new imputation functions (per-column iteration, tag setting, return count)
- `proteomics-types.ts`: `SEMTYPE` constants for semantic type assignment
- `experiment-setup.ts`: `setGroups()`/`getGroups()` for group annotation persistence via DataFrame tags

### Established Patterns
- DataFrame tags for state: `proteomics.normalized`, `proteomics.imputed`, `proteomics.groups`
- Intensity columns use `SEMTYPE.INTENSITY` and `log2()` prefix naming convention
- R scripts in `scripts/` directory with `//language: r` metadata, called via `grok.functions.call()`
- `DG.DataFrame.fromCsv()` with `columnImportOptions` to force column types

### Integration Points
- Parser registered in `package.ts` as file handler (like MaxQuant)
- New functions registered via `//name:` metadata in package source
- Group auto-population feeds into `experiment-setup.ts` via `setGroups()`

</code_context>

<deferred>
## Deferred Ideas

None -- discussion stayed within phase scope

</deferred>

---

*Phase: 10-spectronaut-parser-and-core-algorithms*
*Context gathered: 2026-03-07*
