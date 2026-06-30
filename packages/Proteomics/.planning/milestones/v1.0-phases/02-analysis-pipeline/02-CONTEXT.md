# Phase 2: Analysis Pipeline - Context

**Gathered:** 2026-02-28
**Status:** Ready for planning

<domain>
## Phase Boundary

Scientists can take imported proteomics data through the complete analysis pipeline: assign experimental groups, normalize intensity columns, impute missing values, and run differential expression. Produces log2FC and p-value columns ready for visualization in Phase 3.

</domain>

<decisions>
## Implementation Decisions

### Group Annotation UX
- Two side-by-side multi-select column pickers in a dialog — user selects intensity columns for Group A and Group B
- Two groups only (control vs treatment) — no multi-group support in this phase
- Auto-populate lists with columns that have the Proteomics-Intensity semantic type (from MaxQuant parser)
- User-editable group names via text inputs, defaulting to "Control" and "Treatment"
- Group assignments persisted as DataFrame tags (JSON) via df.setTag()

### DE Method
- Client-side Welch's t-test using @datagrok-libraries/statistics tTest(), not R limma
- Benjamini-Hochberg FDR correction only (no Bonferroni option)
- Dialog with pre-filled defaults: p-value threshold 0.05, |log2FC| threshold 1.0
- Add a boolean "significant" column marking proteins passing both thresholds
- R limma deferred to Phase 4 alongside DEqMS

### Pipeline Workflow
- Separate menu items under Proteomics: Annotate Experiment, Normalize, Impute, Differential Expression
- Each step validates prerequisites: DE checks groups annotated, impute checks normalized, etc.
- Normalization dialog: method selector (only median centering for v1)
- Imputation dialog: downshift (1.8) and width (0.3) parameter inputs with tooltips
- DE dialog: threshold inputs for significance cutoffs
- Info balloon (grok.shell.info) after each operation with summary (e.g., "Normalized 12 columns", "Imputed 847 missing values")

### Column Handling
- Normalization: in-place modification of log2 intensity columns (matches Perseus convention)
- Imputation: in-place fill of NaN/null values in intensity columns
- DE results: add new columns (log2FC, p-value, adj. p-value, significant) to existing DataFrame
- DE result columns get correct semantic types (Proteomics-Log2FC, Proteomics-PValue)

### Operation State Tracking
- Track applied operations via DataFrame tags: proteomics.normalized, proteomics.imputed, proteomics.de_complete
- Use tags for prerequisite validation (prevent double-normalization, ensure correct order)
- Group assignments stored as JSON in tag: proteomics.groups

### Claude's Discretion
- Exact dialog layout and styling
- Error message wording for prerequisite failures
- Whether to warn or block on double-operations
- Tooltip text for imputation parameters

</decisions>

<specifics>
## Specific Ideas

- DE result columns should be named descriptively: "log2FC (Control vs Treatment)", "p-value", "adj. p-value", "significant"
- Imputation parameters 1.8 and 0.3 are Perseus defaults — keep them as defaults but let users adjust
- The prerequisite check pattern should be reusable across steps

</specifics>

<code_context>
## Existing Code Insights

### Reusable Assets
- `@datagrok-libraries/statistics`: tTest() for Welch's t-test, fdrcorrection() for BH FDR correction
- `proteomics-types.ts`: SEMTYPE constants for semantic type assignment (PROTEIN_ID, GENE_SYMBOL, LOG2FC, P_VALUE, INTENSITY)
- `column-detection.ts`: findColumn() and findProteomicsColumns() for robust column lookup with fallbacks
- `maxquant-parser.ts`: Reference for column detection patterns, semantic type assignment, log2 transformation

### Established Patterns
- Column detection: semantic type first, name heuristics fallback (from column-detection.ts)
- Null values: use DG.FLOAT_NULL (not NaN) for missing floats (from maxquant-parser.ts)
- File import: DG.Utils.openFile() + grok.shell.addTableView() (from package.ts)
- Test framework: category/test/expect from @datagrok-libraries/test (from tests/parsers.ts)
- DataFrame tags: df.setTag()/getTag() for persistent metadata

### Integration Points
- package.ts: Skeleton menu handlers already exist for Annotate, Normalize, Impute, DE — currently show info stubs
- analysis/ directory: Skeleton files exist for normalization.ts, imputation.ts, differential-expression.ts
- Menu items wired via //top-menu: Proteomics | ... comment annotations

</code_context>

<deferred>
## Deferred Ideas

- R limma integration — Phase 4 alongside DEqMS
- Multi-group comparisons (ANOVA-style) — future capability
- Quantile normalization and VSN — v2 requirements (ANLY-06, ANLY-07)
- kNN imputation — v2 requirement (ANLY-08)

</deferred>

---

*Phase: 02-analysis-pipeline*
*Context gathered: 2026-02-28*
