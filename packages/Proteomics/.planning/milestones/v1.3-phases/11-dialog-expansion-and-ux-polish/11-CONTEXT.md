# Phase 11: Dialog Expansion and UX Polish - Context

**Gathered:** 2026-03-07
**Status:** Ready for planning

<domain>
## Phase Boundary

Expand normalization, imputation, and DE dialogs with method selectors, conditional parameters, visual feedback, and descriptive viewer titles. Retain imported filenames as DataFrame names. No new analysis algorithms -- Phase 10 provides those; this phase wires them into polished UI.

</domain>

<decisions>
## Implementation Decisions

### Normalization dialog
- Box plot showing one box per sample for intensity distribution visualization
- Reactive preview: box plot updates when user changes normalization method (clone DataFrame, run selected method, update plot)
- Before-state only -- no side-by-side before/after; the plot shows current or previewed state
- Inline yellow/orange banner warning if Spectronaut pre-normalized data detected (`proteomics.preNormalized` tag)
- Warning text: "This data may be pre-normalized (Spectronaut). Additional normalization may distort results."
- Warning is non-blocking -- user can still proceed

### Imputation dialog
- Method selector dropdown with conditional inline parameters below:
  - MinProb: downshift + width inputs
  - kNN: k neighbors input
  - Zero/Mean/Median: no additional parameters
- Per-group minimum valid values filter: protein must have at least N valid values in at least one annotated group
- Proteins failing filter are removed (rows deleted), not hidden
- Live count updates as filter threshold changes: "Will keep X/Y proteins (Z removed)"
- Same reactive pattern as enrichment dialog's live threshold count

### DE dialog
- Single dropdown for comparison direction with auto-generated pairs: "Treatment vs Control", "Control vs Treatment"
- Dynamic hint text below dropdown: "Positive log2FC = higher in [first group], Negative log2FC = higher in [second group]"
- Comparison picker inserted at top of dialog; existing method selector, thresholds, and peptide column remain below
- All methods use consistent conditional visibility: DEqMS shows peptide column, Limma and t-test show no extra params

### Viewer titles
- Format: "Viewer: Context" -- e.g., "Volcano: Treatment vs Control", "PCA: All Groups", "Heatmap: Top 50 DE Proteins"
- Set via `setOptions({title: '...'})` or viewer creation options (proven in QC dashboard)

### DataFrame naming
- All import paths use original filename minus extension: 'proteinGroups.txt' -> 'proteinGroups', 'HYE_mix.tsv' -> 'HYE_mix'
- Consistent across MaxQuant, Spectronaut, and generic import
- Auxiliary DataFrames include source context: 'PCA: proteinGroups', 'Enrichment Results: proteinGroups', 'Heatmap: proteinGroups'

### Claude's Discretion
- Exact box plot sizing and positioning within dialog
- Clone strategy for reactive normalization preview (full clone vs column subset)
- Default minimum valid values threshold (likely 2 or 3)
- Exact wording of info messages and tooltips

</decisions>

<specifics>
## Specific Ideas

- Normalization box plot should feel like QC dashboard's intensity distribution viewer -- same Datagrok box plot viewer, just embedded in a dialog
- Imputation live count follows the enrichment dialog pattern (onChanged.subscribe -> update info div)
- DE comparison dropdown generates all permutations from `getGroups()` -- pairs are "GroupA vs GroupB" strings parsed back to group references on OK
- Keep dialog changes minimal -- insert new controls above/around existing ones rather than restructuring

</specifics>

<code_context>
## Existing Code Insights

### Reusable Assets
- `differential-expression.ts:275-280`: Conditional visibility pattern (`peptideRow.style.display` toggle on `methodInput.onChanged`)
- `enrichment.ts:344-355`: Live count update pattern (`fcInput.onChanged.subscribe(updateCount)`)
- `generic-parser.ts:86-115`: Viewer-in-dialog pattern (preview grid embedded via `.root`)
- `qc-dashboard.ts:65-109`: Viewer title pattern (`title: 'MA Plot'` in viewer options)
- `experiment-setup.ts:30-62`: Group retrieval via `getGroups()` for generating comparison pairs

### Established Patterns
- Dialogs use `ui.dialog().add(...).onOK(() => {...}).show()`
- Reactive updates via `.onChanged.subscribe()` chains
- DataFrame tags track state: `proteomics.normalized`, `proteomics.imputed`, `proteomics.preNormalized`
- Progress indicators via `DG.TaskBarProgressIndicator.create()`
- Input creation: `ui.input.choice()`, `ui.input.float()`, `ui.input.int()`, `ui.input.columns()`

### Integration Points
- Normalization dialog: `analysis/normalization.ts:29-53` -- expand existing `showNormalizationDialog()`
- Imputation dialog: `analysis/imputation.ts:52-82` -- expand existing `showImputationDialog()`
- DE dialog: `analysis/differential-expression.ts:235-360` -- expand existing `showDEDialog()`
- Viewer creation: `viewers/volcano.ts`, `viewers/heatmap.ts`, `viewers/pca-plot.ts` -- add title options
- Import naming: `parsers/maxquant-parser.ts:134`, `parsers/generic-parser.ts:185` -- standardize df.name

</code_context>

<deferred>
## Deferred Ideas

None -- discussion stayed within phase scope

</deferred>

---

*Phase: 11-dialog-expansion-and-ux-polish*
*Context gathered: 2026-03-07*
