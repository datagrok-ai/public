# Phase 4: Annotation and Alternatives - Context

**Gathered:** 2026-03-01
**Status:** Ready for planning

<domain>
## Phase Boundary

UniProt protein info panel triggered by ProteinId semantic type + DEqMS as alternative differential expression method alongside existing limma. No new analysis types, no new visualizations, no new import formats.

</domain>

<decisions>
## Implementation Decisions

### UniProt Panel Content
- Claude's discretion on detail level (essential summary vs rich detail)
- Claude's discretion on GO term display (grouped by category vs flat list)
- Include a clickable link to the full UniProt entry page (Claude decides prominent vs subtle placement)
- Claude's discretion on error state handling (inline message, retry, etc.)

### DEqMS Integration
- Claude decides how DEqMS is offered alongside limma (method selector in existing dialog vs separate menu item)
- Claude decides on peptide count column default and selection UX
- Claude decides on fallback chain (keep t-test fallback or not)
- Same output column format as limma (log2FC, p-value, adj.p-value, significant) -- method stored in DataFrame metadata tag, so existing viewers work identically

### Panel Trigger Behavior
- Claude decides how to handle multi-accession cells (e.g., 'P12345;Q67890')
- Claude decides on caching strategy for UniProt API responses
- Claude decides on ID format parsing robustness (sp|ACC|NAME, tr|ACC|NAME, bare accession)
- Claude decides on loading indicator approach (spinner vs skeleton)

### Demo & Discoverability
- Claude decides whether to extend the demo to showcase UniProt panel
- Claude decides on DEqMS menu placement (inside existing dialog vs separate entry)
- Claude decides on help text / tooltips for method selection guidance

### Claude's Discretion
User deferred all implementation decisions to Claude's judgment across all four discussion areas. Key constraint: DEqMS output columns must match limma's column names so existing volcano, heatmap, and PCA viewers work without modification.

</decisions>

<specifics>
## Specific Ideas

No specific requirements -- open to standard approaches. User trusts Claude to make scientifically and architecturally appropriate choices across all areas.

</specifics>

<code_context>
## Existing Code Insights

### Reusable Assets
- `SEMTYPE.PROTEIN_ID` ('Proteomics-ProteinId'): Already assigned to protein ID columns during MaxQuant import
- `SEMTYPE.GENE_SYMBOL` ('Proteomics-GeneSymbol'): Already assigned to gene name columns
- `findColumn()` in `utils/column-detection.ts`: Finds columns by semantic type with name-hint fallback
- `buildExpressionDf()` in `differential-expression.ts`: Builds clean expression matrix for R scripts
- `grok.decorators.panel`: Decorator available for auto-triggering panels on semantic types
- Existing limma R script (`scripts/limma_de.R`): DEqMS extends this with `spectraCounteBayes(fit)`

### Established Patterns
- Top-menu registration via `@grok.decorators.func({'top-menu': 'Proteomics | ...'})`
- DataFrame tags for state tracking: `df.setTag('proteomics.de_complete', 'true')`
- Dialog pattern: `ui.dialog().add().onOK().show()` with input widgets
- R script invocation: `grok.functions.call('Proteomics:functionName', {...})`
- Limma fallback pattern: try server-side R, catch and fall back to client-side

### Integration Points
- `package.ts`: New panel function and possible menu item registration
- `package.g.ts`: Must be manually updated for panels with `//meta.role: panel` and semType annotation
- `scripts/`: New or extended R script for DEqMS
- `differential-expression.ts`: DEqMS logic or method selection

</code_context>

<deferred>
## Deferred Ideas

None -- discussion stayed within phase scope.

</deferred>

---

*Phase: 04-annotation-and-alternatives*
*Context gathered: 2026-03-01*
