# Phase 8: Gene ID Mapping & Enrichment Analysis - Context

**Gathered:** 2026-03-06
**Status:** Ready for planning

<domain>
## Phase Boundary

Scientists can map protein IDs to gene symbols and run GO/KEGG/Reactome overrepresentation analysis on their DE results, with results displayed in an interactive table. This phase covers ID mapping, enrichment computation, and result table display. Enrichment visualization (dot plots, bar charts) and volcano plot integration (highlight proteins by term) are separate phases.

</domain>

<decisions>
## Implementation Decisions

### Enrichment Workflow
- Menu item: `Proteomics | Enrichment Analysis...` opens a config dialog
- Dialog fields: significance thresholds (FC, p-value), enrichment sources (GO BP/MF/CC, KEGG, Reactome), organism dropdown
- Auto-fill thresholds from existing DE columns (log2FC, adj.p-value) with defaults FC>=1, padj<=0.05
- All three enrichment sources (GO, KEGG, Reactome) enabled by default — user can uncheck any
- Live count of significant proteins updates as thresholds change: "42 of 1,200 proteins are significant with these thresholds"
- Progress dialog blocks while API call runs: "Running enrichment analysis..."

### Result Presentation
- Results open as a new Datagrok table view (follows PCA pattern for different-shaped data)
- Full schema: Source (GO:BP/GO:MF/GO:CC/KEGG/REAC), Term ID, Term Name, P-value, FDR, Gene Count, Gene Ratio, Intersection (member genes)
- Show all results, no pre-filtering — user can sort/filter manually
- Significant terms (FDR < 0.05) highlighted but not hidden

### ID Mapping
- Gene mapping is automatic and transparent — triggered as part of enrichment, not a separate menu item
- If a gene name column with SEMTYPE.GENE_SYMBOL already exists, use it directly — skip API mapping
- If no gene names exist, map UniProt accessions to gene symbols via g:Profiler g:Convert
- Add a 'Gene Symbol (mapped)' column to the main table after mapping (only when no gene name column exists)
- Show mapping statistics in progress dialog: "850/1,200 proteins mapped (70.8%). 350 unmapped."

### Organism Selection
- Organism dropdown in the enrichment dialog, defaulting to Homo sapiens
- Short curated list of 8-10 common proteomics organisms: human, mouse, rat, yeast, E. coli, zebrafish, Drosophila, Arabidopsis, C. elegans

### Claude's Discretion
- Exact highlighting mechanism for significant terms (color coding, row tags, etc.)
- Error handling for API failures (timeout, rate limit, no results)
- Exact wording of mapping statistics message
- How to handle proteins that map to multiple gene symbols

</decisions>

<specifics>
## Specific Ideas

- g:Profiler REST API handles both ID mapping (g:Convert) and enrichment (g:GOSt) — pure TypeScript, no R scripts needed
- Use `parseAccession()` from uniprot-panel.ts to extract clean accessions before sending to g:Profiler
- Background set for enrichment must be all quantified proteins (not whole genome) — use g:Profiler's `domain_scope: "custom"` parameter

</specifics>

<code_context>
## Existing Code Insights

### Reusable Assets
- `parseAccession()` in panels/uniprot-panel.ts: extracts clean UniProt accessions from sp|ACC|NAME, semicolons, CPTAC formats — use to prepare protein list for g:Profiler
- `SEMTYPE` constants in utils/proteomics-types.ts: GENE_SYMBOL, PROTEIN_ID already defined
- `findColumn()` in utils/column-detection.ts: finds columns by semantic type — use to detect existing gene name columns

### Established Patterns
- Top menu entries via `@grok.decorators.func({'top-menu': 'Proteomics | ...'})`
- Dialog construction: `ui.dialog()`, `ui.input.choice()`, `ui.input.bool()`, `ui.input.float()`
- Separate table view for different-shaped results (PCA pattern)
- Progress indicator pattern from DE analysis
- REST API fetch pattern from UniProt panel (async fetch with error handling)

### Integration Points
- New menu entry in `PackageFunctions` class in package.ts
- Reads DE result columns (log2FC, adj.p-value, significant) from current DataFrame
- Reads protein ID column (SEMTYPE.PROTEIN_ID) for accession extraction
- New enrichment DataFrame added as separate table view via `grok.shell.addTableView()`
- Optional new gene symbol column added to source DataFrame

</code_context>

<deferred>
## Deferred Ideas

- Enrichment visualization (dot plots, bar charts) — Phase 9 (VIZ-02, VIZ-03)
- Volcano plot integration (select GO term -> highlight proteins) — Phase 9 (ENRICH-04)
- GSEA (rank-based enrichment) — tracked as ENRICH-05 in future requirements
- Functional annotation clustering — tracked as ENRICH-06 in future requirements
- Standalone ID mapping menu item — add if users request it

</deferred>

---

*Phase: 08-gene-id-mapping-enrichment-analysis*
*Context gathered: 2026-03-06*
