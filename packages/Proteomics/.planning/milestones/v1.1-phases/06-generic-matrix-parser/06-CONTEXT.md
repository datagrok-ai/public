# Phase 6: Generic Matrix Parser - Context

**Gathered:** 2026-03-05
**Status:** Ready for planning

<domain>
## Phase Boundary

Scientists can import any CSV/TSV proteomics matrix into the full analysis pipeline via a column mapping dialog. After import, the DataFrame has the same semantic types as MaxQuant import, enabling the full downstream pipeline (annotation, normalization, imputation, DE) without modification.

</domain>

<decisions>
## Implementation Decisions

### File Detection & Entry Point
- Top menu entry only: "Proteomics | Import | Generic Matrix..." — no file handler sniffing
- Accept CSV and TSV files (.csv, .tsv, .txt) with auto-detected delimiter
- Use Datagrok's built-in dataset input widget (table selector + file picker + folder browser + database icons) — NOT a bare OS file dialog
- Backport the same dataset input widget to the MaxQuant importer for consistency

### Column Mapping Dialog
- Single dialog with all fields visible at once
- Layout: dataset input (top), protein ID column dropdown, gene name column dropdown (optional), intensity columns multi-select, log2 transform toggle, data preview grid (bottom)
- Auto-suggest intensity columns by keyword matching: intensity, lfq, ibaq, tmt, reporter, abundance
- Auto-suggest protein ID column by keyword matching: protein, accession, uniprot
- Live-updating preview: 5 rows, updates as column selections change to show only mapped columns
- Import and Cancel buttons at bottom

### Log2 Transform
- Auto-detect whether data is already log2-transformed by checking value ranges (0-30 = log2, 1e3+ = raw)
- Pre-set the log2 toggle based on detection, with a hint message (e.g., "Data appears to be raw intensities")
- User can override the auto-detection
- Keep both original intensity columns AND log2-prefixed copies (mirror MaxQuant behavior)
- Always rename with log2() prefix for downstream compatibility, even if data was already log2-transformed

### Downstream Compatibility
- Assign SEMTYPE.INTENSITY to both original and log2 columns
- Assign SEMTYPE.PROTEIN_ID to selected protein ID column
- Assign SEMTYPE.GENE_SYMBOL to selected gene name column (if provided)
- Create 'Primary Protein ID' / 'Primary Gene Name' columns only when semicolons detected in source data
- No row filtering — generic import keeps all rows as-is
- Tag source format: `proteomics.source = 'generic'` — backport `proteomics.source = 'maxquant'` to MaxQuant importer
- DataFrame name reflects source dataset/file name — backport this to MaxQuant importer too

### Claude's Discretion
- Exact auto-detection threshold for log2 vs raw intensities
- Column picker filter configuration details
- Preview grid styling and column ordering
- Error handling for malformed files

</decisions>

<specifics>
## Specific Ideas

- Use Datagrok's built-in dataset input widget (the one with table selector dropdown, file picker, folder browser, and database icons) — user provided a screenshot reference
- Unify both import paths (MaxQuant and Generic) to use the same dataset input widget for consistent UX

</specifics>

<code_context>
## Existing Code Insights

### Reusable Assets
- `processIntensityColumns()` in maxquant-parser.ts: log2 transform + semType assignment logic — extract into shared utility
- `addPrimaryColumn()` in maxquant-parser.ts: semicolon-delimited field extraction — reuse for generic parser
- `assignSemanticTypes()` in maxquant-parser.ts: protein/gene semType assignment — generalize for user-selected columns
- `SEMTYPE` constants in proteomics-types.ts: all semantic type strings already defined
- `findColumn()` / `findProteomicsColumns()` in column-detection.ts: column discovery by semType + name hints

### Established Patterns
- Downstream pipeline requires BOTH `semType === SEMTYPE.INTENSITY` AND `name.startsWith('log2(')` for columns to be picked up by normalization, imputation, experiment-setup
- `DG.DataFrame.fromCsv()` handles CSV/TSV parsing with delimiter and columnImportOptions
- `ui.dialog()`, `ui.input.column()`, `ui.input.columns()`, `ui.input.bool()` for dialog construction
- `df.plot.grid()` for data preview grids

### Integration Points
- New menu entry in `PackageFunctions` class in package.ts alongside existing `importMaxQuant()`
- Shared log2 transform utility imported by both maxquant-parser.ts and new generic-parser.ts
- Output DataFrame feeds into existing annotation, normalization, imputation, DE pipeline unchanged

</code_context>

<deferred>
## Deferred Ideas

- Tool-specific parsers (FragPipe, DIA-NN, Spectronaut, Proteome Discoverer) — tracked as IMPORT-07 through IMPORT-10 in future requirements
- Optional row filter column in import dialog — add if requested
- Excel (.xlsx) file support — add if requested

</deferred>

---

*Phase: 06-generic-matrix-parser*
*Context gathered: 2026-03-05*
