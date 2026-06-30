# Phase 1: Data Import and Foundation - Context

**Gathered:** 2026-02-28
**Status:** Ready for planning

<domain>
## Phase Boundary

Import MaxQuant proteinGroups.txt into Datagrok with automatic contaminant/reverse filtering, intensity column auto-detection (LFQ/iBAQ/Reporter), log2 transformation, semantic type assignment, and a bundled demo dataset. This phase delivers the data foundation that all downstream analysis and visualization depends on.

</domain>

<decisions>
## Implementation Decisions

### Import mechanism
- Use menu item pattern (`Proteomics | Import | MaxQuant...`) already scaffolded in package.ts
- Follow Datagrok file-handler pattern for auto-detection of proteinGroups.txt format
- Parser function signature already defined: `parseMaxQuant(file: DG.FileInfo): Promise<DG.DataFrame>`

### Filtering
- Filter contaminants (CON__ prefix or `+` in "Potential contaminant" column), reverse hits (REV__ prefix or `+` in "Reverse" column), and only-identified-by-site rows
- Handle both MaxQuant 1.x and 2.x column naming conventions

### Intensity columns
- Auto-detect intensity type by prefix pattern: `LFQ intensity`, `Intensity`, `iBAQ`, `Reporter intensity`
- Keep raw intensity columns and add log2-transformed columns alongside them
- Assign Proteomics-Intensity semantic type to intensity columns

### Semantic types
- Use existing SEMTYPE constants from `src/utils/proteomics-types.ts`
- Existing detectors in `detectors.js` handle detection for subsequently opened tables
- Parser should explicitly assign semantic types on import (don't rely solely on detectors)

### Demo dataset
- Use a small, publicly available MaxQuant proteinGroups.txt (e.g., PXD000561 HeLa benchmark or similar)
- Bundle in `files/demo/` directory
- Keep small enough to load quickly (<5MB)

### Claude's Discretion
- Exact import dialog UI layout and options
- Column naming convention for log2-transformed columns
- How to display filtering summary (info balloon vs silent)
- Whether to auto-open a table view after import
- Demo dataset selection (any public MaxQuant output that's small and representative)
- Whether to show before/after protein count after filtering
- Error handling for malformed files

</decisions>

<specifics>
## Specific Ideas

- Product spec specifies parsing semicolon-delimited Protein IDs into primary ID (first entry)
- Product spec specifies parsing semicolon-delimited Gene Names into primary gene name
- Scaffold already has the full file structure: parsers/, utils/, detectors.js
- Research identified that the DEResult interface should be deleted in favor of DataFrames — but that's Phase 2 scope

</specifics>

<code_context>
## Existing Code Insights

### Reusable Assets
- `src/utils/proteomics-types.ts`: SEMTYPE constants already defined (PROTEIN_ID, GENE_SYMBOL, LOG2FC, P_VALUE, INTENSITY)
- `src/utils/column-detection.ts`: `findColumn()` and `findProteomicsColumns()` with semtype + name hint fallback
- `detectors.js`: 5 semantic type detectors already implemented (detectProteinId, detectGeneSymbol, detectLog2FC, detectPValue, detectIntensity)
- `src/package.ts`: Menu item stubs already wired for `importMaxQuant()` and `proteomicsDemo()`

### Established Patterns
- Parser pattern: async function taking `DG.FileInfo`, returning `Promise<DG.DataFrame>`
- Function metadata: JSDoc-style `//name:`, `//tags:`, `//top-menu:` comments processed by `grok api`
- Decorator polyfill in package.ts (lines 13-36) — work within this pattern, don't change it
- Files accessed via `_package.files.readAsText('demo/...')` at runtime

### Integration Points
- `src/package.ts` `importMaxQuant()` method — wire to parser
- `src/package.ts` `proteomicsDemo()` method — wire to demo dataset loading
- `src/package.ts` `initProteomics()` — register semantic types if needed
- Menu structure: `Proteomics | Import | MaxQuant...`

</code_context>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 01-data-import-and-foundation*
*Context gathered: 2026-02-28*
