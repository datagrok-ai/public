# Chem Package

## Summary

The Chem package (`@datagrok/chem`, v1.17.4) is the primary cheminformatics extension for Datagrok.
It provides molecule rendering, substructure/similarity search, fingerprint computation, property
calculation, R-group analysis, scaffold trees, activity cliffs, matched molecular pair analysis, and
chemical space visualization — all running in the browser via RDKit WASM with optional Docker-based
server computation. Handles datasets of up to 10 million molecules.

User documentation: [Cheminformatics](../../help/datagrok/solutions/domains/chem/chem.md)

## Architecture

### Key files

| File                                                     | Description                                            |
|----------------------------------------------------------|--------------------------------------------------------|
| `src/package.ts`                                         | Entry point — registers all functions, viewers, panels |
| `src/chem-searches.ts`                                   | Fingerprint & similarity search orchestration          |
| `src/constants.ts`                                       | Enums, tags, configuration constants                   |
| `src/scripts-api.ts`                                     | Namespace wrapper for calling package functions        |
| `src/rendering/rdkit-cell-renderer.ts`                   | RDKit-based grid cell renderer with LRU caching        |
| `src/rendering/rdkit-reaction-renderer.ts`               | Chemical reaction rendering                            |
| `src/widgets/chem-substructure-filter.ts`                | Substructure filter with sketcher integration          |
| `src/widgets/scaffold-tree.ts`                           | Scaffold tree viewer                                   |
| `src/widgets/properties.ts`                              | Chemical properties panel                              |
| `src/widgets/drug-likeness.ts`                           | Lipinski rule compliance panel                         |
| `src/widgets/toxicity.ts`                                | Toxicity risk assessment panel                         |
| `src/widgets/identifiers.ts`                             | Identifier mapping (PubChem, ChEMBL, etc.)             |
| `src/widgets/structure2d.ts`                             | 2D structure viewer                                    |
| `src/widgets/structure3d.ts`                             | 3D structure viewer (NGL)                              |
| `src/widgets/synthon-search.ts`                          | Synthon-based searching                                |
| `src/analysis/activity-cliffs.ts`                        | Activity cliff detection & visualization               |
| `src/analysis/chem-similarity-viewer.ts`                 | Interactive similarity search viewer                   |
| `src/analysis/chem-diversity-viewer.ts`                  | Diversity analysis viewer                              |
| `src/analysis/chem-space.ts`                             | Chemical space (t-SNE/UMAP) dimensionality reduction   |
| `src/analysis/r-group-analysis.ts`                       | R-group decomposition                                  |
| `src/analysis/mpo.ts`                                    | Multi-parameter optimization                           |
| `src/analysis/molecular-matched-pairs/`                  | MMPA: fragments, generations, differences, viewer      |
| `src/analysis/bit-birch/`                                | BitBIRCH clustering (WASM-based)                       |
| `src/rdkit-service/rdkit-service.ts`                     | Web Worker manager for parallel RDKit operations       |
| `src/rdkit-service/rdkit-service-worker-substructure.ts` | Substructure search worker                             |
| `src/rdkit-service/rdkit-service-worker-similarity.ts`   | Fingerprint computation worker                         |
| `src/open-chem/ocl-sketcher.ts`                          | OpenChemLib sketcher wrapper                           |
| `src/open-chem/ocl-cell-renderer.ts`                     | OpenChemLib cell renderer                              |
| `src/open-chem/sdf-importer.ts`                          | SDF file import                                        |
| `src/open-chem/ocl-service/`                             | OCL property calculations in workers                   |
| `src/docker/api.ts`                                      | Communication with chem-chem Docker container          |
| `src/descriptors/descriptors-calculation.ts`             | Docker-based descriptor computation                    |
| `src/utils/chem-common-rdkit.ts`                         | RDKit module initialization and management             |
| `src/utils/convert-notation-utils.ts`                    | SMILES/SMARTS/MolBlock conversion                      |
| `src/utils/sdf-utils.ts`                                 | SDF file handling                                      |
| `src/utils/most-common-subs.ts`                          | Maximum Common Substructure (MCS)                      |
| `src/file-importers/smi-importer.ts`                     | SMILES/SMARTS file import                              |
| `src/file-importers/mol2-importer.ts`                    | Tripos MOL2 format import                              |
| `src/mpo/`                                               | MPO profiles, scoring, context panels                  |


### Web worker architecture

`RdKitService` manages 2–4 web workers (scaled to CPU cores). Work is distributed in striped batches
(worker 1 gets rows 0, n, 2n…; worker 2 gets rows 1, n+1, 2n+1…). Progress is tracked via events
with termination signals.

### Docker integration

The `chem-chem` Docker container provides server-side computation:
- `/chem/descriptors` — chemical descriptor calculation
- `/chem/descriptors/tree` — descriptor metadata
- `/chem/molecules_to_canonical` — SMILES canonicalization
- `/chem/molecules_to_inchi` — SMILES to InChI
- `/chem/molecules_to_inchi_key` — InChI Key generation
- `/chem/inchi_to_smiles` — InChI to SMILES

Requests are gzipped for efficiency. See `src/docker/api.ts`.

### JS API

The `dg.chem` namespace (`public/js-api/src/chem.ts`) exposes:

| Function / Class                  | Description                                            |
|-----------------------------------|--------------------------------------------------------|
| `chem.Sketcher`                   | Molecule sketcher widget (supports OCL, Ketcher, etc.) |
| `chem.getSimilarities(col, mol)`  | Similarity scores for a column against a molecule      |
| `chem.findSimilar(col, mol)`      | Top-N molecules by similarity with cutoff              |
| `chem.diversitySearch(col)`       | Select diverse molecule subset                         |
| `chem.searchSubstructure(col, q)` | Substructure search returning BitSet                   |
| `chem.rGroup(df, col, core)`      | R-group decomposition                                  |
| `chem.mcs(df, col)`               | Maximum Common Substructure                            |
| `chem.descriptors(df, col, desc)` | Chemical descriptor calculation                        |
| `chem.svgMol(smiles, w, h)`       | Render molecule to SVG                                 |
| `chem.canvasMol(...)`             | Render molecule to canvas via RDKit                    |
| `chem.drawMolecule(mol, w, h)`    | Create molecule visualization element                  |
| `chem.convert(s, from, to)`       | Convert between SMILES/SMARTS/MolBlock                 |
| `chem.checkSmiles(s)`             | Validate SMILES string                                 |
| `chem.Notation` enum              | `Smiles`, `Smarts`, `MolBlock`, `V3KMolBlock`, etc.    |

RDKit JS API interfaces are in `@datagrok-libraries/chem-meta` (`libraries/chem-meta/src/rdkit-api.ts`):
`RDModule`, `RDMol`, `MolList`, `RDReaction`, `RDSubstructLibrary`, `RGroupDecomp`.

### Public functions (package.ts)

Core rendering: `chemCellRenderer()`, `rdKitCellRenderer()`, `rdKitReactionRenderer()`,
`rdKitMixtureRenderer()`

Search: `getMorganFingerprints()`, `getMorganFingerprint()`, `getSimilarities()`,
`getDiversities()`, `findSimilar()`

Analysis: `rGroupAnalysis()`, `chemSpace()`, `scaffoldTreeViewer()`, `substructureFilter()`

Conversion: `recalculateCoords()`, `getRdKitModule()`, `getMolFileHandler()`, `canvasMol()`,
`drawMolecule()`, `getCLogP()`

### Key constants and tags

Fingerprint types: Morgan (2048 bits, radius 2), RDKit, MACCS, AtomPair, TopologicalTorsion, Pattern

Substructure filter modes: Exact Match, Contains, Included In, Is Similar, Not Contains, Not Included In

Column tags: `.%scaffold-col`, `.%chem-scaffold-align`, `.%chem-scaffold-highlight`,
`chem-apply-filter-sync`, `.%chem-space-embedding-col`

### Dependencies

| Dependency                       | Purpose                                     |
|----------------------------------|---------------------------------------------|
| `@datagrok-libraries/chem-meta`  | RDKit JS API interfaces, molfile parsing    |
| `@datagrok-libraries/ml`         | ML utilities, activity cliffs               |
| `@datagrok-libraries/math`       | High-performance math                       |
| `@datagrok-libraries/statistics` | Statistics, MPO                             |
| `@datagrok-libraries/utils`      | BitArray, caching                           |
| `openchemlib` (7.2.3)            | Alternative renderer & sketcher             |
| NGL Viewer (2.4.0)               | 3D structure visualization                  |
| RDKit WASM                       | Core cheminformatics engine (via chem-meta) |

### File format support

Input: SMILES, SMARTS, Molblock V2000/V3000, InChI, MOL2, SDF
Output: SMILES, Molblock, SDF, PNG/SVG rendering

## Usage

Samples: [Chemistry API samples](../../packages/ApiSamples/scripts/domains/chem/)

```typescript
// Substructure search
const hits = await grok.chem.searchSubstructure(col, 'c1ccccc1');

// Similarity search
const similar = await grok.chem.findSimilar(col, 'CCO', {limit: 10, cutoff: 0.5});

// R-group decomposition
const rgroups = await grok.chem.rGroup(df, 'smiles', 'c1ccccc1');

// Render molecule
const svg = grok.chem.svgMol('CCO', 200, 150);

// Convert notation
const molblock = grok.chem.convert('CCO', DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock);
```

## Test

```
grok test
```

Test files are in `src/tests/` (20+ files covering rendering, search, filters, MMPA, importers,
notation, scaffolds, and UI).
