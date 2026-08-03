# Chem Package

## Summary

The Chem package (`@datagrok/chem`, v1.17.6) is the primary cheminformatics extension for Datagrok.
It provides molecule rendering, substructure/similarity search, fingerprint computation, property
calculation, R-group analysis, scaffold trees, activity cliffs, matched molecular pair analysis,
multi-parameter optimization (MPO), reaction enumeration and transformation, BitBIRCH/MCS
clustering, pharmacophore features, structural alerts, ChemProp ML, and chemical space
visualization — all running in the browser via RDKit WASM with optional Docker-based server
computation. Handles datasets of up to 10 million molecules.

User documentation: [Cheminformatics](../../help/datagrok/solutions/domains/chem/chem.md)

## Architecture

### Key files

| File                                                     | Description                                            |
|----------------------------------------------------------|--------------------------------------------------------|
| `src/package.ts`                                         | Entry point — registers all functions, viewers, panels via `PackageFunctions` class |
| `src/chem-searches.ts`                                   | Fingerprint & similarity search orchestration          |
| `src/constants.ts`                                       | Enums, tags, configuration constants                   |
| `src/scripts-api.ts`                                     | Namespace wrapper for calling package functions        |
| `src/demo/demo.ts`                                       | Demo entry points (overview, MMPA, similarity, scaffold, etc.) |
| `src/rendering/rdkit-cell-renderer.ts`                   | RDKit-based grid cell renderer with LRU caching        |
| `src/rendering/rdkit-reaction-renderer.ts`               | Reaction rendering — multi-step routes, branched/convergent paths via `BRANCH_DELIMITER` and `parseMultiStepReaction` |
| `src/rendering/render-molecule.ts`                       | High-level molecule rendering helper                   |
| `src/rendering/molecule-label.ts`                        | Molecule label overlays                                |
| `src/rendering/mixture-cell-renderer.ts`                 | Cell renderer for chemical mixtures                    |
| `src/rendering/mol3d-hover-handler.ts`                   | 3D molecule hover/preview handler                      |
| `src/rendering/atom-picker-controller.ts`                | Atom picker controller (interactive atom selection)    |
| `src/widgets/chem-substructure-filter.ts`                | Substructure filter with sketcher integration          |
| `src/widgets/scaffold-tree.ts`                           | Scaffold tree viewer (with node reordering)            |
| `src/widgets/scaffold-tree-filter.ts`                    | Scaffold-tree-driven filter                            |
| `src/widgets/properties.ts`                              | Chemical properties panel                              |
| `src/widgets/biochem-properties-widget.ts`               | Biochem properties panel                               |
| `src/widgets/drug-likeness.ts`                           | Lipinski rule compliance panel                         |
| `src/widgets/toxicity.ts`                                | Toxicity risk assessment panel                         |
| `src/widgets/identifiers.ts`                             | Identifier mapping (PubChem, ChEMBL, etc.)             |
| `src/widgets/structure2d.ts`                             | 2D structure viewer                                    |
| `src/widgets/structure3d.ts`                             | 3D structure viewer (NGL)                              |
| `src/widgets/synthon-search.ts`                          | Synthon-based searching                                |
| `src/widgets/structural-alerts.ts`                       | Structural alerts panel                                |
| `src/widgets/pharmacophore-features.ts`                  | Pharmacophore features panel                           |
| `src/widgets/molfile.ts`                                 | Molfile inspector                                      |
| `src/widgets/col-highlights.ts`                          | Per-column molecule highlight overlay                  |
| `src/panels/chem-column-property-panel.ts`               | Column property panel for chemical columns             |
| `src/panels/inchi.ts`                                    | InChI / InChI Key panel                                |
| `src/panels/structural-alerts.ts`                        | Context panel for structural alerts                    |
| `src/panels/pharmacophore-features.ts`                   | Context panel for pharmacophore features               |
| `src/analysis/activity-cliffs.ts`                        | Activity cliff detection & visualization               |
| `src/analysis/chem-similarity-viewer.ts`                 | Interactive similarity search viewer                   |
| `src/analysis/chem-diversity-viewer.ts`                  | Diversity analysis viewer                              |
| `src/analysis/chem-search-base-viewer.ts`                | Shared base for similarity/diversity viewers           |
| `src/analysis/chem-space.ts`                             | Chemical space (t-SNE/UMAP) dimensionality reduction   |
| `src/analysis/r-group-analysis.ts`                       | R-group decomposition                                  |
| `src/analysis/mpo.ts`                                    | Multi-parameter optimization (data-driven scoring)     |
| `src/analysis/deprotect.ts`                              | Deprotection (functional group removal) analysis       |
| `src/analysis/molecular-matched-pairs/`                  | MMPA: fragments, generations, differences, viewer      |
| `src/analysis/bit-birch/`                                | BitBIRCH clustering (WASM-based)                       |
| `src/mpo/`                                               | MPO profile manager, profiles view, scores viewer, context panel, profile handler/creator |
| `src/utils/reaction-enumeration/`                        | Forward-reaction library enumeration over BBs and SMARTS templates (see "Reaction enumeration" below) |
| `src/utils/reactions/`                                   | Reaction transformation/two-component reactions: browser, editor, storage, UI |
| `src/rdkit.worker.ts`                                    | Top-level RDKit worker entry                           |
| `src/rdkit-service/rdkit-service.ts`                     | Web Worker manager for parallel RDKit operations       |
| `src/rdkit-service/rdkit-service-worker-substructure.ts` | Substructure search worker                             |
| `src/rdkit-service/rdkit-service-worker-similarity.ts`   | Fingerprint computation worker                         |
| `src/open-chem/ocl-sketcher.ts`                          | OpenChemLib sketcher wrapper                           |
| `src/open-chem/sdf-importer.ts`                          | SDF file import                                        |
| `src/open-chem/ocl-service/`                             | OCL property calculations in workers                   |
| `src/docker/api.ts`                                      | Communication with chem-chem Docker container          |
| `src/descriptors/descriptors-calculation.ts`             | Docker-based descriptor computation                    |
| `src/utils/chem-common-rdkit.ts`                         | RDKit module initialization and management             |
| `src/utils/chem-common.ts`                               | Cross-module common helpers                            |
| `src/utils/chem-common-ocl.ts`                           | OpenChemLib common helpers                             |
| `src/utils/convert-notation-utils.ts`                    | SMILES/SMARTS/MolBlock conversion                      |
| `src/utils/mol-creation_rdkit.ts`                        | Molecule creation via RDKit (handles malformed input)  |
| `src/utils/sdf-utils.ts`                                 | SDF file handling                                      |
| `src/utils/most-common-subs.ts`                          | Maximum Common Substructure (MCS)                      |
| `src/utils/elemental-analysis-utils.ts`                  | Elemental analysis utilities                           |
| `src/utils/atom-index-mapper.ts`                         | Atom index mapping helpers                             |
| `src/utils/chem-atom-picker-utils.ts`                    | Atom picker utility functions                          |
| `src/utils/similarity-utils.ts`                          | Similarity computation utilities                       |
| `src/utils/mixfile.ts`                                   | Chemical mixture (mixfile) parsing                     |
| `src/file-importers/smi-importer.ts`                     | SMILES/SMARTS file import                              |
| `src/file-importers/mol2-importer.ts`                    | Tripos MOL2 format import                              |


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
- `/chem/molecules_to_inchi_key` — SMILES → InChI Key
- `/chem/inchi_to_smiles` — InChI to SMILES
- `/chem/inchi_to_inchi_key` — InChI → InChI Key

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

All functions are registered through the `PackageFunctions` class using `@grok.decorators.*`
(`func`, `app`, `init`, `autostart`, `editor`, `fileExporter`). Examples by area:

Core rendering: `chemCellRenderer`, `rdKitCellRenderer`, `rdKitReactionRenderer`,
`rdKitMixtureRenderer`

Search: `getMorganFingerprints`, `getMorganFingerprint`, `getSimilarities`, `getDiversities`,
`findSimilar`, `searchSubstructure`, `synthonSearchFunc`, `sortBySimilarity`,
`similaritySearchTopMenu`, `diversitySearchTopMenu`, `similarityMatrixTopMenu`

Analysis & viewers: `rGroupsAnalysisMenu`, `chemSpaceTopMenu`, `scaffoldTreeViewer`,
`scaffoldTreeFilter`, `substructureFilter`, `activityCliffs`, `clusterMCSTopMenu`,
`bitbirchClusteringTopMenu`, `pharmacophoreFeaturesTopMenu`, `structuralAlertsTopMenu`,
`elementalAnalysis`, `deprotect`

Apps (`#app` role, browseable in Browse tree): `reactionEnumeratorApp` (Reaction Enumerator,
"Chem | Reactions"), `transformationReactionsApp`, `twoComponentReactionsApp`, `mpoProfilesApp`,
`demoChemOverview` and other demo entries

MPO: `mpoCalculate`, `mpoWidget`, `mpoProfileEditor`, `mpoTransformFunction`, `checkJsonMpoProfile`

ChemProp ML: `trainChemprop`, `applyChemprop`, `trainModelChemprop`, `applyModelChemprop`,
`getChempropError`

Conversion / utilities: `convertMolNotation`, `convertNotation`, `recalculateCoords`,
`getRdKitModule`, `getMolFileHandler`, `canvasMol`, `drawMolecule`, `getCLogP`,
`getMolecularFormula`, `addInchisTopMenu`, `addInchisKeysTopMenu`, `freeTextToSmiles`,
`namesToSmiles`, `removeDuplicates`, `removeWaterAndSaltsTopMenu`, `validateMolecule`,
`saveAsSdf` (file exporter)

Editors: `SearchSubstructureEditor`, `ActivityCliffsEditor`, `ChemSpaceEditor`, `MMPEditor`,
`deprotectEditor` — function-call dialogs registered via `@grok.decorators.editor`

### Key constants and tags

Fingerprint types: Morgan (2048 bits, radius 2), RDKit, MACCS, AtomPair, TopologicalTorsion, Pattern

Substructure filter modes: Exact Match, Contains, Included In, Is Similar, Not Contains, Not Included In

Column tags: `.%scaffold-col`, `.%chem-scaffold-align`, `.%chem-scaffold-highlight`,
`chem-apply-filter-sync`, `.%chem-space-embedding-col`

### Dependencies

| Dependency                       | Purpose                                                      |
|----------------------------------|--------------------------------------------------------------|
| `@datagrok-libraries/chem-meta`  | RDKit JS API interfaces, molfile parsing                     |
| `@datagrok-libraries/ml`         | ML utilities, activity cliffs                                |
| `@datagrok-libraries/math`       | High-performance math                                        |
| `@datagrok-libraries/statistics` | Statistics, MPO                                              |
| `@datagrok-libraries/utils`      | BitArray, caching                                            |
| `@datagrok-libraries/gridext`    | Grid extensions used by chem viewers                         |
| `@datagrok-libraries/tutorials` | Tutorial framework integration                               |
| `openchemlib` (^7.2.3)           | Alternative renderer & sketcher                              |
| `ngl` (^2.4.0)                   | 3D structure visualization                                   |
| `js-yaml` (^4.1.1)               | YAML config parse/serialize (Reaction Enumerator)            |
| `pako` / `jszip`                 | Gzip + zip handling for Docker payloads and bundled assets   |
| `typescript-lru-cache`           | LRU caches in cell renderers                                 |
| RDKit WASM                       | Core cheminformatics engine (loaded via chem-meta + worker)  |

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

## Tests

Test files are in `src/tests/` (~40 files). Coverage includes: rendering, substructure/similarity
search and filters, MMPA, SDF / MOL2 importers, notation conversion, scaffold tree (incl. node
reordering), reaction enumeration, atom picker (2D + 3D hover, escape), pharmacophore features,
synthon search, ChemProp, vector functions, save-as-SDF, sketcher, projects/clone-layout, and
top-menu/menu test suites for chem-space, cliffs, R-groups, similarity/diversity, etc.

## Reaction enumeration

Module: `src/utils/reaction-enumeration/`. Exposed as the **Reaction Enumerator** app
(`reactionEnumeratorApp`, browsePath `Chem | Reactions`). Forward-reaction library enumeration:
takes a table of SMARTS reaction templates + a building-block library and produces a product table
with full synthesis routes.

| File                                              | Purpose                                                                    |
|---------------------------------------------------|----------------------------------------------------------------------------|
| `reaction-enumeration/config.ts`                  | YAML-backed `EnumeratorConfig` types, defaults, `cloneConfig`, parse/emit  |
| `reaction-enumeration/enumerate.ts`               | Core algorithm: round loop, template parsing, route reconstruction, `formatRoute` |
| `reaction-enumeration/filters.ts`                 | Product-spec filters (atom counts, charges, isotopes, exclusion SMARTS), `computeMolStats` |
| `reaction-enumeration/config-form.ts`             | Full-config dialog (every YAML field) — opened via "Edit full config…"     |
| `reaction-enumeration/enumerator-app.ts`          | App view: data + config inputs, debounced preview grid, run/cancel         |
| `files/enumerations/bb.csv` / `reactions.csv` / `ex_smarts.csv` | Bundled default building blocks, reactions, exclusion SMARTS  |
| `files/enumerations/aa_bb.csv` / `aa_reactions.csv` | Amino-acid library + peptide-coupling/disulfide reactions for tests        |

Key concepts:
- **depth_first** vs **breadth_first**: in depth_first round R > 1, each step combines EXACTLY ONE
  round-(R-1) product with original BBs (linear chain extension, no merging two complex products).
  Breadth-first allows any combo from rounds 0..R-1 (convergent routes possible).
- **Route format**: every step is emitted as its own `reactants>>product` segment, joined by
  `BRANCH_DELIMITER` (`--**--`) defined in `src/rendering/rdkit-reaction-renderer.ts`. The renderer's
  `parseMultiStepReaction(string): string[][]` consumes this format and the cell renderer draws each
  step with a step-arrow separator.
- Tests: `src/tests/reaction-enumeration-tests.ts` — covers parsing helpers, route formatting,
  cross-round route reconstruction, depth-first linear-extension constraint, AA library end-to-end.
