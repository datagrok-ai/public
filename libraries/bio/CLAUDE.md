# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

**@datagrok-libraries/bio** is the core shared library for all bioinformatics functionality in the Datagrok platform. It provides types, utilities, renderers, and service interfaces for working with macromolecules (peptides, DNA, RNA, HELM), 3D molecular structures (PDB/PDBQT/mmCIF), phylogenetic trees, and monomer libraries.

This is a **library** (not a package) — it has no `package.ts`, no webpack, no platform registration. It is compiled with `tsc` to `.js`/`.d.ts` and consumed by multiple packages (`Bio`, `Helm`, `Peptides`, `BiostructureViewer`, `Dendrogram`, `HitTriage`, etc.).

## Build Commands

```bash
npm install
npm run build              # grok check --soft && tsc
npm run build-all          # Builds chem-meta → js-api → gridext → utils → ml → this library
npm run lint               # eslint "./src/**/*.ts"
npm run lint-fix
npm run link-all           # Links chem-meta, datagrok-api, gridext, utils, ml
```

## Key Dependencies

- `datagrok-api` — platform API (`DG`, `grok`, `ui`)
- `@datagrok-libraries/chem-meta` — RDKit API types, chemistry metadata (MolNotation, RDModule)
- `@datagrok-libraries/utils` — BitArray, vector operations, type declarations, SVG utilities
- `@datagrok-libraries/ml` — distance metrics, dimensionality reduction
- `@datagrok-libraries/gridext` — grid extension types (IGridNeighbor)
- `@datagrok-libraries/helm-web-editor` — JSDraw2 + Pistoia HELM Web Editor types (devDependency, types only)
- `@datagrok-libraries/js-draw-lite` — JSDraw2 Lite types (devDependency, types only)
- `lru-cache`, `rxjs`, `wu`, `cash-dom`, `uuid`, `fastest-levenshtein`

## Architecture — Service Pattern

Many capabilities are defined as **interfaces here** but **implemented in packages**. Factory functions dynamically resolve implementations at runtime via `DG.Func.find()`:

| Interface | Factory function | Implementing package |
|---|---|---|
| `ISeqHelper` | `getSeqHelper()` | Bio |
| `IMonomerLibHelper` | `getMonomerLibHelper()` | Bio |
| `IHelmHelper` | `getHelmHelper()` | Helm |
| `HelmServiceBase` | `getHelmService()` | Helm |
| `IPdbHelper` | `getPdbHelper()` | BiostructureViewer |
| `IAutoDockService` | `getAutoDockService()` | Docking |
| `NglGlServiceBase` | `getNglGlService()` | BiostructureViewer |
| `PhylocanvasGlServiceBase` | `getPhylocanvasGlService()` | PhyloTreeViewer |
| `ITreeHelper` | `getTreeHelper()` | Dendrogram |
| `IDendrogramService` | `getDendrogramService()` | Dendrogram |
| `RDModule` | `getRdKitModule()` | Chem |

This decouples the library from concrete package implementations, allowing packages to be installed independently.

## Source Structure

### Top-Level Modules (`src/`)

| File | Purpose |
|---|---|
| `aminoacids.ts` | `Aminoacids` class (names, codes, semantic types) + `AminoacidsPalettes` (Lesk, GrokGroups, RasMol color schemes) |
| `nucleotides.ts` | `Nucleotides` class (names, codes) + `NucleotidesPalettes` (Chromatogram colors) |
| `seq-palettes.ts` | `SeqPalette` interface + `SeqPaletteBase` base class with shared color palette definitions |
| `sequence-encoder.ts` | `AlignedSequenceEncoder` — categorical/numerical encoding of amino acid sequences (Wimley-White, categorial scales) |
| `unknown.ts` | `UnknownSeqPalettes` — hash-based color assignment for non-standard monomers, reads custom colors from monomer library |
| `utils.ts` | `getColumnSeparator()` — detects separator character in sequence columns |

### `types/` — Core Type Definitions

| File | Purpose |
|---|---|
| `monomer-library.ts` | **Key file.** `Monomer` type (HELM JSON schema), `IMonomerLib` / `IMonomerLibBase` (monomer library interface), `IMonomerLibHelper` (library management singleton), `IMonomerLibProvider` (pluggable library sources — files, DBs, etc.), `MonomerSet` / `MonomerSetPlaceholder` |
| `renderer.ts` | `IRenderer` interface — `onRendered` observable + `invalidate()` + `awaitRendered()` for test synchronization |
| `input.ts` | `InputColumnBase` — extended column input with `setColumnInputTable()`, augments `ui.input` namespace |
| `dojo.ts` | Dojo Toolkit module declaration (for HELM Web Editor) |
| `ngl.ts` | NGL.js library type declarations (Stage, Component, Viewer, Colormaker, etc.) |
| `index.ts` | Re-exports from `monomer-library.ts` |

### `helm/` — HELM Notation Types and Helpers

| File | Purpose |
|---|---|
| `types.ts` | Re-exports all HELM types from `@datagrok-libraries/helm-web-editor` and `js-draw-lite`: `HelmType`, `HelmAtom`, `HelmBond`, `HelmMol`, `HelmEditor`, `PolymerType`, `MonomerType`, `IHelmWebEditor`, `ISeqMonomer`, `MonomerSetType`, etc. |
| `consts.ts` | Re-exports `HelmTypes`, `MonomerTypes`, `PolymerTypes`, `MonomerNumberingTypes`, `HelmTabKeys` |
| `helm-helper.ts` | `IHelmHelper` interface (parse, removeGaps, getMolfiles, createHelmInput, createWebEditorApp), `HelmInputBase` abstract class, `HelmNotSupportedError`, `getHelmHelper()` factory, `getMonomerHandleArgs()` helper. Augments `ui.input` with `helmAsync()` |
| `utils.ts` | `cleanupHelmSymbol()` — strips brackets from HELM monomer symbols |

### `molecule/` — 2D Molecule Column Handling

| File | Purpose |
|---|---|
| `types.ts` | `MoleculeBase` data class, `MoleculeBuildDataFunc` type |
| `molecule-build-data.ts` | `buildDataMolV2000` / `buildDataMolV3000` — parsers for MOL V2000/V3000 |
| `molecule-units-handler.ts` | `MoleculeUnitsHandler` — units handler for `Molecule` semtype columns (molV2000/molV3000), `getAsPdb()` conversion |

### `molecule-3d/` — 3D Molecular Structure Column Handling

| File | Purpose |
|---|---|
| `types.ts` | `Molecule3DBase` data class, `Molecule3DBuildDataFunc` type |
| `molecule-3d-build-data.ts` | `buildDataPdb` / `buildDataPdbqt` / `buildDataMmcif` — parsers for PDB/PDBQT/mmCIF |
| `molecule-3d-units-handler.ts` | `Molecule3DUnitsHandler` — units handler for `Molecule3D` semtype columns (pdb/pdbqt/mmcif), `getAsPdb()` conversion |
| `index.ts` | Barrel re-export |

### `monomer-works/` — Sequence-to-Molfile Conversion Engine

The core engine for converting macromolecule sequences into V3000 molfiles.

| File | Purpose |
|---|---|
| `types.ts` | Core types: `MonomerGraph` (atoms, bonds, R-groups), `MolGraph` (assembled molecule graph), `MonomerMapValue` (atom/bond indices per position), `MonomerSequenceDict`, `SeqToMolfileWorkerInput/Output`, `LibSettings` |
| `consts.ts` | V2K/V3K parsing tokens, precision factors, canonical nucleotide component symbols (ribose, deoxyribose, phosphate) |
| `monomer-works.ts` | `MonomerWorks` class — wraps `IMonomerLib`, provides `getMonomerMolfile()`. `helmTypeToPolymerType()` converter |
| `to-atomic-level.ts` | **Core engine.** `seqColToMolFileColumn()` — converts macromolecule column to molfile column. Handles molfile parsing, V2K→V3K conversion, spatial adjustment (rotation, flipping, coordinate shifting), R-group capping, stereo-center handling, and final V3K assembly |
| `to-atomic-level-utils.ts` | Molfile assembly engine: `buildMolfileFromSeq()` — chains monomers into a complete V3K molfile. Handles peptide bonds, nucleotide assembly (sugar+phosphate+base), terminal capping |
| `seq-to-molfile.ts` | Parallel orchestration: `seqToMolfile()` — spawns Web Workers for parallel conversion. `getHighlightMoleculeData()` — builds per-monomer atom/bond coloring for highlights |
| `seq-to-molfile-worker.ts` | Web Worker entry point — receives sequence chunks, calls `buildMolfileFromSeq()`, posts results back |
| `monomer-hover.ts` | `addMonomerHoverLink()` — creates hover links between sequence cells and molecule cells with LRU-cached monomer maps. `executeHoverLinks()` / `getHoverLinks()` / `addHoverLink()` |
| `monomer-utils.ts` | Analytics: `encodeMonomerSymbol()` (Unicode encoding), `getMonomerMolfiles()`, `sdfToJson()` (SDF→HELM library conversion), `calculatePositionChemSimilarity()`, `monomerPairwiseSimilarity()`, `getSubstitutionMatrix()` |
| `lib-settings.ts` | `loadLibSettings()` / `saveLibSettings()` — user-specific monomer library preferences with chunked storage |
| `utils.ts` | Small shared helpers: `getAtomicLevelColName()`, `helmTypeToNotation()`, `hexToPercentRgba()` |

**Data flow for sequence → molfile conversion:**
1. `seqColToMolFileColumn()` in `to-atomic-level.ts` is the entry point
2. It parses monomer molfiles from the library into `MonomerGraph` objects
3. `seqToMolfile()` in `seq-to-molfile.ts` spawns Web Workers
4. Each worker (`seq-to-molfile-worker.ts`) calls `buildMolfileFromSeq()` from `to-atomic-level-utils.ts`

### `utils/` — Core Utilities

#### Macromolecule subsystem (`utils/macromolecule/`)

| File | Purpose |
|---|---|
| `types.ts` | **Key file.** `ISplitted` (split sequence interface), `INotationProvider` (notation-specific behavior), `SplitterFunc` type, `CandidateType` (alphabet detection), `SeqColStats` / `SeqColStatsCached` |
| `consts.ts` | **Key file.** `NOTATION` enum (FASTA, SEPARATOR, HELM), `ALIGNMENT` enum, `ALPHABET` enum (DNA, RNA, PT, UN), `TAGS` (column tag names for units, aligned, alphabet, separator, region), `GAP_SYMBOL`, `ALPHABET_CHARS`, regex patterns for FASTA/HELM parsing |
| `utils.ts` | **Largest utility file.** All `Splitter` implementations: `SplitterBase`, `SplitterFasta`, `SplitterHelm`, `SplitterBiln`. Splitter factory functions. Alphabet detection via cosine similarity. Palette selection. Monomer abbreviation. `getJoiner()`, `getAlphabetSimilarity()`, `getStats()`, `candidateStats()` |
| `seq-handler.ts` | `ISeqHandler` interface — central abstraction for a macromolecule column: notation, alphabet, separator, splitter/joiner, HELM conversion, region extraction, distance function selection. `SeqValueBase` — wraps a single row value with `getOriginal()`, `getCanonical()`, `helm` |
| `scoring.ts` | `ScoringMethod` enum + `calculateIdentityScoring()` / `calculateChemSimilarityScoring()` — sequence scoring vs reference |
| `alignment.ts` | `pairwiseAlignmentWithEmptyPositions()` — Needleman-Wunsch alignment with gap penalties |
| `monomers.ts` | Precomputed SMILES/fingerprints for 20 amino acids + 5 nucleotides. `calcMonomerFps()`, `matchMonomerToNatural()` |
| `index.ts` | Barrel re-exports for the macromolecule module |

#### Cell Rendering

| File | Purpose |
|---|---|
| `cell-renderer-consts.ts` | Temp/tag name constants for renderer settings (font size, color coding, max monomer length, multiline) |
| `cell-renderer.ts` | Low-level `drawMonomer()` function — paints individual monomers on canvas with MSA/classic styles, transparency, color coding |
| `cell-renderer-monomer-placer.ts` | `MonomerPlacer` — the core rendering engine for macromolecule sequences in grid cells. Computes per-position widths, handles single-line/multiline layout, MSA alignment, reference sequence comparison, position scrolling, hit-testing for hover/click, monomer tooltips |
| `cell-renderer-back-base.ts` | `CellRendererBackBase` — lifecycle management base class (dataframe subscriptions, dirty tracking, `IRenderer` implementation) |
| `cell-renderer-async-base.ts` | `GridCellRendererBackAsyncBase` — async rendering with priority queue, LRU image cache, timeout/retry, stale task sweeping |

#### Other Utilities

| File | Purpose |
|---|---|
| `seq-helper.ts` | `ISeqHelper` interface + `getSeqHelper()` factory. `ISeqHelper` is the primary entry point for sequence operations: `getSeqHandler()`, `setNotationProvider()`, `getMonomersList()`, `helmToSmiles()`, `helmToMolfile()` |
| `const.ts` | HELM monomer library JSON schema field names (`HELM_REQUIRED_FIELD`, `HELM_OPTIONAL_FIELDS`, `HELM_RGROUP_FIELDS`), SDF-to-JSON mapping, dummy monomer template, encoding ranges |
| `splitter.ts` | `joinDataFrames()` — splits aligned sequences into per-position columns |
| `composition-table.ts` | `getCompositionTable()` — builds HTML monomer composition bar chart |
| `fasta-handler.ts` | `FastaFileHandler` — parses FASTA files into DataFrames |
| `generator.ts` | Test data generators for synthetic macromolecule columns |
| `sequence-position-scroller.ts` | `SequencePositionScroller` — interactive MSA header with WebLogo tracks, conservation scoring, position slider, composition tooltips (~1600 lines) |
| `data-provider.ts` | `getDataProviderFuncs()` — discovers registered data provider functions by semtype |
| `docker.ts` | `awaitStatus()` — polls Docker container status until target state |
| `syncer.ts` | `Syncer` — serializes concurrent async operations into sequential execution |
| `units-handler-base.ts` | `UnitsHandlerBase<S, D>` — generic base class for semtype column units handling with lazy data parsing |
| `err-info.ts` | `errInfo()` / `toErrorString()` — error message extraction from various error types |
| `logger.ts` | `ILogger` / `BioLogger` — structured logging with prefix and debug suppression |
| `monomer-ui.ts` | Interfaces for monomer editor UI (`ICreateMonomerForm`, `IMonomerManager`) and DG property schema for monomer input fields |

### `viewers/` — Viewer Interface Definitions

All viewer types follow `IXxxViewer` + property defaults pattern. Implementations live in their respective packages.

| File | Purpose |
|---|---|
| `viewer.ts` | `IBioViewer` base interface (extends `IRenderer`) |
| `web-logo.ts` | `IWebLogoViewer` — sequence logo viewer properties and interface |
| `helm-service.ts` | `HelmServiceBase` — off-screen HELM rendering service abstraction |
| `molecule3d.ts` | `Molecule3DData` type, `DockingRole` enum, `IMolecule3DViewer` interface |
| `molstar-viewer.ts` | `IMolstarViewer` — Mol* biostructure viewer properties |
| `ngl-gl-service.ts` | `NglGlServiceBase` — NGL cell-level rendering service |
| `ngl-gl-viewer.ts` | `INglGlViewer` — full NGL 3D viewer properties |
| `phylocanvas-gl-viewer.ts` | `IPhylocanvasGlViewer` + `PhylocanvasGlServiceBase` — phylogenetic tree viewer/service |
| `vd-regions.ts` | `IVdRegionsViewer` — antibody V-Domain region viewer (CDR/FR with IMGT numbering) |
| `biotrack.ts` | `IBiotrackViewer` — biotrack viewer interface |

### `pdb/` — PDB/PDBQT File Format Handling

| File | Purpose |
|---|---|
| `types.ts` | `BiostructureData` type (binary/text with extension), `BiostructureDataJson` namespace for serialization |
| `pdb-helper.ts` | `IPdbHelper` interface (parse, convert PDB/PDBQT/mol), `PdbAtomDataFrame` typed DataFrame |
| `auto-dock-service.ts` | `IAutoDockService` interface (docking operations, Docker container management) |
| `index.ts` | `PdbTag` constant export |
| `format/` | PDB/PDBQT atom record parsing classes: `AtomRecordBase` → `PdbAtomRecord` / `PdbqtAtomRecord`, with coordinate parsing, chain/residue sorting, TER records |

### `trees/` — Phylogenetic Tree Handling

| File | Purpose |
|---|---|
| `types.ts` | `NodeType` (name, branch_length, children), `NodeCuttedType`, `ClusterMatrix` |
| `consts.ts` | Newick tags, tree color palette, distance metrics (Euclidean, Manhattan), linkage methods (single, complete, average, Ward, etc.) |
| `phylocanvas.ts` | PhylocanvasGL integration: tree layouts, node shapes, Newick parser |
| `tree-helper.ts` | `ITreeHelper` — comprehensive tree manipulation (Newick↔DataFrame, clustering, filtering, cutting, grid reordering) |
| `dendrogram.ts` | `IDendrogramService` — grid-adjacent dendrogram visualization |
| `utils.ts` | `isLeaf()` helper |
| `index.ts` | Barrel re-export |

### `tests/` — Library-Level Tests

| File | Purpose |
|---|---|
| `monomer-lib-tests.ts` | `testExpectedMonomerLib()` — validates IMonomerLib against expected polymer type counts |
| `palettes-tests.ts` | Tests for nucleotide and amino acid palette instantiation |

### Other

| File | Purpose |
|---|---|
| `chem/rdkit-module.ts` | `getRdKitModule()` — dynamically loads RDKit WASM from Chem package |
| `ntseq/ntseq.js` | High-performance nucleotide sequence library (4-bit encoding, complement, translation, alignment/mapping) |
| `substructure-filter/bio-substructure-filter-types.ts` | `IBioSubstructureFilter` abstract class — framework for FASTA/HELM substructure filtering |

## Key Concepts

### Semantic Types & Notations

Macromolecule columns have `semType = 'Macromolecule'` with metadata tags:
- `units` — notation: `fasta`, `separator`, `helm`
- `alphabet` — `PT` (peptide), `DNA`, `RNA`, `UN` (unknown)
- `separator` — monomer delimiter for separator notation (`.`, `-`, ` `)
- `aligned` — `SEQ` (unaligned), `SEQ.MSA` (multiple sequence alignment)

### Monomer Library System

Monomers follow the Pistoia HELM JSON schema. The library system has three layers:
1. **`IMonomerLibProvider`** — pluggable data sources (files, databases)
2. **`IMonomerLibHelper`** — singleton manager, loads/reloads libraries from configured providers
3. **`IMonomerLib`** — the in-memory monomer dictionary, queryable by polymer type + symbol

### Sequence Splitting

Sequences are split into per-position monomers via `SplitterFunc` implementations:
- `SplitterFasta` — character-level splitting for FASTA notation
- `SplitterBase` — separator-based splitting
- `SplitterHelm` — HELM parsing with graph/connection info
- `SplitterBiln` — BILN notation with cyclization marks

The split result implements `ISplitted` with `getOriginal(pos)`, `getCanonical(pos)`, `isGap(pos)`, `length`.

### Cell Rendering Pipeline

Macromolecule grid cell rendering uses a layered async architecture:
1. `CellRendererBackBase` — lifecycle management, dirty tracking
2. `GridCellRendererBackAsyncBase` — async rendering with LRU cache
3. `MonomerPlacer` — layout engine (single/multiline, position widths, scrolling)
4. `drawMonomer()` — low-level canvas drawing with palette colors

## Quick Lookups

| Looking for... | Check first |
|---|---|
| Monomer library types (`Monomer`, `IMonomerLib`) | `src/types/monomer-library.ts` |
| Monomer library management (`IMonomerLibHelper`) | `src/types/monomer-library.ts` |
| HELM types (`HelmType`, `HelmAtom`, `HelmMol`) | `src/helm/types.ts` |
| HELM helper interface (`IHelmHelper`) | `src/helm/helm-helper.ts` |
| Sequence helper interface (`ISeqHelper`) | `src/utils/seq-helper.ts` |
| Macromolecule handler interface (`ISeqHandler`) | `src/utils/macromolecule/seq-handler.ts` |
| Notation/alphabet/tag constants | `src/utils/macromolecule/consts.ts` |
| Sequence splitter implementations | `src/utils/macromolecule/utils.ts` |
| Macromolecule type system (`ISplitted`, `SplitterFunc`) | `src/utils/macromolecule/types.ts` |
| Amino acid palettes & names | `src/aminoacids.ts` |
| Nucleotide palettes & names | `src/nucleotides.ts` |
| Non-standard monomer palettes | `src/unknown.ts` |
| Cell renderer for sequences | `src/utils/cell-renderer-monomer-placer.ts` |
| Low-level monomer drawing | `src/utils/cell-renderer.ts` |
| Async cell renderer base classes | `src/utils/cell-renderer-async-base.ts` |
| Sequence-to-molfile conversion | `src/monomer-works/to-atomic-level.ts` |
| Molfile assembly from monomers | `src/monomer-works/to-atomic-level-utils.ts` |
| Parallel molfile conversion | `src/monomer-works/seq-to-molfile.ts` |
| Monomer hover highlighting | `src/monomer-works/monomer-hover.ts` |
| Monomer similarity/scoring matrices | `src/monomer-works/monomer-utils.ts` |
| HELM JSON schema constants | `src/utils/const.ts` |
| PDB/PDBQT parsing | `src/pdb/format/` |
| PDB helper interface | `src/pdb/pdb-helper.ts` |
| AutoDock service interface | `src/pdb/auto-dock-service.ts` |
| 3D structure data types | `src/viewers/molecule3d.ts` + `src/pdb/types.ts` |
| Tree types & Newick parsing | `src/trees/phylocanvas.ts` |
| Tree helper interface | `src/trees/tree-helper.ts` |
| Dendrogram service interface | `src/trees/dendrogram.ts` |
| WebLogo viewer interface | `src/viewers/web-logo.ts` |
| Molstar viewer interface | `src/viewers/molstar-viewer.ts` |
| NGL viewer interface | `src/viewers/ngl-gl-viewer.ts` |
| FASTA file parsing | `src/utils/fasta-handler.ts` |
| Sequence alignment (Needleman-Wunsch) | `src/utils/macromolecule/alignment.ts` |
| Sequence scoring (identity/similarity) | `src/utils/macromolecule/scoring.ts` |
| Natural monomer fingerprints | `src/utils/macromolecule/monomers.ts` |
| MSA position scroller/header | `src/utils/sequence-position-scroller.ts` |
| Substructure filter framework | `src/substructure-filter/bio-substructure-filter-types.ts` |
| Nucleotide sequence library (fast) | `src/ntseq/ntseq.js` |
| RDKit module accessor | `src/chem/rdkit-module.ts` |
| Test data generation | `src/utils/generator.ts` |
| Error handling utilities | `src/utils/err-info.ts` |
| Docker container polling | `src/utils/docker.ts` |
