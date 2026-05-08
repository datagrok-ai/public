# CLAUDE.md

## Overview

**@datagrok-libraries/bio** is the core shared library for all bioinformatics functionality in the Datagrok platform. It provides types, utilities, renderers, and service interfaces for working with macromolecules (peptides, DNA, RNA, HELM), 3D molecular structures (PDB/PDBQT/mmCIF), phylogenetic trees, and monomer libraries.

This is a **library** (not a package) — it has no `package.ts`, no webpack, no platform registration. It is compiled with `tsc` to `.js`/`.d.ts` and consumed by multiple packages (`Bio`, `Helm`, `Peptides`, `BiostructureViewer`, `Dendrogram`, `HitTriage`, etc.).

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
| `types.ts` | Re-exports HELM types from `@datagrok-libraries/helm-web-editor` and `js-draw-lite`: `HelmType`, `HelmAtom` / `Atom` / `IJsAtom`, `HelmBond` / `Bond`, `HelmMol` / `Mol`, `HelmEditor` / `Editor`, `MonomerExplorer`, `PolymerType`, `MonomerType`, `IHelmWebEditor`, `ISeqMonomer`, `MonomerSetType`, `WebEditorRGroups`, `MonomersFuncs`, `TabDescType`, `IBio` / `IHelmBio`, `HelmString`, `Point`, browser-side singletons `HweWindow` / `ScilModuleType` / `JSDraw2ModuleType` / `OrgType` / `DojoType` / `DojoxType` |
| `consts.ts` | Re-exports `HelmTypes`, `MonomerTypes`, `PolymerTypes`, `MonomerNumberingTypes`, `HelmTabKeys` |
| `helm-helper.ts` | `IHelmHelper` interface — `parse()`, `removeGaps()`, `getMolfiles()`, `createHelmInput()`, `createHelmWebEditor()`, `createWebEditorApp()`, `getHoveredAtom()`, plus monomer-funcs hooks (`originalMonomersFuncs` getter, `buildMonomersFuncsFromLib()`, `overrideMonomersFuncs()` / `revertOriginalMonomersFuncs()`) and `seqHelper` getter. `HelmInputBase` abstract class. `HelmNotSupportedError` (+ `HelmNotSupportedErrorType`). `getHelmHelper()` factory. `getMonomerHandleArgs()` helper. `IHelmInputInitOptions`, `HelmConvertRes` exported. Augments `ui.input` with `helmAsync()` |
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
| `types.ts` | Core types: `Atoms` / `Bonds` / `MonomerMetadata` (parsed-molfile data), `MolGraph` (assembled molecule graph) + `MonomerMolGraphMap` + `getMolGraph` / `hasMolGraph` / `setMolGraph` helpers, `LibMonomerKey`, `MonomerMap` / `MonomerMapValue` / `MolfileWithMap`, `LoopVariables` / `LoopConstants` / `NumberWrapper`, `Point`, `ITypedArray`, `UserLibSettings`, `NucleotideRole` enum (`SUGAR` / `BASE` / `PHOSPHATE` / `TERMINAL_5P` / `TERMINAL_3P`), `SeqToMolfileWorkerData` / `SeqToMolfileWorkerRes` |
| `consts.ts` | `monomerWorksConsts` object — V2K/V3K parsing tokens, `PRECISION_FACTOR`, canonical nucleotide component keys (`RIBOSE` / `DEOXYRIBOSE` / `PHOSPHATE`), `OXYGEN` / `HYDROGEN` element symbols |
| `monomer-works.ts` | `MonomerWorks` class — wraps `IMonomerLib`, provides `getMonomerMolfile()`. `helmTypeToPolymerType()` converter. `IMolfileConverter` interface |
| `to-atomic-level.ts` | **Core engine.** `_toAtomicLevel()` — converts macromolecule column to V3000 molfile column. HELM RNA path (triples mode) keeps per-position sugar/base/phosphate (via `getMonomerSequencesArray`, `getMonomersDictFromLib`, `buildRolesForHelmRna`); legacy path uses default ribose/deoxyribose/phosphate. Handles V2K→V3K conversion (`convertMolfileToV3K`, `fixV2000MolfileRAtomLines`), atom/bond block parsing (`parseAtomBlock` is internal; `parseBondBlock`, `parseAtomAndBondCounts`, `parseCapGroups`, `parseCapGroupIdxMap` / `parseCapGroupIdxMapV2K` / `parseCapGroupIdxMapV3K` are exported), R-group capping, stereo-center handling, V3K serialization (`convertMolGraphToMolfileV3K`), capped-monomer rendering (`getSymbolToCappedMolfileMap`, `capPeptideMonomer`). The geometric pipeline (`adjustSugarMonomerGraph`, `adjustPhosphateMonomerGraph`, `adjustBaseMonomerGraph`, `adjustPeptideMonomerGraph`) places R-attached atoms (`terminalNodes`) on the OX axis with R3 / branch up; the sugar branch is overridden so the base is placed above the topmost atom for non-canonical sugars (e.g. LNA's 2,4-bridge). `setShiftsAndTerminalNodes` injects an O at the R1 cap position when the cap is `H` so phosphorothioate-style linkers (sp, en, …) keep a real C-O-P bridging atom |
| `to-atomic-level-utils.ts` | Molfile assembly engine: `monomerSeqToMolfile()` — chains monomers into a complete V3K molfile (peptide bonds, nucleotide assembly: sugar+base+phosphate, terminal capping). `runTriplesAssembly()` (internal) handles HELM RNA per-row triples including 5'/3' terminal modifiers (Chol / GalNAc / Bio) and missing-trailing-phosphate. `getFormattedMonomerLib()` (lib → per-symbol object map). `keepPrecision()` (coordinate rounding helper) |
| `seq-to-molfile.ts` | `seqToMolFileWorker()` — orchestrates the molfile column build (single-threaded today; the `Worker` plumbing is wired but disabled — runs synchronously per row). `getMolHighlight()` — builds an `ISubstruct` from per-monomer atom/bond ranges for hover highlights. `SeqToMolfileResult` type |
| `seq-to-molfile-worker.ts` | Web Worker entry point — receives `SeqToMolfileWorkerData`, calls `monomerSeqToMolfile()`, posts `SeqToMolfileWorkerRes` back |
| `monomer-hover.ts` | `addMonomerHoverLink()` — wires hover from a sequence cell to a molecule cell with LRU-cached monomer maps. `execMonomerHoverLinks()` / `getMonomerHoverLinks()`. `MonomerHoverLinksTemp` column-tag constant |
| `monomer-utils.ts` | Analytics: `encodeMonomers()` (Unicode encoding of monomer symbols for cosine-distance use), `getMolfilesFromSeq()` / `getMolfilesFromSingleSeq()` (per-position molfile assembly), `createMomomersMolDict()`, `createJsonMonomerLibFromSdf()` (SDF table → HELM JSON library) |
| `lib-settings.ts` | `getUserLibSettings()` / `setUserLibSettings()` — user-specific monomer library preferences with chunked storage |
| `utils.ts` | Small shared helpers: `getUnusedColName()`, `getMolColName()`, `alphabetToPolymerType()`, `hexToPercentRgb()`. `MonomerHoverLink` type |

**Data flow for sequence → molfile conversion:**
1. `_toAtomicLevel()` in `to-atomic-level.ts` is the entry point. HELM RNA columns stay in triples mode (per-position sugar / base / phosphate); other notations convert to separator first.
2. It parses monomer molfiles from the library into `MolGraph` objects (`getMonomersDictFromLib` → `getMolGraph` per symbol).
3. `seqToMolFileWorker()` in `seq-to-molfile.ts` iterates rows.
4. Each row calls `monomerSeqToMolfile()` from `to-atomic-level-utils.ts`, which threads the row through the geometric assembly loop. The chirality engine (`Chem:convertToV3KViaOCL`) is applied to the resulting column to add `STEABS` blocks.

### `utils/` — Core Utilities

#### Macromolecule subsystem (`utils/macromolecule/`)

| File | Purpose |
|---|---|
| `types.ts` | **Key file.** `ISeqSplitted` (split sequence interface) + `SeqSplittedBase`, `INotationProvider` (notation-specific behavior) + `NotationProviderBase`, `SplitterFunc` type, `IMonomerCanonicalizer`, `CandidateType` / `CandidateSimType` (alphabet detection), `MonomerFreqs`, `SeqColStats`, graph types `ISeqConnection` / `ISeqGraphInfo` |
| `consts.ts` | **Key file.** `NOTATION` enum (FASTA, SEPARATOR, HELM, BILN, CUSTOM), `ALIGNMENT` enum, `ALPHABET` enum (DNA, RNA, PT, UN), `TAGS` enum (column tag names — `units`, `aligned`, `alphabet`, `separator`, `region`, plus also exported as `BioTags`), `Alphabets` namespace (per-alphabet character sets), `GAP_SYMBOL`, `GapOriginals`, `candidateAlphabets`, `monomerRe` / `helmRe` / `helmPp1Re`, `MONOMER_MOTIF_SPLITTER`, `MONOMER_CANONICALIZER_FUNC_TAG` / `MONOMER_CANONICALIZER_TEMP`, `NOTATION_PROVIDER_CONSTRUCTOR_ROLE`, `positionSeparator` |
| `utils.ts` | **Largest utility file.** Splitter implementations: `StringListSeqSplitted` (separator), `FastaSimpleSeqSplitted`, `HelmSplitted` (with graph info), `BilnSeqSplitted` (cycles/connections). Splitter factory funcs: `splitterAsFasta`, `splitterAsFastaSimple`, `splitterAsHelm`, `splitterAsBiln`, `getSplitterWithSeparator`, `getSplitter`. Alphabet detection: `detectAlphabet`, `detectHelmAlphabet`, `getAlphabet`, `getAlphabetSimilarity`, `getStatsForCol`, `pickUpPalette` / `pickUpSeqCol` / `getPaletteByType`. Monomer abbreviation (`monomerToShort`). HELM helpers: `polymerTypeToHelmType`, `RNA_HELM_TRIPLET_MONOMER_REG`, `RNA_HELM_TERMINAL_PHOSPHATELESS_MONOMER_REG` |
| `seq-handler.ts` | `ISeqHandler` interface — central abstraction for a macromolecule column: notation, alphabet, separator, splitter/joiner, HELM conversion, region extraction, distance function selection. `SeqValueBase` — wraps a single row value with `getOriginal()`, `getCanonical()`, `helm` |
| `scoring.ts` | `ScoringMethod` enum + `calculateIdentityScoring()` / `calculateChemSimilarityScoring()` — sequence scoring vs reference |
| `alignment.ts` | `pairwiseAlignmentWithEmptyPositions()` — Needleman-Wunsch alignment with gap penalties |
| `monomers.ts` | Precomputed SMILES + fingerprints for the 20 amino acids and 5 nucleotides: `naturalMonomers`, `naturalMonomerFps`, `getMorganFingerprint()`, `mostSimilarNaturalAnalog()` |
| `annotations.ts` | Sequence-annotation type system: `AnnotationVisualType` / `AnnotationCategory` / `LiabilitySeverity` enums, `SeqAnnotation` / `SeqAnnotationHit` interfaces, `RowAnnotationData` row payload, default colors/opacity constants used by the annotation cell-renderer overlay |
| `numbering-schemes.ts` | Antibody numbering scheme region definitions: `NumberingScheme` / `ChainType` enums, `SchemeRegionDef`, and `IMGT_REGIONS` / `KABAT_REGIONS` / `CHOTHIA_REGIONS` / `AHO_REGIONS` plus `SCHEME_REGIONS` lookup |
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
| `seq-helper.ts` | `ISeqHelper` interface + `getSeqHelper()` factory. `ISeqHelper` is the primary entry point for sequence operations: `helmToAtomicLevel()` / `helmToAtomicLevelSingle()` (column-wide and one-shot HELM → V3K molfile), `getSeqHandler()`, `getSeqMonomers()`, `setUnitsToFastaColumn()` / `setUnitsToSeparatorColumn()` / `setUnitsToHelmColumn()`, `getHelmToMolfileConverter()`. Also exports `ToAtomicLevelRes`, `IHelmToMolfileConverter` |
| `const.ts` | HELM monomer library JSON schema field names (`HELM_REQUIRED_FIELD`, `HELM_OPTIONAL_FIELDS`, `HELM_RGROUP_FIELDS`), SDF-to-JSON mapping, dummy monomer template, encoding ranges |
| `splitter.ts` | `splitAlignedSequences()` — splits aligned sequences into per-position columns |
| `composition-table.ts` | `getCompositionTable()` — builds HTML monomer composition bar chart |
| `fasta-handler.ts` | `FastaFileHandler` — parses FASTA files into DataFrames |
| `generator.ts` | Test data generators for synthetic macromolecule columns |
| `sequence-position-scroller.ts` | MSA header / scroll bar engine: `MSAScrollingHeader`, abstract `MSAHeaderTrack`, plus `WebLogoTrack` and `ConservationTrack` track implementations. Renders sequence-logo + conservation tracks above an MSA grid with a scrollable position slider and composition tooltips |
| `annotation-track.ts` | `AnnotationTrack` — `MSAHeaderTrack` subclass that overlays region/liability annotations on the MSA scroller |
| `cell-renderer-annotations.ts` | `AnnotationRenderer` — paints region rectangles and liability underlines on the macromolecule cell renderer |
| `macromolecule-highlight.ts` | Cross-package monomer highlight bus: `fireMacromoleculeHighlight()`, `setMacromoleculeMonomerHighlight()` / `setMacromoleculeMonomerHighlights()`, `clear*` siblings, `MACROMOLECULE_HIGHLIGHT_EVENT_ID` / `…_TEMP` / `…_ALPHA` constants, `MacromoleculeHighlightEntry` / `EventArgs` / `Colors` / `Target` / `RowSpec` types, `macromoleculeHighlightColorToCss()` |
| `monomer-selection-dialog.ts` | `MonomerSelectionWidget` + `showMonomerSelectionDialog()` — picker dialog for choosing monomer symbols. `parseMonomerSymbolList()` parses comma/space-separated symbol input |
| `sequence-column-input.ts` | `ISequenceColumnInput` / `SequenceColumnInputBase` + `createSequenceColumnInput()` — DG input that filters its choices to `Macromolecule` columns and exposes notation/alphabet metadata |
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
| `index.ts` | `TAGS` enum (currently just `TAGS.PDB` = `'.pdb'`) — column-tag identifier used for PDB-bearing data |
| `format/` | PDB/PDBQT atom record parsing classes. `types-base.ts` defines the abstract `LineBase` and `AtomBase` (with shared coord/residue/chain parsing); `types-pdb.ts` and `types-pdbqt.ts` carry the `Pdb*` / `Pdbqt*` line implementations including TER records and chain/residue sorting; `types.ts` re-exports the public surface |

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
| `monomer-lib-tests.ts` | `expectMonomerLib()` — validates `IMonomerLib` against expected polymer type counts |
| `palettes-tests.ts` | `_testPaletteN()` / `_testPaletteAA()` — palette instantiation tests for nucleotide and amino acid palettes |

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

Sequences are split into per-position monomers via `SplitterFunc` factories returning `ISeqSplitted`:
- `splitterAsFastaSimple` / `splitterAsFasta` → `FastaSimpleSeqSplitted` — character-level splitting for FASTA notation (the multichar variant handles `[multi]` brackets)
- `getSplitterWithSeparator` → `StringListSeqSplitted` — separator-based splitting
- `splitterAsHelm` → `HelmSplitted` — HELM parsing with graph/connection info (`ISeqGraphInfo` exposes per-monomer `polymerTypes`, `cycles`, etc.)
- `splitterAsBiln` → `BilnSeqSplitted` — BILN notation with cyclization marks
- `getSplitter(notation, sep)` is the dispatcher used everywhere else

The split result implements `ISeqSplitted` with `getOriginal(pos)`, `getCanonical(pos)`, `isGap(pos)`, `length`, plus `graphInfo` on HELM. `SeqSplittedBase` is the abstract base for custom implementations.

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
| Macromolecule type system (`ISeqSplitted`, `SplitterFunc`) | `src/utils/macromolecule/types.ts` |
| Antibody numbering scheme region maps | `src/utils/macromolecule/numbering-schemes.ts` |
| Sequence-annotation type system | `src/utils/macromolecule/annotations.ts` |
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
| MSA annotation track | `src/utils/annotation-track.ts` |
| Annotation cell-renderer overlay | `src/utils/cell-renderer-annotations.ts` |
| Cross-package monomer highlight bus | `src/utils/macromolecule-highlight.ts` |
| Monomer picker dialog | `src/utils/monomer-selection-dialog.ts` |
| Sequence column input (`Macromolecule`-only) | `src/utils/sequence-column-input.ts` |
| Substructure filter framework | `src/substructure-filter/bio-substructure-filter-types.ts` |
| Nucleotide sequence library (fast) | `src/ntseq/ntseq.js` |
| RDKit module accessor | `src/chem/rdkit-module.ts` |
| Test data generation | `src/utils/generator.ts` |
| Error handling utilities | `src/utils/err-info.ts` |
| Docker container polling | `src/utils/docker.ts` |
