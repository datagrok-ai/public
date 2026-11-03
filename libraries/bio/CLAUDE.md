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
| `types.ts` | Re-exports HELM types from `@datagrok-libraries/helm-web-editor` and `js-draw-lite`: `HelmType`, `HelmAtom` / `Atom`, `Bond`, `HelmMol` / `Mol`, `HelmEditor` / `Editor`, `PolymerType`, `MonomerType`, `IHelmWebEditor`, `ISeqMonomer`, `MonomerSetType`, `WebEditorRGroups`, `MonomersFuncs`, `TabDescType`, `IHelmBio`, `HelmString`, `Point`, browser-side singletons `HweWindow` / `ScilModuleType` / `JSDraw2ModuleType` / `OrgType` / `DojoType` / `DojoxType` |
| `consts.ts` | Re-exports `HelmTypes`, `MonomerTypes`, `PolymerTypes`, `MonomerNumberingTypes` |
| `helm-helper.ts` | `IHelmHelper` interface — `getMolfiles()`, `createHelmInput()`, `createHelmWebEditor()`, `createWebEditorApp()`, `getHoveredAtom()`, `buildMonomersFuncsFromLib()` and `seqHelper` getter. `HelmInputBase` abstract class. `HelmNotSupportedError` (+ `HelmNotSupportedErrorType`). `getHelmHelper()` factory. `getMonomerHandleArgs()` helper. `IHelmInputInitOptions` exported. Augments `ui.input` with `helmAsync()` |
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
| `types.ts` | `MonomerMap` / `MonomerMapValue` / `MolfileWithMap` (a molecule with its per-position monomer atoms and bonds), `UserLibSettings` |
| `monomer-works.ts` | `helmTypeToPolymerType()` converter |
| `assembler.ts` | **Sequence → molecule.** Every conversion runs on Chem's macromolecule assembler (Crux; `IMacromoleculeAssembler` and its types in chem-meta's `macromolecule-assembler.ts`): `getMacromoleculeAssembler()` (session singleton on `window.$macromoleculeAssemblerPromise`, from `Chem:getMacromoleculeAssembler`), `setAssemblerLibrary()` (the effective library once per key), `assembleSeqColumn()` (any notation, `layout: 'linear' \| 'nonlinear'`, the linear one falling back per row; column `molfile(<seq>)`, Molecule/molblock, tags `.sequence-src-col` and `SEQUENCE_SRC_LAYOUT_TAG` `.sequence-src-layout`; per-row errors as warnings), `assembleRowWithMap()` (one row on the main-thread module with its map), `readyAssemblerLibrary()` + `assembleSeqColumnSync()` / `assembleRowWithMapSync()` (the synchronous twins for the legacy synchronous APIs: `IHelmToMolfileConverter`, `helmToAtomicLevelSingle`), `asmMapToMonomerMap()` (the assembler's atom/bond provenance → `MonomerMap`), `MonomerMapProvider` (per-row maps on demand for the hover: LRU 100, a sequential queue, the hovered row first), `getSeqSrcLayout()` (every layout but 'linear' is non-linear) |
| `to-atomic-level.ts` | `_toAtomicLevel()` — the linear conversion (`assembleSeqColumn(..., 'linear')`; the `rdKitModule` param is unused) |
| `seq-to-molfile.ts` | `getMolHighlight()` — builds an `ISubstruct` from per-monomer atom/bond lists for hover highlights |
| `assembly-input.ts` | Input of the macromolecule assembler (`IMacromoleculeAssembler`, types in chem-meta's `macromolecule-assembler.ts`; used by `assembler.ts`): `toAsmSequence()` turns one split sequence of any notation (HELM, BILN, custom via its HELM, FASTA/separator with DNA/RNA expanded to sugar/base/phosphate) into an `AsmSequence` with exact library keys (`MonomerResolver`, through `IMonomerLibBase.getMonomer`, quirks included), inline-SMILES monomers inline, gaps kept; per-row errors for unknown monomers, ambiguity, repeats/ranges, `?`, BLOB and the splitter's `graphInfo.errors` (HELM hwe cannot read). `buildSeqBatch()` packs rows into an `AsmBatch`. `IAsmGraphInfoExt` is the hook a splitter fills (chain ids, biotypes, hbond kind, repeats, errors); `expandAsmRepeats()` |
| `library-json.ts` | The effective monomer library for the assembler: `getEffectiveMonomers()`, `getMonomerLibJson()` (Pistoia JSON array, cached per key), `getMonomerLibKey()` (source + a counter bumped on `onChanged` + a hash of the effective sources and duplicate winners), `getSeqColMonomerLib()` (the column's `MmcrTemps.overriddenLibrary` or the system library), `isLibraryMonomer()` |
| `monomer-hover.ts` | `buildMonomerHoverLink()` — hover from a sequence cell to a molecule cell; the maps come per row from the assembler (`MonomerMapProvider`, the layout of the molecule column's `.sequence-src-layout` tag): the first hover requests the map and redraws the cell when it arrives; with `.sequence-src-highlight-monomers` on, the rows the renderer draws are requested and their grids redrawn. `addMonomerHoverLink()`, `execMonomerHoverLinks()` / `getMonomerHoverLinks()`. `MonomerHoverLinksTemp` column-tag constant |
| `monomer-utils.ts` | `createJsonMonomerLibFromSdf()` (SDF table → HELM JSON library), `sequenceChemSimilarity()` (per-position monomer similarity to a reference), `getMonomerSubstitutionMatrix()` / `calculateMonomerSimilarity()` (monomer similarity matrices from fingerprints) |
| `lib-settings.ts` | `getUserLibSettings()` / `setUserLibSettings()` — user-specific monomer library preferences with chunked storage |
| `utils.ts` | Small shared helpers: `getUnusedColName()`, `getMolColName()`, `hexToPercentRgb()`. `MonomerHoverLink` type |

**Data flow for sequence → molfile conversion:**
1. `_toAtomicLevel()` (linear) and Bio's `sequenceToMolfile()` / `SeqHelper.helmToAtomicLevel()` call `assembleSeqColumn()` in `assembler.ts`.
2. `buildSeqBatch()` (`assembly-input.ts`) turns every row into the parsed sequence; `getMonomerLibJson()` (`library-json.ts`) gives the effective library and its key, registered once per key with `setLibrary`.
3. `assembler.assemble(batch, key, {layout})` builds the V3000 molblocks in Chem's Crux workers (stereo absolute, no OpenChemLib, no RDKit); row errors become warnings.
4. The hover asks for one row's map at a time (`assembleOne(..., {layout, map: true})`, `MonomerMapProvider`).

### `utils/` — Core Utilities

#### Macromolecule subsystem (`utils/macromolecule/`)

| File | Purpose |
|---|---|
| `types.ts` | **Key file.** `ISeqSplitted` (split sequence interface) + `SeqSplittedBase`, `INotationProvider` (notation-specific behavior) + `NotationProviderBase`, `SplitterFunc` type, `IMonomerCanonicalizer`, `CandidateType` / `CandidateSimType` (alphabet detection), `MonomerFreqs`, `SeqColStats`, graph types `ISeqConnection` / `ISeqGraphInfo` |
| `consts.ts` | **Key file.** `NOTATION` enum (FASTA, SEPARATOR, HELM, BILN, CUSTOM), `ALIGNMENT` enum, `ALPHABET` enum (DNA, RNA, PT, UN), `TAGS` enum (column tag names — `units`, `aligned`, `alphabet`, `separator`, `region`, plus also exported as `BioTags`), `Alphabets` namespace (per-alphabet character sets), `GAP_SYMBOL`, `GapOriginals`, `candidateAlphabets`, `monomerRe`, `MONOMER_MOTIF_SPLITTER`, `MONOMER_CANONICALIZER_FUNC_TAG` / `MONOMER_CANONICALIZER_TEMP`, `NOTATION_PROVIDER_CONSTRUCTOR_ROLE`, `positionSeparator` |
| `utils.ts` | **Largest utility file.** Splitter implementations: `StringListSeqSplitted` (separator), `FastaSimpleSeqSplitted`, `HelmSplitted` (abstract, with graph info), `BilnSeqSplitted` (cycles/connections). Splitter factory funcs: `splitterAsFasta`, `splitterAsFastaSimple`, `splitterAsHelm`, `splitterAsBiln`, `getSplitterWithSeparator`, `getSplitter`. Alphabet detection: `detectAlphabet`, `detectHelmAlphabet`, `getAlphabet`, `getAlphabetSimilarity`, `getStatsForCol`, `pickUpPalette` / `pickUpSeqCol` / `getPaletteByType`. Monomer abbreviation (`monomerToShort`). HELM helper: `polymerTypeToHelmType` |
| `seq-handler.ts` | `ISeqHandler` interface — central abstraction for a macromolecule column: notation, alphabet, separator, splitter/joiner, HELM conversion, region extraction, distance function selection. `SeqValueBase` — wraps a single row value with `getOriginal()`, `getCanonical()`, `helm` |
| `scoring.ts` | `ScoringMethod` enum + `calculateIdentityScoring()` / `calculateChemSimilarityScoring()` — sequence scoring vs reference |
| `alignment.ts` | `pairwiseAlignmentWithEmptyPositions()` — Needleman-Wunsch alignment with gap penalties |
| `monomers.ts` | Precomputed Morgan fingerprints for the 20 amino acids and 5 nucleobases: `naturalMonomerFps`, `getMorganFingerprint()`, `mostSimilarNaturalAnalog()` |
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
| `const.ts` | HELM monomer library JSON schema field names (`HELM_REQUIRED_FIELD`, `HELM_OPTIONAL_FIELDS`, `HELM_RGROUP_FIELDS`), SDF-to-JSON mapping, dummy monomer template, RNA sugar/phosphate symbols |
| `splitter.ts` | `splitAlignedSequences()` — splits aligned sequences into per-position columns |
| `composition-table.ts` | `getCompositionTable()` — builds HTML monomer composition bar chart |
| `fasta-handler.ts` | `FastaFileHandler` — parses FASTA files into DataFrames |
| `generator.ts` | Test data generators for synthetic macromolecule columns |
| `sequence-position-scroller.ts` | MSA header / scroll bar engine: `MSAScrollingHeader`, abstract `MSAHeaderTrack`, plus `WebLogoTrack` and `ConservationTrack` track implementations. Renders sequence-logo + conservation tracks above an MSA grid with a scrollable position slider and composition tooltips. The ruler and the tooltips name positions by the column's `.positionLabels` / `.positionNames` (`getPositionName`, passed to tracks by `addTrack`); the header shows its event element only while the grid paints it and adds `header tracks/positions of <col>` readings and `position <name> of <col> header` areas to the grid's status |
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
| `monomer-ui.ts` | Interfaces for monomer editor UI (`INewMonomerForm` + `MonomerFormFields`, `IMonomerManager`) |

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
- `splitterAsHelm` → `HelmSplitted` — HELM parsed by hwe (lenient) with graph/connection info (`ISeqGraphInfoHelm`: per-monomer `polymerTypes`, chain ids, hydrogen bonds, repeats, row `errors`; invalid HELM keeps what hwe read, with `recovered` set)
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
| Assembler input (parsed sequences) | `src/monomer-works/assembly-input.ts` |
| Sequence → molecule (assembler) | `src/monomer-works/assembler.ts` |
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
