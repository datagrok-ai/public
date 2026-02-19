# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

**Bio** (`@datagrok/bio`) is the main Datagrok plugin for **bioinformatics** — automatic detection, rendering, editing, analysis, and conversion of macromolecule sequences (peptides, DNA, RNA, HELM, BILN). It provides semantic type detection, custom cell renderers, viewers (WebLogo, VD Regions), analysis tools (sequence space, activity cliffs, MSA, similarity/diversity search), monomer library management, substructure filtering, and atomic-level conversion.

Category: **Bioinformatics**. Top menu: `Bio | ...`.

## Build Commands

```bash
npm install
npm run build              # grok api && grok check --soft && webpack
npm run test               # grok test
npm run lint               # eslint src --ext .ts
npm run lint-fix
npm run build-all          # Builds chem-meta → js-api → utils → bio-lib → this package
npm run link-all           # Links local datagrok-api and @datagrok-libraries/*
```

## Key Dependencies

- `@datagrok-libraries/bio` — **core shared library**: macromolecule types, HELM types, splitters, palettes, MonomerWorks (seq→molfile), cell renderer engines, monomer hover, seq-helper/seq-handler interfaces
- `@datagrok-libraries/ml` — distance matrices, dimensionality reduction (UMAP/tSNE), activity cliffs, KNN, MCL clustering, macromolecule distance functions
- `@datagrok-libraries/math` — DBSCAN worker, WebGPU utilities
- `@datagrok-libraries/utils` — BitArray, type declarations, SVG utilities
- `@datagrok-libraries/chem-meta` — RDKit API types (`RDModule`)
- `@datagrok-libraries/tutorials` — tutorial framework
- `@biowasm/aioli` — WebAssembly runtime for kalign (MSA)
- `ajv` / `ajv-errors` — JSON schema validation for monomer libraries
- `openchemlib` — molfile conversion, chirality engine

### Relationship to `@datagrok-libraries/bio`

This **package** implements interfaces defined in the **library**. Many capabilities follow a service pattern: the library defines an interface + `getXxxHelper()` factory, and this package registers the concrete implementation via `DG.Func.find()`.

| Interface (in bio lib) | Implementing class (this package) |
|---|---|
| `ISeqHelper` | `SeqHelper` (in `src/utils/seq-helper/seq-helper.ts`) |
| `IMonomerLibHelper` | `MonomerLibManager` (in `src/utils/monomer-lib/lib-manager.ts`) |

Other packages depend on these implementations: **Helm**, **Peptides**, **BiostructureViewer**, **Dendrogram**, **HitTriage**, etc.

## Architecture

### Entry Point — `src/package.ts`

`PackageFunctions` class registers all platform-visible functions via `@grok.decorators`:

**Initialization:**
- `initBio` — `@init`: loads RDKit module, package properties, monomer library (via `MonomerLibManager`), monomer sets, creates `SeqHelper` singleton

**Registered Functions (Top Menu `Bio | ...`):**

| Menu Path | Function | Purpose |
|---|---|---|
| `Bio \| Analyze \| Activity Cliffs...` | `activityCliffs` | Detects sequence pairs with similar structure but significant activity difference |
| `Bio \| Analyze \| Sequence Space...` | `sequenceSpaceTopMenu` | UMAP/tSNE 2D projection of sequences by pairwise distance |
| `Bio \| Analyze \| MSA...` | `multipleSequenceAlignmentDialog` | Multiple sequence alignment via kalign (WASM) or PepSeA (Docker) |
| `Bio \| Analyze \| Composition` | `compositionAnalysis` | Docks a WebLogo viewer for sequence composition |
| `Bio \| Transform \| Convert Sequence Notation...` | `convertDialog` | FASTA ↔ SEPARATOR ↔ HELM ↔ BILN conversion |
| `Bio \| Transform \| To Atomic Level...` | `toAtomicLevel` | Converts sequences to V3000 molfiles |
| `Bio \| Transform \| Split to Monomers...` | `splitToMonomersTopMenu` | Splits aligned sequences into per-position columns |
| `Bio \| Transform \| Molecules to HELM...` | `moleculesToHelmTopMenu` | Converts peptide molecules to HELM via Python script |
| `Bio \| Calculate \| Get Region...` | `getRegionTopMenu` | Extracts sub-region from macromolecule column |
| `Bio \| Calculate \| Identity...` | `sequenceIdentityScoring` | Fraction of matching monomers vs reference |
| `Bio \| Calculate \| Similarity...` | `sequenceSimilarityScoring` | Sum of monomer fingerprint similarities vs reference |
| `Bio \| Search \| Similarity Search` | `similaritySearchTopMenu` | K-nearest neighbor sequence search |
| `Bio \| Search \| Diversity Search` | `diversitySearchTopMenu` | Maximally diverse subset selection |
| `Bio \| Search \| Subsequence Search...` | `SubsequenceSearchTopMenu` | Substructure filter (regex or RDKit) |
| `Bio \| Manage \| Monomer Libraries` | `manageLibrariesView` | Full monomer library management UI |
| `Bio \| Manage \| Monomers` | `manageMonomersView` | Individual monomer CRUD editor |
| `Bio \| Manage \| Match with Monomer Library...` | `matchWithMonomerLibrary` | Matches molecules to library monomers |

**Registered Viewers:**
- `WebLogo` → `WebLogoViewer` — sequence logo with entropy/full-height modes
- `VdRegions` → `VdRegionsViewer` — antibody V-domain region viewer (FR1-4, CDR1-3)
- `Sequence Similarity Search` → `SequenceSimilarityViewer`
- `Sequence Diversity Search` → `SequenceDiversityViewer`

**Cell Renderers:**
- `fastaSequenceCellRenderer`, `separatorSequenceCellRenderer`, `bilnSequenceCellRenderer`, `customSequenceCellRenderer` → `MacromoleculeSequenceCellRenderer` (for `Macromolecule` columns, units: fasta/separator/biln/custom)
- `MacromoleculeDifferenceCellRenderer` — renders sequence diffs
- `monomerCellRenderer` → `MonomerCellRenderer` — renders individual monomer cells

**Property Panels:**
- `macroMolColumnPropertyPanel` — renderer settings (font size, color coding, reference seq, multiline)
- `compositionAnalysisWidget` — monomer composition table for a cell
- `toAtomicLevelPanel` — 2D molecule widget from sequence
- `sequence3dStructureWidget` — 3D NGL viewer widget from sequence
- `libraryPanel` — link to monomer library manager
- `getRegionPanel` — extract region UI
- `sequenceTooltip` — WebLogo column tooltip

**File Handlers:**
- `importFasta` — opens `.fasta`, `.fna`, `.fa`, etc. files
- `importBam` — stub for BAM files

**Other:**
- `bioSubstructureFilter` — macromolecule substructure filter
- `saveAsFasta` — file exporter (As FASTA...)
- `Encode Sequences` / `Helm Fingerprints` — preprocessing functions for dimensionality reduction
- `getSeqHelper`, `getMonomerLibHelper`, `getBioLib`, `getSeqHandler` — service accessors for other packages
- `refineNotationProviderForBiln` — BILN notation refiner
- `addCopyMenu` — context menu "Copy as..." for different notations

**Exports**: `_package` (instance of `BioPackage`), auto-generated function wrappers

### Package Singleton — `src/package-types.ts`

`BioPackage` (extends `DG.Package`) — the package singleton managing:
- `seqHelper: ISeqHelper` — sequence operations facade
- `monomerLib: IMonomerLib` — current merged monomer library
- `monomerSets: IMonomerSet` — monomer set preferences
- `rdKitModule: RDModule` — RDKit WebAssembly module
- `properties: BioPackageProperties` — observable package settings (FontSize, MaxMonomerLength, TooltipWebLogo, DefaultSeparator)
- `completeInit()` — wires up all singletons after loading

### Semantic Type Detector — `detectors.js`

`BioPackageDetectors.detectMacromolecule(col)` — **synchronous** detector (714 lines) that classifies string columns as `Macromolecule`:
1. Rejects non-string columns, URL-like data, and columns resembling SMILES/SMARTS/numbers
2. Tests for **HELM** notation first (starts with `PEPTIDE1{`, `RNA1{`, etc.)
3. Detects **separators** by analyzing character frequency distributions
4. Detects **BILN** notation (separator `-` with connection patterns like `(1,2)`)
5. Classifies **alphabet** (DNA, RNA, PT=peptide, UN=unknown) via cosine similarity to known monomer sets
6. Determines **alignment** (SEQ vs SEQ.MSA) by checking if all sequences have the same length
7. Sets column tags: `units` (fasta/separator/helm), `aligned`, `alphabet`, `separator`, `.alphabetIsMultichar`
8. Column name heuristics (`peptide`, `sequence`, `oligo`, etc.) lower detection thresholds

Key methods: `detectSeparator()`, `detectAlphabet()`, `getAlphabetSimilarity()`, `getStats()`, `checkBadMultichar()`

## Source Structure

### Viewers (`src/viewers/`)

| File | Class | Purpose |
|---|---|---|
| `web-logo-viewer.ts` | `WebLogoViewer` | Canvas-based sequence logo viewer (~1430 lines). Renders monomer frequency/composition per position. Modes: Entropy (information content height) and Full (equal height). Supports aggregation by column, horizontal range slider, tooltips, click-to-filter. Properties: sequenceColumnName, positionWidth, startPosition, endPosition, mode, etc. |
| `vd-regions-viewer.ts` | `VdRegionsViewer` | Antibody V-domain regions viewer (~530 lines). Displays FR1-4, CDR1-3 regions for Heavy/Light chains in a grid, each cell containing a child WebLogo. Supports IMGT-style nomenclature, fit-to-width, chain/region filtering. |
| `utils.ts` | — | Single helper: `isReallyChanged(a, b, eps)` |

### Analysis (`src/analysis/`)

| File | Purpose |
|---|---|
| `sequence-space.ts` | `getEncodedSeqSpaceCol()` — encodes sequences to UTF chars for distance computation, builds scoring matrix/gap penalty options for Needleman-Wunsch, Hamming, Levenshtein, Monomer Chemical Distance |
| `sequence-activity-cliffs.ts` | Activity cliff utilities: `getSequenceMoleculeDistances()`, `getSequenceMoleculeSimilarities()`, `createTooltipElement()`, `createPropPanelElement()`, `drawDifference()` (renders visual sequence diff), `createLinesGrid()` |
| `sequence-similarity-viewer.ts` | `SequenceSimilarityViewer` — KNN-based viewer finding K most similar sequences to current row. Caches sparse KNN matrix. |
| `sequence-diversity-viewer.ts` | `SequenceDiversityViewer` — selects maximally diverse sequence subset via full distance matrix. Random sampling for large datasets (>10k rows). |
| `sequence-search-base-viewer.ts` | `SequenceSearchBaseViewer` — abstract base for similarity/diversity viewers. Properties: distanceMetric, fingerprintType, gapOpen, gapExtend. |

### Widgets (`src/widgets/`)

| File | Purpose |
|---|---|
| `bio-substructure-filter.ts` | `BioSubstructureFilter` — collaborative filter for macromolecule columns. Routes to FASTA text filter, separator filter, or HELM editor filter based on notation. |
| `bio-substructure-filter-helm.ts` | `HelmSubstructureFilterEditor` — HELM-specific filter with embedded HELM web editor for query drawing |
| `composition-analysis-widget.ts` | `getCompositionAnalysisWidget()` — monomer composition table with color-coded counts for a single cell value |
| `representations.ts` | `getMacromoleculeColumnPropertyPanel()` — UI for renderer settings: font size, max monomer length, gap length, color coding scheme, reference sequence, multiline mode |
| `to-atomic-level-widget.ts` | `toAtomicLevelSingle()` (seq→molfile), `toAtomicLevelWidget()` (2D molecule drawing), `molecular3DStructureWidget()` (3D NGL viewer) |
| `sequence-scrolling-widget.ts` | `handleSequenceHeaderRendering()` — MSA column header with WebLogo + conservation tracks, viewport-aware lazy caching (50-position chunks), click to dock position statistics viewer |
| `package-settings-editor-widget.ts` | `PackageSettingsEditorWidget` — Bio package global settings form |

### Utilities (`src/utils/`)

#### Core Utilities

| File | Purpose |
|---|---|
| `cell-renderer.ts` | `MacromoleculeSequenceCellRenderer` — main cell renderer for macromolecule columns. Delegates to `MonomerPlacer` back-end from bio library. Also `MacromoleculeDifferenceCellRenderer` / `MacromoleculeDoubleCellRenderer` for showing diffs between two aligned sequences. |
| `monomer-cell-renderer.ts` | `MonomerCellRenderer` — renders individual monomer cells with color from library. Extends `MonomerCellRendererBase`. |
| `monomer-cell-renderer-base.ts` | `MonomerCellRendererBase` — abstract base that async-loads monomer library, invalidates grid on lib change |
| `convert.ts` | `convert()` — dialog for notation conversion (FASTA ↔ SEPARATOR ↔ HELM ↔ BILN). `convertDo()` performs actual conversion via `ISeqHandler`. |
| `get-region.ts` | `getRegionDo()` — extracts positional sub-region from macromolecule column via `ISeqHandler.getRegion()` |
| `split-to-monomers.ts` | `splitToMonomersUI()` — splits aligned sequences into per-position `Monomer` columns using `joinDataFrames` |
| `sequence-to-mol.ts` | `sequenceToMolfile()` — converts macromolecule column to atomic-level molfile column (linear via `_toAtomicLevel`, nonlinear via HELM converter) |
| `calculate-scores.ts` | `calculateScoresWithEmptyValues()` — wraps `calculateIdentityScoring`/`calculateChemSimilarityScoring` from bio lib, handles empty values |
| `save-as-fasta.ts` | `saveAsFastaUI()` — dialog + download for FASTA export. `buildFasta()` builds the string. |
| `biln.ts` | `BilnNotationProvider` — notation provider for BILN sequences with splitter, HELM converter, cell renderer back-end |
| `context-menu.ts` | `addCopyMenuUI()` — adds "Copy as FASTA/SEPARATOR/HELM/BILN" to cell context menu |
| `check-input-column.ts` | `checkInputColumnUI()` / `checkInputColumn()` — validates column has `Macromolecule` semtype |
| `ui-utils.ts` | `getMacromoleculeColumns()`, `safeReplace()`, `setGridColWidth()` |
| `types.ts` | Shared types: `DfPair`, `ActivityCliffsData`, `MsaOptions` (kalign/pepsea params), `AARDict` |
| `constants.ts` | Constants: `SEM_TYPES`, `PEPSEA_VERSION`, `DEFAULT_MSA_PARAMETERS`, amino acid groupings |
| `agg.ts` | `getAggregatedValue()` — generic aggregation helper using `DG.DataFrame.aggregate()` |

#### Multiple Sequence Alignment

| File | Purpose |
|---|---|
| `multiple-sequence-alignment.ts` | `multipleSequenceAlignment()` — core MSA via **kalign** (WebAssembly/Aioli). Supports per-cluster alignment, gap penalties, selected-rows-only mode. |
| `multiple-sequence-alignment-ui.ts` | `multipleSequenceAlignmentUI()` — MSA dialog with column selection, alignment method (kalign vs PepSeA), gap penalties |
| `pepsea.ts` | `pepseaAlignSequences()` — MSA for HELM peptides via **PepSeA Docker container** (mafft/linsi/ginsi methods) |

#### Seq Helper — `src/utils/seq-helper/`

| File | Purpose |
|---|---|
| `seq-helper.ts` | `SeqHelper` (implements `ISeqHelper`) — central facade for sequence operations: `getSeqHandler()`, `getSeqMonomers()`, `helmToMolfileV3K()`, `helmToSmiles()`, column tag management |
| `seq-handler.ts` | `SeqHandler` (implements `ISeqHandler`) — per-column handler for splitting, notation conversion (FASTA/SEP/HELM/BILN), region extraction, stats, alphabet detection. Contains joiner functions. |
| `index.ts` | Barrel re-export of `SeqHelper` |

#### Monomer Library — `src/utils/monomer-lib/`

| File | Purpose |
|---|---|
| `lib-manager.ts` | `MonomerLibManager` (implements `IMonomerLibHelper`) — **singleton** managing the entire monomer library system. Discovers `IMonomerLibProvider` instances (file-based + custom), loads/merges libraries per user settings, exposes `onChanged` observables. Singleton via `getInstance()`. |
| `monomer-lib.ts` | `MonomerLib` (extends `MonomerLibBase`, implements `IMonomerLib`) — merged library from multiple sources. Handles duplicate tracking, user duplicate preferences, summary stats. |
| `monomer-lib-base.ts` | `MonomerLibBase` — base implementation: monomer lookup, symbol listing, missing monomer creation, R-group extraction from SMILES, tooltip rendering, color computation with contrast logic |
| `monomer-colors.ts` | Static color mappings for natural monomers (nucleotide chromatogram palette, amino acid GrokGroups palette) |
| `consts.ts` | Test constants: `LIB_SETTINGS_FOR_TESTS`, `LIB_MONOMER_COUNTS` |
| `smiles2Monomer.ts` | `smiles2Monomer()` — converts inline SMILES (CX-SMILES with R-group labels) to `Monomer` objects with auto-derived R-groups |
| `web-editor-monomer-dummy.ts` | Placeholder `WebEditorMonomer` implementations (inline SMILES, gap, ambiguous, missing, broken) |
| `web-editor-monomer-of-library.ts` | `getWebEditorMonomerOfLib()` — constructs HELM web editor monomer from real `Monomer` object in library |

##### Library File Manager — `src/utils/monomer-lib/library-file-manager/`

| File | Purpose |
|---|---|
| `ui.ts` | UI for library management: `showManageLibrariesDialog()` (checkbox per library with edit/delete), `showManageLibrariesView()` (full view with duplicate manager), `getMonomerLibraryManagerLink()` (panel widget) |
| `MonomerLibFromFilesProvider` | Reads/writes monomer library JSON files from `System:AppData/Bio/` file shares. CRUD operations, HELM JSON schema validation. |
| `validator` | AJV-based JSON schema validation for monomer library files |

##### Monomer Manager — `src/utils/monomer-lib/monomer-manager/`

| File | Purpose |
|---|---|
| `monomer-manager.ts` | `MonomerManager` — full CRUD UI for monomers within a library. Grid view with context menus, create/edit dialog with SMILES↔molfile standardization. `matchMoleculesWithMonomers()` — matches molecules against library by canonical SMILES. `standardizeMonomerLibrary()` — normalizes a library JSON. |
| `duplicate-manager.ts` | `DuplicateManager` — UI for resolving duplicate monomer conflicts across libraries. Monomer cards grouped by symbol, user picks preferred source. |
| `default-r-groups.ts` | Default R-group definitions (R1-R6 with cap group SMILES) |

#### HELM to Molfile — `src/utils/helm-to-molfile/`

| File | Purpose |
|---|---|
| `utils.ts` | `getMolColumnFromHelm()` (HELM → beautified V3K molfile column), `getSmilesColumnFromHelm()` (HELM → SMILES column) |
| `converter/` | Full pipeline (~20 files) for HELM → V3K molfile conversion. `HelmToMolfileConverter` class handles HELM decomposition → monomer assembly → composite V3K with monomer-position maps and chirality (via OpenChemLib). |

#### Substructure Search — `src/substructure-search/`

| File | Purpose |
|---|---|
| `substructure-search.ts` | `SubstructureSearchDialog` — UI for subsequence/substructure search. `linearSearch()` — regex match on FASTA/separator. `helmSearch()` — converts to monomer-level pseudo-molfiles, delegates to `Chem:searchSubstructure`. `invalidateMols()` — caches monomeric molfile column. |

#### Calculations — `src/calculations/`

| File | Purpose |
|---|---|
| `monomerLevelMols.ts` | `getMonomerLevelMols()` — creates pseudo-molfile column where each monomer becomes a unique pseudo-atom (mass-labeled `At`), for use in substructure search and fingerprinting |

### Other Source Files

| File | Purpose |
|---|---|
| `seq_align.ts` | `SequenceAlignment` — Needleman-Wunsch / Smith-Waterman pairwise alignment with BLOSUM matrices (45/50/62/80/90) |
| `function-edtiors/split-to-monomers-editor.ts` | `SplitToMonomersFunctionEditor` — custom dialog for Split to Monomers function |
| `apps/web-logo-app.ts` | `WebLogoApp` — test app that opens a table with a docked WebLogo viewer |
| `apps/get-region-app.ts` | `GetRegionApp` — test/demo app for region extraction with named positions |

### Demos (`src/demo/`)

| File | Demo |
|---|---|
| `bio01-similarity-diversity.ts` | Similarity & diversity search side-by-side |
| `bio01a-hierarchical-clustering-and-sequence-space.ts` | UMAP sequence space + hierarchical clustering |
| `bio01b-hierarchical-clustering-and-activity-cliffs.ts` | Activity cliffs + UMAP + clustering |
| `bio03-atomic-level.ts` | HELM → atomic level conversion with monomer hover highlighting |
| `bio05-helm-msa-sequence-space.ts` | HELM MSA via PepSeA + sequence space scatter plot |
| `utils.ts` | Shared: `demoSequenceSpace()` (UMAP/tSNE + scatter), constants |

## Tests (`src/tests/`)

Test entry point: `src/package-test.ts` — imports all test files, exports `test()` and `initAutoTests()`.

| File | What it tests |
|---|---|
| `_first-tests.ts` | Package import sanity |
| `detectors-tests.ts` | Macromolecule detection on various data |
| `detectors-weak-and-likely-tests.ts` | Edge cases for detector heuristics |
| `detectors-benchmark-tests.ts` | Detector performance benchmarks |
| `splitters-test.ts` | Sequence splitting for different notations |
| `seq-handler-tests.ts` | SeqHandler operations |
| `seq-handler-splitted-tests.ts` | Splitted sequence handling |
| `seq-handler-get-region-tests.ts` | Region extraction |
| `seq-handler-get-helm-tests.ts` | HELM conversion from other notations |
| `renderers-test.ts` | Cell renderer behavior |
| `renderers-monomer-placer-tests.ts` | MonomerPlacer rendering engine |
| `converters-test.ts` | Notation conversion |
| `bio-tests.ts` | General bio functionality |
| `helm-tests.ts` | HELM-specific operations |
| `msa-tests.ts` | Multiple sequence alignment |
| `pepsea-tests.ts` | PepSeA Docker MSA |
| `to-atomic-level-tests.ts` | Sequence → molfile conversion |
| `to-atomic-level-ui-tests.ts` | Atomic level UI widgets |
| `activity-cliffs-tests.ts` | Activity cliff detection |
| `sequence-space-test.ts` | Sequence space UMAP/tSNE |
| `similarity-diversity-tests.ts` | Similarity/diversity viewers |
| `substructure-filters-tests.ts` | Substructure filter |
| `monomer-libraries-tests.ts` | Monomer library loading/merging |
| `lib-tests.ts` | Library management |
| `scoring.ts` | Identity/similarity scoring |
| `WebLogo-positions-test.ts` | WebLogo position handling |
| `WebLogo-project-tests.ts` | WebLogo save/load projects |
| `WebLogo-layout-tests.ts` | WebLogo layout/resize |
| `viewers.ts` | Viewer integration tests |
| `fasta-handler-test.ts` | FASTA file import |
| `fasta-export-tests.ts` | FASTA file export |
| `Palettes-test.ts` | Color palette tests |
| `biln-tests.ts` | BILN notation |
| `mm-distance-tests.ts` | Macromolecule distance functions |
| `checkInputColumn-tests.ts` | Input validation |

## Initialization Flow

1. Platform loads `detectors.js` → `detectMacromolecule()` classifies string columns as `Macromolecule` (synchronous, ~714 lines)
2. `initBio()` called on package init:
   a. Loads RDKit WebAssembly module
   b. Loads package properties (`BioPackageProperties`)
   c. Creates `MonomerLibManager` singleton → discovers library providers → loads/merges monomer libraries from `System:AppData/Bio/` files
   d. Loads monomer sets (async, non-blocking)
   e. Creates `SeqHelper` singleton
   f. Calls `_package.completeInit()` to wire everything together
   g. Calls `handleSequenceHeaderRendering()` to set up MSA column header tracks

## Package Properties (Settings)

| Property | Type | Default | Description |
|---|---|---|---|
| `FontSize` | int | 12 | Font size for monomer symbols in renderer |
| `MaxMonomerLength` | string | "4" | Max displayed monomer symbol length ("long" = unlimited) |
| `TooltipWebLogo` | bool | true | Show WebLogo in column header tooltip |
| `DefaultSeparator` | string | "." | Default separator for notation conversion |

## Data Files

- `files/monomer-libraries/` — HELM monomer library JSON files (HELMCoreLibrary, polytool-lib, sample-lib)
- `files/monomer-sets/` — Monomer set definitions (PEPTIDE, RNA)
- `files/samples/` — Sample data files for all notation types (FASTA, SEPARATOR, HELM, BILN, MSA)
- `files/tests/` — Test data files
- `files/schemas/` — JSON schemas for monomer sets
- `files/icons/` — Viewer icons

## Quick Lookups

| Looking for... | Check first |
|---|---|
| Function/panel/viewer registration | `src/package.ts` (`PackageFunctions` class) |
| Package singleton, init flow | `src/package-types.ts` (`BioPackage`) + `initBioInt()` in `src/package.ts` |
| Semantic type detection | `detectors.js` (`detectMacromolecule`) |
| Sequence operations facade | `src/utils/seq-helper/seq-helper.ts` (`SeqHelper`) |
| Per-column sequence handling | `src/utils/seq-helper/seq-handler.ts` (`SeqHandler`) |
| Monomer library management | `src/utils/monomer-lib/lib-manager.ts` (`MonomerLibManager`) |
| Monomer library data model | `src/utils/monomer-lib/monomer-lib.ts` (`MonomerLib`) |
| Monomer library UI | `src/utils/monomer-lib/library-file-manager/ui.ts` |
| Monomer CRUD editor | `src/utils/monomer-lib/monomer-manager/monomer-manager.ts` |
| Cell renderer (sequence columns) | `src/utils/cell-renderer.ts` |
| Cell renderer (monomer columns) | `src/utils/monomer-cell-renderer.ts` |
| WebLogo viewer | `src/viewers/web-logo-viewer.ts` |
| VD Regions viewer | `src/viewers/vd-regions-viewer.ts` |
| Substructure search/filter | `src/substructure-search/substructure-search.ts` + `src/widgets/bio-substructure-filter.ts` |
| Activity cliffs | `src/analysis/sequence-activity-cliffs.ts` |
| Sequence space (UMAP/tSNE) | `src/analysis/sequence-space.ts` |
| Similarity/diversity viewers | `src/analysis/sequence-similarity-viewer.ts` / `sequence-diversity-viewer.ts` |
| MSA (kalign) | `src/utils/multiple-sequence-alignment.ts` |
| MSA (PepSeA Docker) | `src/utils/pepsea.ts` |
| Notation conversion | `src/utils/convert.ts` |
| Seq → molfile conversion | `src/utils/sequence-to-mol.ts` |
| HELM → molfile pipeline | `src/utils/helm-to-molfile/converter/` |
| FASTA import/export | `src/utils/save-as-fasta.ts` + bio lib `FastaFileHandler` |
| Pairwise alignment | `src/seq_align.ts` (`SequenceAlignment`) |
| BILN notation support | `src/utils/biln.ts` |
| MSA header with WebLogo | `src/widgets/sequence-scrolling-widget.ts` |
| Atomic level widgets | `src/widgets/to-atomic-level-widget.ts` |
| Demo scripts | `src/demo/` |
| Test entry point | `src/package-test.ts` |
| Auto-generated wrappers | `src/package.g.ts` / `src/package-api.ts` |
| Python scripts | `scripts/` (mol-to-helm.py, sequence_generator.py, embed.py) |
