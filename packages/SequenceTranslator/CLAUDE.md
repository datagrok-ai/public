# SequenceTranslator Package - CLAUDE.md

## Overview

**Package:** `@datagrok/sequence-translator` (v1.10.16)
**Category:** Bioinformatics
**Author:** Davit Rizhinashvili
**Description:** Translates oligonucleotide sequences between different representations (Axolabs, BioSpring, GCRS, Mermade, HELM, nucleotides, etc.). Provides pattern design, structure visualization, advanced polymer conversion (PolyTool), and a dedicated **OligoNucleotide** semantic type / cell renderer for siRNA / ASO duplex visualization on top of HELM.

## Dependencies

| Dependency | Purpose |
|------------|---------|
| `datagrok-api` | Core Datagrok platform API |
| `@datagrok-libraries/bio` | HELM helper, sequence handler, macromolecule utils, central monomer library |
| `@datagrok-libraries/chem-meta` | RDKit JS API, molfile parsing |
| `@datagrok-libraries/utils` | Error handling, common utilities |
| `@datagrok-libraries/tutorials` | Tutorial creation helpers |
| `openchemlib` | Chemical structure rendering |
| `lodash` | Utility functions |
| `save-svg-as-png` | SVG-to-PNG export in pattern app |
| `typeahead-standalone` | Autocomplete for monomer selection |
| `object-hash` | SHA1 hashing for pattern identity |

**Dev-only:** `@datagrok/bio`, `@datagrok/chem`, `@datagrok/helm` (for testing), `@datagrok-libraries/helm-web-editor`, `@datagrok-libraries/js-draw-lite`

## Package Property

- **MonomersPath** (string): Path to additional monomer libraries. Default: `System:AppData/SequenceTranslator/monomers`. Fallback: `System:AppData/SequenceTranslator/monomers-sample`.

> The OligoNucleotide subsystem (under `src/oligo-renderer/`) deliberately uses the **central Bio monomer library** (`HELMCoreLibrary.json` + the new `oligo-conjugates.json` shipped in `packages/Bio/files/monomer-libraries/`), not the ST-internal one. ST's monomer-lib stays scoped to translator/pattern/structure apps.

---

## Architecture

The package is organized into **four main applications**, a **PolyTool subsystem**, the **OligoNucleotide renderer subsystem**, and shared infrastructure:

```
src/
├── package.ts                    # Main entry point, exports all Datagrok functions
├── package.g.ts                  # AUTO-GENERATED - do not edit
├── package-api.ts                # AUTO-GENERATED - typed function wrappers
├── package-test.ts               # Test orchestration
├── consts.ts                     # Package-level constants (PolyToolTags, PolyToolDataRole)
├── types.ts                      # ITranslationHelper interface
│
├── apps/
│   ├── common/                   # Shared model, data loading, validation, UI base
│   ├── translator/               # Format conversion app
│   ├── pattern/                  # Pattern design app
│   └── structure/                # Molecular structure app
│
├── oligo-renderer/               # OligoNucleotide semType, cell renderer, panels, converters
│   ├── types.ts                  # Constants, parsed-model types, modification dictionary, DuplexAlignment
│   ├── helm-parser.ts            # Lightweight RNA/DNA HELM duplex parser + canonicalizer + pair/strandtype parsing
│   ├── alignment.ts              # Sense/antisense register: explicit HELM pairs or auto-align by complementarity
│   ├── canvas-renderer.ts        # Pure-fn canvas drawing of duplex (shift-aware overhangs), hit-testing
│   ├── cell-renderer.ts          # OligoNucleotideCellRenderer (extends DG.GridCellRenderer)
│   ├── tooltip.ts                # Hover tooltip with cached RDKit structures
│   ├── legend-panel.ts           # Per-cell summary + filtered color legend
│   ├── structures-panel.ts       # Sense/antisense full molecular structures (Bio:toAtomicLevelPanel)
│   └── converters.ts             # convertHelmColumnToOligo / combineSenseAntisenseToOligo / tagAsOligoNucleotide
│
├── polytool/                     # Advanced polymer conversion subsystem
│   ├── conversion/               # Core algorithms (chain, rules, reactions)
│   └── *.ts                      # Dialogs, enumeration, handlers
│
├── plugins/                      # External plugin integration (MerMade)
├── utils/                        # Cyclized/dimerized notation providers, error handling
├── demo/                         # Demo entry points
└── tests/                        # Test suite

prototypes/                       # Off-tree dev tooling for the OligoNucleotide subsystem
detectors.js                      # Autostart bootstrap + context-menu wiring (no semType detector)
```

---

## Entry Point: package.ts

Central `PackageFunctions` class exposes all Datagrok-registered functions. Key instance:

```typescript
export const _package: OligoToolkitPackage = new OligoToolkitPackage({debug: true});
```

### Initialization Flow

1. `init()` (decorator `@init`) → gets HELM helper from bio library, calls `completeInit()`
2. `completeInit()` → sets up monomer library wrapper
3. `initLibData()` → lazy-loads JSON data files (called on app open)

### Registered Functions

The table below lists every function/app/panel registered by `PackageFunctions`. **Method name** is the TypeScript identifier; **registered name** is what the platform sees (defaults to method name when no `name:` decorator option is given).

#### Apps

| Method | Registered name | browsePath | Description |
|---|---|---|---|
| `oligoToolkitApp` | Oligo Toolkit | `Peptides \| Oligo Toolkit` | Combined tabbed app (translator + pattern + structure + plugins) |
| `oligoTranslatorApp` | Oligo Translator | `Peptides \| Oligo Toolkit` | Standalone translator |
| `oligoPatternApp` | Oligo Pattern | `Peptides \| Oligo Toolkit` | Standalone pattern designer |
| `oligoStructureApp` | Oligo Structure | `Peptides \| Oligo Toolkit` | Standalone structure viewer |
| `ptEnumeratorHelmApp` | HELM Enumerator | `Peptides \| PolyTool` | HELM enumerator app |
| `ptEnumeratorChemApp` | Chem Enumerator | `Chem \| PolyTool` | Chem enumerator app |

#### Top-Menu Functions

| Method | Registered name | top-menu | Description |
|---|---|---|---|
| `polyToolConvertTopMenu` | polyToolConvert | `Bio \| PolyTool \| Convert...` | Open convert dialog |
| `polyToolEnumerateHelmTopMenu` | polyToolEnumerateHelm | `Bio \| PolyTool \| Enumerate HELM...` | HELM enumerate dialog |
| `polyToolEnumerateChemTopMenu` | polyToolEnumerateChem | `Bio \| PolyTool \| Enumerate Chem...` | Chem enumerate dialog |
| `chemEnumerateReactionsTopMenu` | chemEnumerateReactions | `Chem \| Transform \| Reactions \| Enumerate...` | Alternate menu for chem enumeration |
| `getPolyToolCombineDialog` | Combine Sequences | `Bio \| PolyTool \| Combine Sequences...` | Cartesian-product combine |

#### Cell Renderer

| Method | cellType / columnTags | Description |
|---|---|---|
| `oligoNucleotideCellRenderer` | `cellType: OligoNucleotide`, `columnTags: quality=OligoNucleotide` | Duplex view for HELM-backed OligoNucleotide columns |

#### Panels (`tags: panel, widgets`)

| Method | Registered name | Input semType | Description |
|---|---|---|---|
| `oligoNucleotidePanel` | Oligo-Nucleotide | OligoNucleotide | Per-cell summary, conjugate counts, color legend filtered to the cell |
| `oligoNucleotideStructuresPanel` | Oligo Structures | OligoNucleotide | Lazy accordion: sense / antisense full structures via `Bio:toAtomicLevelPanel` |

#### Public API & Helpers

| Method | Registered name | Notes |
|---|---|---|
| `getTranslationHelper` | getTranslationHelper | Returns `ITranslationHelper` (calls `_package.initLibData()` first) |
| `getCodeToWeightsMap` | getCodeToWeightsMap | Code → MW lookup map |
| `validateSequence` | validateSequence | Boolean validation against detected format |
| `getMolfileFromGcrsSequence` | **`validateSequence`** ⚠️ | GCRS → V3000 molfile. **Note:** the `name:` decorator currently mis-registers this as `validateSequence`, colliding with the previous entry. Method name still works for in-package callers. |
| `linkStrands` | linkStrands | Multi-strand assembly (V3000) |
| `translateOligonucleotideSequence` | translateOligonucleotideSequence | Format-to-format conversion |
| `polyToolConvert2` | polyToolConvert2 | Bulk conversion (uses editor `getPolyToolConvertEditor`) |
| `getPolyToolConvertEditor` | getPolyToolConvertEditor | Function-call editor for `polyToolConvert2` (`@editor` decorator) |
| `polyToolColumnChoice` | polyToolColumnChoice | Mark macromolecule column + redetect semantic types |
| `createMonomerLibraryForPolyTool` | createMonomerLibraryForPolyTool | CSV → JSON monomer library |
| `enumerateSingleHelmSequence` | Enumerate Single HELM Sequence | Programmatic single-HELM enumeration → DataFrame |
| `enumerateSingleHelmSequenceWithNaturalAAs` | Enumerate Single HELM Sequence with natural amino acids | All-positions natural-AA enumeration |

#### Context-Menu Wrappers (invoked from `detectors.js`)

| Method | Registered name | Trigger | Description |
|---|---|---|---|
| `getPtHelmEnumeratorDialog` | Polytool Helm Enumerator dialog | Right-click Macromolecule cell | HELM enumerator dialog seeded by clicked cell |
| `getPtChemEnumeratorDialog` | Polytool Chem Enumerator dialog | Right-click Molecule cell | Chem enumerator dialog seeded by clicked cell |
| `getPtOligoEnumeratorDialog` | Polytool Oligo Enumerator dialog | Right-click OligoNucleotide cell | Wraps oligo HELM in a temp Macromolecule cell, runs HELM enumerator with `outputAsOligo=true` so the result column is tagged OligoNucleotide |
| `convertHelmToOligoNucleotide` | convertHelmToOligoNucleotide | Submenu `Oligo` on HELM cell | Clones a HELM column and tags it as OligoNucleotide |
| `combineSenseAntisenseToOligoNucleotide` | combineSenseAntisenseToOligoNucleotide | Submenu `Oligo` on Macromolecule cell | Opens `.prepare(...).edit()` to pick antisense; produces combined Oligo column |

#### Notation Provider Glue

| Method | Registered name | Notes |
|---|---|---|
| `applyNotationProviderForCyclized` | **`applyNotationProviderForHarmonizedSequence`** | Tags column + attaches `CyclizedNotationProvider` |
| `harmonizedSequenceNotationProviderConstructor` | harmonizedSequenceNotationProviderConstructor | `meta.role: notationProviderConstructor`; returns the provider class |

#### Demo Entries

| Method | Registered name | demoPath |
|---|---|---|
| `demoTranslateSequence` | demoOligoTranslator | `Bioinformatics \| Oligo Toolkit \| Translator` |
| `demoOligoPattern` | demoOligoPattern | `Bioinformatics \| Oligo Toolkit \| Pattern` |
| `demoOligoStructure` | demoOligoStructure | `Bioinformatics \| Oligo Toolkit \| Structure` |

---

## detectors.js — Bootstrap & Context-Menu Wiring

`SequenceTranslatorPackageDetectors extends DG.Package`. **No live `semTypeDetector` is registered** — only an autostart bootstrap and a `notationRefiner`.

### `autostart()` (`tags: autostart, meta.role: autostart`)
Subscribes once to `grok.events.onContextMenu` and dispatches by `tableColumn.semType`.

For `DG.GridCell` items with a `tableColumn`:

- **`semType === Macromolecule`**:
  - top-level item **`PolyTool-Enumerate`** → calls `${pkg}:getPtHelmEnumeratorDialog`
  - submenu **`Oligo`**:
    - `Convert HELM to Oligo` (only when `col.meta.units === 'helm'`) → `${pkg}:convertHelmToOligoNucleotide`
    - `Combine sense+antisense to Oligo...` → `DG.Func.find(...).prepare({table, senseCol: col}).edit()` for `combineSenseAntisenseToOligoNucleotide`
- **`semType === Molecule`**:
  - `PolyTool-Enumerate` → `${pkg}:getPtChemEnumeratorDialog`
- **`semType === 'OligoNucleotide'`**:
  - `PolyTool-Enumerate` → `${pkg}:getPtOligoEnumeratorDialog`

Errors are stashed on `window.$sequenceTranslator.contextMenuError` for postmortem.

### `refineNotationProviderForHarmonizedSequence` (`meta.role: notationRefiner, tags: notationRefiner`)
Auto-detects cyclized (`^.+\(\d+\)$`) and dimerized (`^\(#\d\).+$`) sequences from `stats.freq` and applies the cyclized notation provider via `${pkg}:applyNotationProviderForCyclized`.

A commented-out `detectHarmonizedSequence` is a placeholder; not active.

---

## Module Details

### `src/oligo-renderer/` — OligoNucleotide subsystem

A self-contained renderer + panels + converters for siRNA/ASO duplex HELM cells. Cell value is HELM under the hood; the `OligoNucleotide` semType + `quality=OligoNucleotide` tag selects the duplex renderer.

| File | Purpose |
|---|---|
| `types.ts` | `OLIGO_SEM_TYPE` / `OLIGO_UNITS` constants, `ParsedNucleotide` / `ParsedConjugate` / `ParsedStrand` / `ParsedDuplex` interfaces, modification dictionary (`SUGAR_MODS`, `PHOSPHATE_MODS`, `CONJUGATE_MODS`), HELMCore-canonical alias maps (`SUGAR_ALIASES`, `PHOSPHATE_ALIASES`), color resolvers, hash-color fallback for unknowns |
| `helm-parser.ts` | `parseHelmDuplex(helm)` returns `{sense, antisense, raw, pairs?, strandTypes?, alignment?}`. Parses `RNA1{...}\|RNA2{...}$conn$$ann$V2.0` shape with bracketed/unbracketed monomers; base-less monomers are classified as **`linker`** (known phosphate symbols, via `isLinkerSymbol`) or **`conjugate`**. Reads the `$connections$` `pair` section (`parseHelmConnections`) and extended `strandtype` annotations (`parseHelmAnnotations`): swaps chains so `sense` is the real sense strand, maps HELM monomer positions → nucleotide indices, and derives the explicit alignment shift. Also `canonicalizeHelm()` / `canonicalizeChainBody()` for rewriting aliased symbols (`mR`, `fR`, `LR`, `sP`, …) to HELMCore canonical. Custom (not bio's `HelmHelper.parse`) for hot-path performance |
| `alignment.ts` | Resolves the sense/antisense display register. `resolveDuplexAlignment(model)` returns `{shift, source}` — explicit HELM pairs win; else `computeAutoShift` slides the reversed antisense across sense and picks the column shift with the most Watson-Crick complementary pairs (`bestComplementaryShift`), reducing custom bases via natural analog. A fraction/min-count/gain guard keeps weakly-complementary inputs blunt. `columnOffsets(shift)`, `describeOverhangs(model)` (5'/3' overhang lengths per strand) for layout + panels. `isLinkerMonomer(m)` classifies a monomer as a backbone linker (parser `linker` tag, or a conjugate whose natural analog is `p`) |
| `canvas-renderer.ts` | Pure functions: `computeLayout(cellW, cellH, model, opts)`, `drawDuplex()`, `hitTest()`. `buildGrid` lays sense/antisense on a shared **base**-column grid; the alignment `shift` offsets one strand so **overhangs on either end render true** (non-zero shift forces uniform chip widths so columns stay exact). Standalone linkers render as **arcs (no chip)** and widen the inter-base gap they sit in — a run on one strand the other lacks opens a gap (bases stay paired) while the partner draws one wider gap-spanning arc (`distributeArcs`). `ChipPos.col` is the shared column index; base pairings drawn per shared column (overhang columns simply unpaired). Variable chip widths (conjugates), antisense reversal, leading-extra shift, sugar-mod stripe |
| `cell-renderer.ts` | `OligoNucleotideCellRenderer extends DG.GridCellRenderer`. Per-value parse cache, per-cell layout cache keyed by `colName@col.version::rowIdx` (column-version invalidates on edits), `onMouseMove` resolves hover via `hitTest` → tooltip |
| `tooltip.ts` | `showMonomerTooltip(hit, x, y)`. Builds details synchronously, async-loads sugar / base / 3'-linkage RDKit structures from the **central Bio monomer library** (`getMonomerLibHelper().getMonomerLib()`), with alias resolution and structure caching (`Map<kind:canonicalSymbol, HTMLElement>`, reused across hovers) |
| `legend-panel.ts` | `buildOligoPanel(value)` widget — per-cell sense/antisense lengths, modification counts (collapsed by canonical symbol), conjugates, color legend filtered to mods actually present in the cell. Standalone linkers are counted as linkages (via `isLinkerMonomer`), not conjugates |
| `structures-panel.ts` | `buildOligoStructuresPanel(value)` widget — `splitDuplexToStrandStructures(helm)` splits the duplex into per-strand standalone HELMs **ordered + labeled by the parsed model** (so `strandtype` swaps put the real sense strand in the "Sense" pane, not just chain-0), runs `canonicalizeHelm` on each, places into a temp Macromolecule DataFrame, wraps each row's cell as `DG.SemanticValue.fromTableCell(cell)` and forwards to `Bio:toAtomicLevelPanel` inside lazy `ui.accordion()` panes ("Sense" / "Antisense") |
| `converters.ts` | `tagAsOligoNucleotide(col)`, `convertHelmColumnToOligo(table, helmCol)` (creates a tagged clone), `combineSenseAntisenseToOligo(table, senseCol, antiCol)` (renumbers each chain to RNA1{…}, joins with `\|`) |

### `apps/common/model/` — Core Domain Model

**`oligo-toolkit-package.ts`** — `OligoToolkitPackage` (extends `DG.Package`, implements `ITranslationHelper`)
- Central facade managing initialization, monomer libraries, JSON data
- Factory for `SequenceValidator`, `FormatConverter`, `FormatDetector`
- Lazy init with promise caching

**`data-loader/json-loader.ts`** — `JsonData` class — Loads (in parallel from `MonomersPath`):
- `monomer-lib.json`, `formats-to-helm.json`, `codes-to-symbols.json`, `pattern-app-data.json`, `linkers.json`

**`monomer-lib/lib-wrapper.ts`** — `MonomerLibWrapper` — Adapter between `IMonomerLib` (from bio library) and the application; lookups by symbol/format/code; molecular-weight aggregation.

**`parsing-validation/format-detector.ts`** — `FormatDetector` — infers format from content (HELM prefix first, then code-scan).

**`parsing-validation/format-handler.ts`** — `FormatHandler` — bidirectional format ↔ HELM mapping with regex/negative-lookaround.

**`parsing-validation/sequence-validator.ts`** — `SequenceValidator` — greedy longest-match validation.

**`helpers.ts`** — `sortByReverseLength()`, `download()`, `tryCatch()`.
**`const.ts`** — `NUCLEOTIDES`, `TECHNOLOGIES`, `DEFAULT_FORMATS` enum.

### `apps/common/view/` — Shared UI

`AppUIBase`, `CombinedAppUI`, `IsolatedAppUIBase`, `MonomerLibViewer`, `ColoredTextInput` (textarea with syntax highlighting via `input-painters.ts`), `MoleculeImage` (canvas V3000 renderer), `draw-molecule.ts`, URL `router.ts`, `app-info-dialog.ts`.

### `apps/translator/` — Format Conversion App

**Model:** `format-converter.ts` (`FormatConverter` — HELM-pivot conversion), `conversion-utils.ts` (`getTranslatedSequences`, `getNucleotidesSequence`), `const.ts`.

**View:** `ui.ts` (`OligoTranslatorUI`, `TranslatorAppLayout`) — single + bulk modes, debounced input, EventBus state, SDF/SMILES export.

### `apps/pattern/` — Pattern Design App

**Model:** `event-bus.ts` (RxJS BehaviorSubjects), `data-manager.ts` (singleton CRUD via user storage, SHA1 identity), `translator.ts` (`bulkTranslate`, `applyPatternToRawSequence`), `router.ts` (URL ↔ EventBus), `subscription-manager.ts`, `types.ts`, `const.ts` (strands, termini, `MAX_SEQUENCE_LENGTH = 34`), `utils.ts`.

**View:** `ui.ts` (`OligoPatternUI`), `components/` (left-section, right-section, edit-block, load-block, strand-editor dialog, terminal-modification editor, bulk-convert, numeric-label visibility, translation-examples).

**SVG Rendering** (`view/svg-utils/`): `svg-block-base.ts` (`SVGBlockBase` — base class for composable SVG blocks, see Design Pattern #7), `svg-renderer.ts` (`NucleotidePatternSVGRenderer`), `strands-block.ts`, `title-block.ts`, `legend-block.ts`, `svg-element-factory.ts`, `svg-display-manager.ts` (debounced 100ms + PNG export), `text-dimensions-calculator.ts`, `utils.ts`, `const.ts`.

### `apps/structure/` — Molecular Structure App

**Model:** `sequence-to-molfile.ts` (`SequenceToMolfileConverter`), `monomer-code-parser.ts` (greedy longest-match, auto-inserts phosphate linkers), `mol-transformations.ts` (`linkStrandsV3000`, `getNucleotidesMol`, rotation/reflection helpers), `oligo-structure.ts` (`getMolfileForStrand`, `getLinkedMolfile`, `saveSdf`).

**View:** `ui.ts` (`OligoStructureUI`, `StructureAppLayout`) — three strand inputs, chirality toggle, SDF save, 650×150 canvas with 300ms debounce.

### `polytool/` — Advanced Polymer Conversion Subsystem

#### Top-level

`const.ts`, `types.ts` (`PolyToolEnumeratorParams`, `PolyToolEnumeratorTypes`), `utils.ts` (`_setPeptideColumn`, helpers).

#### Core Conversion (`conversion/`)

`pt-chain.ts` (`Chain`), `pt-conversion.ts` (`doPolyToolConvert`), `pt-tools-parse.ts` (`parseSeparator`, `parseHelm`, `fromObjectsToHelm`, `handleDuplicated`, `handleLinkRules`, `handleReactionRules`), `pt-tools-helmmol.ts` (`getHelmMol`, `helmMolToNotation`), `pt-atomic.ts` (`helmToMol`), `pt-synthetic.ts` (`getOverriddenLibrary`, RDKit reactions), `pt-rules.ts` (`Rules`, `RuleInputs`, `RuleLink`, `RuleReaction`), `pt-rule-cards.ts`, `rule-manager.ts`, `rule-reaction-editor.ts`, `pt-misc.ts` (`Linkage`, index helpers), `style.css`.

#### Dialogs

| File | Purpose |
|---|---|
| `pt-dialog.ts` | Main convert dialog; `polyToolConvertUI()`, `polyToolConvert()` |
| `pt-enumerate-seq-dialog.ts` | HELM enumeration dialog; `polyToolEnumerateHelmUI(cell?, outputAsOligo=false)`, `getPolyToolEnumerateDialog`, `polyToolEnumerateSeq`. The `outputAsOligo` flag tags the result column as OligoNucleotide so the duplex renderer picks it up |
| `pt-chem-enum.ts` | Chem enumeration core (formerly `pt-enumeration-chem.ts` in older docs). Pure logic — RDKit module only, no Datagrok/UI deps. Normalizes 5 R-label spellings (`[N*]`, `[*:N]`, `[*N]`, `[RN]`, `[R:N]`) to `[*:N]`. `Zip` and `Cartesian` modes; `CHEM_ENUM_MAX_RESULTS = 1_000_000` cap. **Two-path SMILES assembly** (`buildJoinedSmiles`): atom-substitution R-groups (no `[*:N]` in their SMILES — single-atom mode) are spliced into the core via plain string replace; labeled R-groups are then joined via shared ring-closure digits across a disconnected SMILES. **`isSingleAtom`** on `ChemEnumRGroup` is set when `makeRGroup` finds 0 R-labels and `trySingleAtomCanonical` (RDKit `get_num_atoms(true) === 1` + `remove_hs_in_place()`) returns a canonical atom token. Pipeline is two-stage: `enumerateRaw` returns parseable-but-uncanonical SMILES with **zero RDKit calls**; `executeEnumeration` then canonicalizes the whole column in one batched `Chem:convertNotation` call, which also fixes aromaticity case after a single-atom splice. `enumerate` / `enumerateSample` are sync-RDKit variants for tests. Helpers: `moveStartRLabelToBranch`, `pickFreeRingDigits`, `formatRingDigit`, `substituteRLabelWithRingDigit`, `substituteRLabelWithAtom` |
| `pt-chem-enum-dialog.ts` | Chem enumeration dialog UI; `polyToolEnumerateChemUI(cell?)` (dialog; cell preload → first core when it has ≥1 R-label, else `grok.shell.info`) and `polyToolEnumerateChemApp` (`DG.View` entry; both share `buildChemEnumPanel` + `bindActionButton`). Card UI 110×104 with 320×260 hover tooltip, three hover-revealed action icons: edit + duplicate (`var(--blue-3)`) + delete (`var(--red-3)`). **Duplicate** reopens the sketcher pre-loaded with the molecule and pushes a NEW card on OK (original untouched). Cores rendered in a single horizontal `virtualView`; R-groups in one horizontal `virtualView` per R# inside a capped-height vertical scroll. Live preview: up to 12 reservoir-sampled molecules via `enumerateSampleRaw` (uncanonical SMILES, `grok.chem.drawMolecule` parses lazily per visible card). `openCoreSketchDialog` accepts when ≥1 R-label is present; `openRGroupSketchDialog` accepts when **(exactly 1 R-label) OR (single atom + target R# set)**. Import wizard (`openImportWizard`) picks Molecule column from any open table, infers R# from column name via `/r[\s_\-:]*(\d+)/i`, dedup checkbox, single-atom rows labeled `R{n} · atom`. Mode tooltip / status messages explain Zip vs Cartesian and the single-atom shortcut |
| `pt-enumeration-helm.ts` | `doPolyToolEnumerateHelm` engine — Single / Parallel / Matrix / Breadth strategies |
| `pt-combine-dialog.ts` | `getPTCombineDialog` — Cartesian combine |
| `pt-placeholders-input.ts` | Grid input for point placeholders |
| `pt-placeholders-breadth-input.ts` | Grid input for range placeholders |
| `pt-unrule.ts` / `pt-unrule-dialog.ts` | HELM → harmonized notation reversal |
| `pt-convert-editor.ts` | `PolyToolConvertFuncEditor` |

#### Handlers

`monomer-lib-handler.ts` (`PolyToolMonomerLibHandler`), `csv-to-json-monomer-lib-converter.ts` (`PolyToolCsvLibHandler`).

### `utils/` — Notation Providers & Error Handling

`cyclized.ts` (`CyclizedNotationProvider` implements `INotationProvider`), `dimerized.ts` (extends cyclized), `cell-renderer-cyclized.ts` (`CyclizedCellRendererBack` for the cyclized notation), `err-info.ts` (`defaultErrorHandler`).

### `plugins/` — External Plugin Integration

`mermade.ts` — `getExternalAppViewFactories()` loads MerMade synthesis plugin via `grok.functions.call()`.

### `demo/`

`demo-st-ui.ts` — demo functions for translator / pattern / structure.

---

## Tests

Test entry: `package-test.ts` → imports each test file under `src/tests/`. Test framework: `@datagrok-libraries/test`.

| File | Category | What it tests |
|---|---|---|
| `formats-to-helm.ts` | Formats to HELM / HELM to Formats | Bidirectional format ↔ HELM conversion |
| `formats-support.ts` | Formats support | All expected formats available |
| `helm-to-nucleotides.ts` | HELM to Nucleotides | HELM → bare nucleotide string |
| `files-tests.ts` | files | CSV file-based integration tests |
| `polytool-chain-from-notation-tests.ts` | PolyTool: Chain | Chain parsing, HELM conversion, rule application |
| `polytool-chain-parse-notation-tests.ts` | PolyTool: Chain: parseNotation | Index translation, monomer counts |
| `polytool-convert-tests.ts` | PolyTool: Convert | End-to-end conversion pipeline |
| `polytool-detectors-custom-notation-test.ts` | PolyTool: detectors | Semantic-type detection |
| `polytool-enumerate-tests.ts` | PolyTool: Enumerate | Single / Parallel / Matrix enumeration |
| `polytool-enumerate-breadth-tests.ts` | PolyTool: Enumerate | Breadth range enumeration |
| `polytool-enumerate-chem-tests.ts` | PolyTool: Enumerate | Chem enumeration |
| `polytool-unrule-tests.ts` | PolyTool: Unrule | HELM → harmonized notation |
| `toAtomicLevel-tests.ts` | toAtomicLevel | Synthetic monomer generation via RDKit |
| `oligo-renderer-tests.ts` | OligoRenderer: parser / dictionary / layout / hit testing / drawing smoke / HELM connections & strand types / auto-alignment / alignment layout / complex HELMs | HELM parse, alias resolution, layout invariants, hit-test correctness, smoke drawing, explicit pair + strandtype parsing & strand swap, complementary auto-shift + overhang detection, shift-aware layout / column pairing, complex real-world duplexes (both-end conjugates, custom bases, ± / large / asymmetric shifts) |
| `oligo-fixtures.ts` | — | Shared HELM fixtures for the OligoNucleotide test suites (blunt, explicit-pair, overhang, swap, complex, and standalone-linker duplexes) |
| `oligo-linkers-tests.ts` | OligoRenderer: standalone linkers — parsing / single-strand rendering / duplex alignment & gaps | Linker vs conjugate classification, leading/middle/trailing/consecutive linkers render as arcs (no chips), and a sense-only linker bulge that widens the gap + draws a wider antisense arc while keeping bases aligned |
| `const.ts` | — | Test data: `formatsToHelm`, `helmToNucleotides` |
| `utils/` | — | Test helpers (`detect-macromolecule-utils.ts`, `index.ts`) |

---

## Data Files

```
files/
├── monomers-sample/
│   ├── monomer-lib.json
│   ├── codes-to-symbols.json
│   ├── formats-to-helm.json
│   ├── linkers.json
│   ├── pattern-app-data.json
│   └── README.md
├── polytool-rules/
│   └── rules_example.json
├── samples/
│   ├── HELM.csv
│   ├── cyclized.csv
│   ├── cyclized_MSA.csv
│   ├── bulk-translation-axolabs.csv
│   └── sirna-demo.csv             # 40+ siRNA / ASO duplexes: 2'-OMe, 2'-F, GalNAc-L3, Chol, LNA/MOE gapmers, plus 3' overhangs (auto-aligned), rows with explicit HELM pair/strandtype info, and standalone-linker rows (5' phosphate cap + mid-strand linker bulge). Overhang/explicit/linker rows patched in by prototypes/patch-sirna-demo.mjs
└── tests/
    ├── axolabs1.csv
    ├── polytool-reaction-lib.json
    ├── chem_enum_cores.csv
    ├── chem_enum_rgroups.csv
    └── README.md
```

> **Note**: the OligoNucleotide subsystem also relies on `packages/Bio/files/monomer-libraries/oligo-conjugates.json` (GalNAc / L3 / Chol / Bio / Toc / Pal / DBCO from PubChem). That file lives in the **Bio** package because it's a global monomer library merged by `MonomerLibFromFilesProvider` — not an ST-internal asset.

---

## `prototypes/` — Off-Tree Dev Tooling

Helper / diagnostic scripts for the OligoNucleotide subsystem. Pure Node (no deps), **not bundled** in the published package, register no Datagrok functions. Run with `node prototypes/<file>`.

| File | Purpose |
|---|---|
| `patch-sirna-demo.mjs` | Idempotently patches `files/samples/sirna-demo.csv`: adds 3' dTdT overhangs to a couple of existing blunt duplexes (auto-aligned) and appends rows with overhangs + explicit HELM `$connections$` pair info / `strandtype` annotations, plus standalone-linker rows (5' phosphate cap + mid-strand linker bulge). Quotes fields containing commas (HELM connection sections) |
| `verify-alignment.mjs` | Parses every duplex in the demo CSV, resolves the alignment, and prints an ASCII sense/pairing/antisense duplex so shifts + overhangs can be eyeballed |
| `scan-test-fixtures.mjs` | Reports the auto-align shift the renderer would pick for every inline `RNA1{…}\|RNA2{…}` fixture in the tests — used to confirm the auto-aligner leaves existing (blunt) fixtures untouched |
| `gen-complex-test-helms.mjs` | Builds the complex-HELM unit-test fixtures (both-end conjugates, custom bases, positive/large/asymmetric shifts, strand swap) — emits exact `$connections$` pair sections; pasted into `src/tests/oligo-fixtures.ts`. Also exports `FIXTURES` for probing |
| `probe-complex.mjs` | Runs the real TS over the `gen-complex-test-helms` fixtures and prints each one's resolved properties (shift, source, overhangs, lengths, conjugates) — so the unit-test expectations are written from observed reality |
| `run-real-ts.mjs` | Transpiles the real `types.ts` / `alignment.ts` / `helm-parser.ts` / `canvas-renderer.ts` (DG deps stubbed) and asserts parse + alignment + `computeLayout` shift behavior against the demo CSV and fixtures |

---

## Key Design Patterns

1. **HELM as Pivot Format** — All format conversions go through HELM as intermediate representation
2. **EventBus (RxJS BehaviorSubjects)** — Central state management in pattern app; decouples UI components
3. **Singleton Managers** — `DataManager`, `RulesManager`, `TextDimensionsCalculator`
4. **IsolatedAppUIBase** — Each sub-app is a standalone view that can be used in combined or isolated mode
5. **Lazy Initialization** — Promise caching for monomer library and JSON data loading
6. **Greedy Longest-Match** — Sequence parsing always tries longest codes first
7. **SVG Composable Blocks** — `TitleBlock`, `StrandsBlock`, `LegendBlock` extend `SVGBlockBase`
8. **Dialog Singletons** — Strand/terminal editors prevent multiple instances
9. **HELMCore-canonical aliases** (oligo-renderer) — input HELM may use vendor / Pistoia / legacy symbols (`mR`, `fR`, `sP`); the renderer resolves through alias maps and serializes canonical (`m`, `fl2r`, `sp`) before forwarding to Bio's library-driven pipelines (so monomer lookups always hit)
10. **Column-version-keyed cache** — the OligoNucleotide cell renderer's per-cell layout cache is keyed by `${colName}@${col.version}::${rowIdx}`, so column edits orphan old entries automatically
11. **Two-stage chem-enum pipeline** — `enumerateRaw` does pure-string assembly (ring-closure trick + atom splice for single-atom R-groups), produces parseable but uncanonical SMILES with **zero per-row RDKit calls**; `executeEnumeration` then canonicalizes the whole output column in a single batched `Chem:convertNotation` call (`overwrite: true`, `kekulize: false`). This also handles aromaticity case-fixup after a single-atom splice (e.g. uppercase `N` spliced into an aromatic ring). The bulk path is parallelized across Chem workers; per-row sync RDKit (`enumerate` / `enumerateSample`) is reserved for tests
12. **Single-atom R-group shortcut** (chem-enum) — drawing a bare atom (`N`, `O`, `Cl`) without any `[*:N]` label is accepted as long as RDKit confirms exactly one heavy atom. Such groups bypass ring-closure joining and are spliced into the core's `[*:N]` slot via plain string replace; mixing single-atom and labeled R-groups in the same enumeration is supported per R-number
13. **Alignment source-of-truth precedence** (oligo-renderer) — strand roles and base-pair register come from the HELM itself when present: `strandtype` annotations decide which chain is sense (chains are swapped so `sense` is always sense), and `$connections$` `pair` entries fix the alignment shift. Only when that info is absent does the renderer auto-align — a single rigid shift (not a gapped alignment) over `[-(antiLen-1), senseLen-1]` maximizing Watson-Crick complementary pairs, with custom bases reduced to their natural analog. A guard (≥3 pairs, ≥60% of the overlap, ≥2 more than blunt) keeps weakly-complementary / synthetic inputs blunt. The shift is a column offset on a shared nucleotide grid, so terminal overhangs render true on either end
14. **Standalone backbone linkers as arcs** (oligo-renderer) — base-less phosphate-type monomers (`p`, `[sp]`, or anything whose natural analog is `p`) are parsed as `linker` (not conjugate) and drawn as linkage arcs with no chip, anywhere they occur (5'/3' caps, mid-strand, consecutive runs). Only **bases** anchor the shared column grid; standalone linkers widen the inter-base gap they sit in. When one strand has a linker run the other lacks, both reserve the wider gap (paired bases stay aligned) and the strand without the run draws a single wider arc spanning the gap. Arcs bridge toward the chain core (a leading-linker arc sits in the gap to the right of its slot, a trailing one to its left, and an inter-base gap places the base phosphate tight after the chip then the linkers in a row) so the backbone reads continuously and arcs land where a normal phosphate would. The truncation "…" marker counts only chip-producing monomers (`chipCount`), so standalone linkers never trip it. Reduces to the plain chip grid when no standalone linkers are present (all prior behavior preserved)

---

## Key Data Flow

```
User Input (sequence string)
  → FormatDetector.getFormat()          # Infer source format
  → SequenceValidator.isValidSequence() # Validate
  → FormatConverter.convertTo(target)   # Convert via HELM pivot
  → Output (translated sequence)

PolyTool:
  Separator notation → parseSeparator() → Chain
  → Chain.applyRules() → handleLinkRules / handleReactionRules
  → Chain.getHelm() → HELM string
  → helmToMol() → V3000 molfile (optional)

Pattern:
  Raw sequence + PatternConfig
  → applyPatternToRawSequence() → modified nucleotides
  → SVGRenderer.renderPattern() → visual output

OligoNucleotide cell rendering:
  HELM cell value
  → parseHelmDuplex() → ParsedDuplex {sense, antisense, pairs?, strandTypes?, alignment?}
                        (explicit HELM pairs/strandtype → swap chains + alignment shift)
  → computeLayout(w, h, model) → resolveDuplexAlignment (explicit shift, else auto by complementarity)
                               → ChipPos[] on a shared column grid (shift offsets overhangs)
  → drawDuplex() → canvas (base pairs drawn per shared column; overhang columns unpaired)
  → onMouseMove → hitTest → showMonomerTooltip
                            → getBioLib() → findMonomerMolfile (canonical alias)
                            → grok.chem.drawMolecule → cached HTMLElement

OligoNucleotide structures panel:
  HELM cell value
  → split chains by '|', renumber each to RNA1{…}, canonicalizeHelm
  → temp DataFrame[Macromolecule helm column]
  → DG.SemanticValue.fromTableCell → Bio:toAtomicLevelPanel widget
  → ui.accordion(lazy panes) per strand
```

---

## Scripts

- `scripts/build-monomer-lib.py` — Python (RDKit) builds monomer-library JSON from source files
- `link-bio` — bash script to build and link `@datagrok-libraries/bio` locally

---

## Quick Lookups

| Looking for... | Check first |
|---|---|
| Function/panel/cell-renderer registration | `src/package.ts` (`PackageFunctions` class) |
| Context-menu wiring + `Oligo` submenu | `detectors.js` (`autostartContextMenu`) |
| OligoNucleotide cell renderer | `src/oligo-renderer/cell-renderer.ts` |
| OligoNucleotide layout / drawing | `src/oligo-renderer/canvas-renderer.ts` |
| Sense/antisense alignment, shift, overhangs | `src/oligo-renderer/alignment.ts` (`resolveDuplexAlignment`, `computeAutoShift`, `describeOverhangs`) |
| HELM pair / strandtype parsing | `src/oligo-renderer/helm-parser.ts` (`parseHelmConnections`, `parseHelmAnnotations`) |
| HELM canonicalization (alias → HELMCore) | `src/oligo-renderer/helm-parser.ts` (`canonicalizeHelm`) |
| Modification colors / aliases | `src/oligo-renderer/types.ts` |
| Sense/antisense structure panel | `src/oligo-renderer/structures-panel.ts` |
| Convert HELM → Oligo / combine sense+AS | `src/oligo-renderer/converters.ts` |
| Translator app | `src/apps/translator/` |
| Pattern designer app | `src/apps/pattern/` |
| Structure app | `src/apps/structure/` |
| PolyTool convert pipeline | `src/polytool/conversion/`, `pt-dialog.ts` |
| HELM enumeration (with `outputAsOligo`) | `src/polytool/pt-enumerate-seq-dialog.ts` |
| Chem enumeration (logic) | `src/polytool/pt-chem-enum.ts` (`buildJoinedSmiles` two-path assembly, `trySingleAtomCanonical`, `substituteRLabelWithAtom`) |
| Chem enumeration (UI) | `src/polytool/pt-chem-enum-dialog.ts` (`buildChemEnumPanel`, sketch dialogs, `openImportWizard`, duplicate/edit/delete card icons) |
| Demo / sample siRNA data | `files/samples/sirna-demo.csv` |
| Off-tree dev tooling | `prototypes/` (patch demo / verify alignment / scan fixtures / run real TS) |
| Auto-generated wrappers | `src/package.g.ts`, `src/package-api.ts` |
