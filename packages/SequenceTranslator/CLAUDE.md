# SequenceTranslator Package - CLAUDE.md

## Overview

**Package:** `@datagrok/sequence-translator` (v1.10.8)
**Category:** Bioinformatics
**Author:** Davit Rizhinashvili
**Description:** Translates oligonucleotide sequences between different representations (Axolabs, BioSpring, GCRS, Mermade, HELM, nucleotides, etc.). Also provides pattern design, structure visualization, and advanced polymer conversion (PolyTool).

## Quick Reference

```bash
npm run build          # grok api && grok check --soft && webpack
npm run build-all      # Build all deps (chem-meta → js-api → utils → bio → tutorials → this)
npm run link-all       # npm link chem-meta datagrok-api utils bio tutorials
npm run test           # grok test
npm run test-local     # grok test --host localhost
npm run lint           # eslint ./src/**/*.ts
npm run lint-fix       # eslint --fix
```

## Dependencies

| Dependency | Purpose |
|------------|---------|
| `datagrok-api` | Core Datagrok platform API |
| `@datagrok-libraries/bio` | HELM helper, sequence handler, macromolecule utils, monomer libraries |
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

---

## Architecture

The package is organized into **four main applications** + a **PolyTool subsystem** + shared infrastructure:

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
│   │   ├── model/                # Core domain logic
│   │   └── view/                 # Shared UI components
│   ├── translator/               # Format conversion app
│   │   ├── model/                # Conversion logic
│   │   └── view/                 # Translator UI
│   ├── pattern/                  # Pattern design app
│   │   ├── model/                # Pattern state, persistence, translation
│   │   └── view/                 # Pattern editor UI + SVG rendering
│   └── structure/                # Molecular structure app
│       ├── model/                # Molfile generation, V3000 manipulation
│       └── view/                 # Structure builder UI
│
├── polytool/                     # Advanced polymer conversion subsystem
│   ├── conversion/               # Core algorithms (chain, rules, reactions)
│   ├── *.ts                      # Dialogs, enumeration, handlers
│
├── plugins/                      # External plugin integration (MerMade)
├── utils/                        # Cyclized/dimerized notation providers, error handling
├── demo/                         # Demo entry points
└── tests/                        # Test suite
```

---

## Entry Point: package.ts

Central `PackageFunctions` class exposes all Datagrok-registered functions. Key instance:

```typescript
export const _package: OligoToolkitPackage = new OligoToolkitPackage({debug: true});
```

### Initialization Flow

1. `init()` → gets HELM helper from bio library, calls `completeInit()`
2. `completeInit()` → sets up monomer library wrapper
3. `initLibData()` → lazy-loads JSON data files (called on app open)

### Registered Functions

| Function | Type | Description |
|----------|------|-------------|
| `oligoToolkitApp()` | app | Combined tabbed app (all 3 sub-apps + plugins) |
| `oligoTranslatorApp()` | app | Standalone translator |
| `oligoPatternApp()` | app | Standalone pattern designer |
| `oligoStructureApp()` | app | Standalone structure viewer |
| `getTranslationHelper()` | func | Returns `ITranslationHelper` for external consumers |
| `validateSequence(seq)` | func | Boolean validation |
| `translateOligonucleotideSequence(seq, src, tgt)` | func | Format-to-format conversion |
| `getMolfileFromGcrsSequence(seq, invert)` | func | GCRS → V3000 molfile |
| `linkStrands(strands)` | func | Multi-strand assembly |
| `polyToolConvert2()` | func | Bulk PolyTool conversion |
| `polyToolEnumerateHelmTopMenu()` | dialog | HELM enumeration |
| `polyToolEnumerateChemTopMenu()` | dialog | Chemical enumeration |
| `enumerateSingleHelmSequence()` | func | Single HELM enumeration API |
| `enumerateSingleHelmSequenceWithNaturalAAs()` | func | Enumerate with natural amino acids |
| `getPolyToolCombineDialog()` | dialog | Combine sequences from columns |
| `applyNotationProviderForCyclized()` | func | Register custom notation |
| `createMonomerLibraryForPolyTool()` | func | CSV → JSON monomer library |

---

## Module Details

### apps/common/model/ — Core Domain Model

**`oligo-toolkit-package.ts`** — `OligoToolkitPackage` class (extends `DG.Package`, implements `ITranslationHelper`)
- Central facade managing initialization, monomer libraries, JSON data
- Factory for `SequenceValidator`, `FormatConverter`, `FormatDetector`
- Lazy init with promise caching

**`data-loader/json-loader.ts`** — `JsonData` class
- Loads 4 JSON files in parallel from `MonomersPath`:
  - `monomer-lib.json` — full monomer library with molfiles
  - `formats-to-helm.json` — format code → HELM mappings
  - `codes-to-symbols.json` — format code → monomer symbol mappings
  - `pattern-app-data.json` — color maps for pattern app
  - `linkers.json` — monomers with phosphate groups

**`monomer-lib/lib-wrapper.ts`** — `MonomerLibWrapper` class
- Adapter between `IMonomerLib` (from bio library) and the application
- Lookups: by symbol, by format, by code
- Molecular weight aggregation, DataFrame generation for viewer

**`parsing-validation/format-detector.ts`** — `FormatDetector` class
- Infers sequence format from content (checks HELM prefix, then scans for matching codes)
- Validates candidates using `SequenceValidator`

**`parsing-validation/format-handler.ts`** — `FormatHandler` class
- Bidirectional mapping between source formats and HELM notation
- Regex generation with negative lookaround for accurate code matching
- Phosphate linkage extraction

**`parsing-validation/sequence-validator.ts`** — `SequenceValidator` class
- Greedy longest-match-first validation
- Returns index of first invalid code or -1 for valid

**`helpers.ts`** — Utility functions: `sortByReverseLength()`, `download()`, `tryCatch()`

**`const.ts`** — `NUCLEOTIDES` array `['A','G','C','U']`, `TECHNOLOGIES`, `DEFAULT_FORMATS` enum

### apps/common/view/ — Shared UI

**`app-ui-base.ts`** — `AppUIBase` abstract class with loading progress indicators

**`combined-app-ui.ts`** — `CombinedAppUI` — Multi-tab view (Translator + Pattern + Structure + external plugins)

**`isolated-app-ui.ts`** — `IsolatedAppUIBase` — Container for standalone app UIs

**`monomer-lib-viewer.ts`** — `MonomerLibViewer` — Interactive monomer library table with molecule rendering

**`components/colored-input/`** — `ColoredTextInput` — Textarea with syntax highlighting overlay
- `input-painters.ts` — `highlightInvalidSubsequence()` painter (red for invalid portions)

**`components/molecule-img.ts`** — `MoleculeImage` — Canvas-based V3000 molfile renderer

**`components/draw-molecule.ts`** — Molecule drawing utilities

**`components/router.ts`** — URL routing for combined app tabs

**`components/app-info-dialog.ts`** — App metadata dialog

---

### apps/translator/ — Format Conversion App

**Model:**
- `format-converter.ts` — `FormatConverter` class — Converts between formats using HELM as pivot
  - `convertTo(targetFormat)` → helmToFormat / formatToHelm / format→HELM→format
- `conversion-utils.ts` — `getTranslatedSequences()` — Translates to ALL supported formats at once
  - `getNucleotidesSequence()` — Extracts nucleotides from HELM
  - `convert()` — Legacy single-conversion wrapper
- `const.ts` — `GROUP_TYPE` (NUCLEOSIDE/LINKAGE), `PHOSPHATE_SYMBOL` ('p'), `UNKNOWN_SYMBOL` ('<?>')

**View:**
- `ui.ts` — `OligoTranslatorUI` + `TranslatorAppLayout`
  - **Single mode:** Text input → format detection → all-format output table with copy buttons
  - **Bulk mode:** Table/column selection → batch conversion → new column added
  - `EventBus` (local) — RxJS BehaviorSubjects for table/column/format state
  - Debounced input (300ms), real-time validation via `ColoredTextInput`
  - SDF export and SMILES conversion support

---

### apps/pattern/ — Pattern Design App

**Model:**
- `event-bus.ts` — `EventBus` — Central RxJS state manager for entire pattern app
  - BehaviorSubjects: pattern name, nucleotide sequences, PTO flags, terminal modifications, table selection, etc.
  - Computed observables: `patternStateChanged$`, `strandsUpdated$`, `uniqueNucleotidesChanged$()`
  - Change tracking (last-loaded config vs current)
- `data-manager.ts` — `DataManager` (singleton) — Pattern CRUD via Datagrok user storage
  - SHA1 hashing for pattern identity
  - Current user vs other users pattern management
  - Pattern uniqueness validation
- `translator.ts` — `bulkTranslate()` — Applies pattern to table columns
  - `applyPatternToRawSequence()` — Single-sequence pattern application
- `router.ts` — `URLRouter` — URL ↔ EventBus synchronization (shareable pattern links)
- `subscription-manager.ts` — `SubscriptionManager` — RxJS subscription cleanup
- `types.ts` — `PatternConfiguration`, `PatternConfigRecord`, error classes
- `const.ts` — Strands (SENSE/ANTISENSE), termini (5'/3'), `MAX_SEQUENCE_LENGTH = 34`, storage name
- `utils.ts` — Nucleotide analysis, strand truncation/extension helpers

**View:**
- `ui.ts` — `OligoPatternUI` — Main orchestrator (DI wiring, layout composition)
  - Left section: Load/Edit controls, table selection for bulk conversion
  - Right section: SVG visualization, save/download/share buttons
- `components/left-section.ts` — Load controls, edit controls, bulk convert controls
- `components/edit-block-controls.ts` — Strand length, nucleobase choice, name/comment editing
- `components/load-block-controls.ts` — Author/pattern dropdowns, delete button
- `components/strand-editor/dialog.ts` — Per-position nucleotide + PTO editing
- `components/terminal-modification-editor.ts` — 5'/3' terminal modification editing
- `components/bulk-convert/` — Table + column selection for batch pattern application
- `components/numeric-label-visibility-controls.ts` — Toggle numeric labels per nucleotide
- `components/translation-examples-block.ts` — Live input→output examples

**SVG Rendering (`view/svg-utils/`):**
- `svg-renderer.ts` — `NucleotidePatternSVGRenderer` — Orchestrates SVG generation
  - Composes: `TitleBlock` + `StrandsBlock` + `LegendBlock`
- `strands-block.ts` — Nucleotide circles, PTO stars, terminal mods, strand labels
- `title-block.ts` — Pattern name + strand lengths
- `legend-block.ts` — Color-coded legend with PTO indicator
- `svg-element-factory.ts` — SVG element creation (circles, text, stars, rectangles)
- `svg-display-manager.ts` — Debounced (100ms) SVG updates + PNG export
- `text-dimensions-calculator.ts` — Canvas-based text measurement
- `const.ts` — Radii, fonts, colors, Y-positions, dimension constants

---

### apps/structure/ — Molecular Structure App

**Model:**
- `sequence-to-molfile.ts` — `SequenceToMolfileConverter` — Parses sequence → retrieves monomer molfiles → links → V3000 molfile
- `monomer-code-parser.ts` — `MonomerSequenceParser` — Greedy longest-match parsing, auto-inserts phosphate linkers
- `mol-transformations.ts` — V3000 molfile manipulation:
  - `linkStrandsV3000()` — Links sense/antisense strands with coordinate transformation
  - `getNucleotidesMol()` — Combines nucleotide molblocks with linkage
  - Rotation, reflection, inversion functions for strand alignment
  - Atom/bond renumbering, stereo configuration handling
- `oligo-structure.ts` — High-level API:
  - `getMolfileForStrand()` — Single strand molfile
  - `getLinkedMolfile()` — Multi-strand assembly
  - `saveSdf()` — SDF file download with metadata

**View:**
- `ui.ts` — `OligoStructureUI` + `StructureAppLayout`
  - Three strand inputs (sense, antisense, antisense2) with direction choice (5'→3' / 3'→5')
  - Chirality toggle, SDF save, real-time molecule visualization
  - Canvas rendering (650x150) with 300ms debounce

---

### polytool/ — Advanced Polymer Conversion Subsystem

The PolyTool handles complex polymer sequence operations: notation conversion, rules-based transformations, enumeration, and chemical synthesis.

#### Core Conversion (`conversion/`)

**`pt-chain.ts`** — `Chain` class — Central data model
- `fromSeparator(notation)` — Parse separator notation (e.g., "A-B-{C}-D")
- `fromHelm(helm)` — Parse HELM notation
- `getNotation()` — Output harmonized notation
- `getHelm()` — Output HELM string
- `applyRules(rules)` — Apply link/reaction rules, create position mappings

**`pt-conversion.ts`** — `doPolyToolConvert()` — Main conversion engine
- Parses separator notation → applies rules → outputs HELM
- Returns `[helms[], isLinear[], positionMaps[]]`

**`pt-tools-parse.ts`** — Parsing functions
- `parseSeparator()` — Handles inline chains with curly-brace fragments
- `parseHelm()` — Extracts peptide definitions and linkages
- `fromObjectsToHelm()` — Serializes monomers + linkages → HELM
- `handleDuplicated()` — Homo/hetero dimer expansion
- `handleLinkRules()` / `handleReactionRules()` — Rule application

**`pt-tools-helmmol.ts`** — `getHelmMol()` / `helmMolToNotation()` — Bidirectional HELM mol ↔ notation

**`pt-atomic.ts`** — `helmToMol()` — HELM → molfile with linearization and chirality

**`pt-synthetic.ts`** — `getOverriddenLibrary()` — RDKit-based reaction execution for synthetic monomers

**`pt-rules.ts`** — `Rules` class, `RuleInputs` class — Rule data model and file I/O
- `RuleLink` type — Linkage rules (monomer pairs, R-groups)
- `RuleReaction` type — Synthesis rules (monomer pairs, SMARTS reactions)

**`pt-rule-cards.ts`** — `RuleCards` — Visual rule preview with monomer card gallery

**`rule-manager.ts`** — `RulesManager` (singleton) — UI-driven rule editing and persistence

**`rule-reaction-editor.ts`** — Reaction editing dialog with molecule sketchers

**`pt-misc.ts`** — `Linkage` type, `getInnerIdx()` / `getOuterIdx()` — Index translation for multi-chain

#### Dialogs

**`pt-dialog.ts`** — Main conversion dialog and chem enumeration dialog
- `getPolyToolConvertDialog()` — Column selector, flags, rule file, history
- `polyToolConvert()` — Execution: validates → applies rules → generates HELM → converts to molfiles

**`pt-enumerate-seq-dialog.ts`** — HELM enumeration dialog
- `getPolyToolEnumerateDialog()` — HELM editor, placeholder grids, enumeration type selector
- `polyToolEnumerateSeq()` — Execution with single/parallel/matrix/breadth strategies

**`pt-enumeration-helm.ts`** — `doPolyToolEnumerateHelm()` — Core enumeration engine
- Single, Parallel, Matrix, Breadth strategies

**`pt-enumeration-chem.ts`** — `getEnumerationChem()` — Chemical library enumeration via RDKit

**`pt-combine-dialog.ts`** — `getPTCombineDialog()` — Cartesian product of sequences from columns

**`pt-monomer-selection-dialog.ts`** — Autocomplete monomer picker with tag-based display

**`pt-placeholders-input.ts`** — Grid input for point placeholders (position + monomers)

**`pt-placeholders-breadth-input.ts`** — Grid input for range placeholders (start/end + monomers)

**`pt-unrule.ts` / `pt-unrule-dialog.ts`** — HELM → harmonized notation reversal

**`pt-convert-editor.ts`** — `PolyToolConvertFuncEditor` — Function call editor

#### Handlers

**`monomer-lib-handler.ts`** — `PolyToolMonomerLibHandler` — Validates and converts monomer libraries
**`csv-to-json-monomer-lib-converter.ts`** — `PolyToolCsvLibHandler` — CSV → JSON conversion

---

### utils/ — Notation Providers & Error Handling

**`cyclized.ts`** — `CyclizedNotationProvider` (implements `INotationProvider`)
- Custom separator-based notation for cyclized peptides
- Converts to HELM via `Chain.fromSeparator()`

**`dimerized.ts`** — `DimerizedNotationProvider` (extends `CyclizedNotationProvider`)

**`cell-renderer-cyclized.ts`** — `CyclizedCellRendererBack` — Grid cell rendering for cyclized sequences

**`err-info.ts`** — `defaultErrorHandler()` — Centralized error logging

### plugins/ — External Plugin Integration

**`mermade.ts`** — `getExternalAppViewFactories()` — Loads MerMade synthesis plugin dynamically via `grok.functions.call()`

### demo/

**`demo-st-ui.ts`** — Demo functions for all three apps (translator, pattern, structure)

---

## Tests

Test entry: `package-test.ts` → `tests/` directory

| Test File | Category | What It Tests |
|-----------|----------|---------------|
| `formats-to-helm.ts` | Formats to HELM / HELM to Formats | Bidirectional format↔HELM conversion |
| `formats-support.ts` | Formats support | All expected formats available |
| `helm-to-nucleotides.ts` | HELM to Nucleotides | HELM → nucleotide string |
| `files-tests.ts` | files | CSV file-based integration tests |
| `polytool-chain-from-notation-tests.ts` | PolyTool: Chain | Chain parsing, HELM conversion, rule application |
| `polytool-chain-parse-notation-tests.ts` | PolyTool: Chain: parseNotation | Index translation, monomer counts |
| `polytool-convert-tests.ts` | PolyTool: Convert | End-to-end conversion pipeline |
| `polytool-detectors-custom-notation-test.ts` | PolyTool: detectors | Semantic type detection |
| `polytool-enumerate-tests.ts` | PolyTool: Enumerate | Single/Parallel/Matrix enumeration |
| `polytool-enumerate-breadth-tests.ts` | PolyTool: Enumerate | Breadth range enumeration |
| `polytool-unrule-tests.ts` | PolyTool: Unrule | HELM → harmonized notation |
| `toAtomicLevel-tests.ts` | toAtomicLevel | Synthetic monomer generation via RDKit |
| `const.ts` | — | Test data: formatsToHelm, helmToNucleotides |
| `utils/` | — | Test helpers: detection utils, format conversion |

---

## Data Files

```
files/
├── monomers-sample/
│   ├── monomer-lib.json          # Monomer library (molfiles, SMILES, metadata)
│   ├── codes-to-symbols.json     # Format code → monomer symbol mappings
│   ├── formats-to-helm.json      # Format code → HELM notation mappings
│   ├── linkers.json              # Monomers with phosphate groups
│   ├── pattern-app-data.json     # Color maps for pattern app SVG
│   └── README.md
├── polytool-rules/
│   └── rules_example.json        # Link/reaction rule definitions
├── samples/
│   ├── HELM.csv                  # HELM format examples
│   ├── cyclized.csv              # Cyclized peptide samples
│   ├── cyclized_MSA.csv          # MSA of cyclized variants
│   └── bulk-translation-axolabs.csv
└── tests/
    ├── axolabs1.csv              # Axolabs format test data
    ├── polytool-reaction-lib.json # Reaction library for tests
    └── README.md
```

## Key Design Patterns

1. **HELM as Pivot Format** — All format conversions go through HELM as intermediate representation
2. **EventBus (RxJS BehaviorSubjects)** — Central state management in pattern app; decouples UI components
3. **Singleton Managers** — `DataManager`, `RulesManager`, `TextDimensionsCalculator`
4. **IsolatedAppUIBase** — Each sub-app is a standalone view that can be used in combined or isolated mode
5. **Lazy Initialization** — Promise caching for monomer library and JSON data loading
6. **Greedy Longest-Match** — Sequence parsing always tries longest codes first
7. **SVG Composable Blocks** — `TitleBlock`, `StrandsBlock`, `LegendBlock` extend `SVGBlockBase`
8. **Dialog Singletons** — Strand/terminal editors prevent multiple instances

## Key Data Flow

```
User Input (sequence string)
  → FormatDetector.getFormat()          # Infer source format
  → SequenceValidator.isValidSequence() # Validate
  → FormatConverter.convertTo(target)   # Convert via HELM pivot
  → Output (translated sequence)

PolyTool:
  Separator notation → parseSeparator() → Chain
  → Chain.applyRules() → handleLinkRules/handleReactionRules
  → Chain.getHelm() → HELM string
  → helmToMol() → V3000 molfile (optional)

Pattern:
  Raw sequence + PatternConfig
  → applyPatternToRawSequence() → modified nucleotides
  → SVGRenderer.renderPattern() → visual output
```

## Scripts

- `scripts/build-monomer-lib.py` — Python script (RDKit) to build monomer library JSON from source files
- `link-bio` — Bash script to build and link `@datagrok-libraries/bio` locally
