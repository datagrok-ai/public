# CLAUDE.md

## Overview

**Helm** (`@datagrok/helm`) is a Datagrok plugin providing full support for **HELM notation** (Hierarchical Editing Language for Macromolecules) — importing, detecting, rendering, editing, and conversion. It is built on the **standalone Datagrok HELM Web Editor** (`@datagrok-libraries/hwe`), integrating it into the platform with a custom grid cell renderer, an input widget, a property panel, and HELM↔other-notation conversion.

The legacy forked Pistoia HELM Web Editor + Scilligence JSDraw2 + Dojo toolkit stack has been **fully removed** (HWE migration, package v3.0). There are no more bundled `vendor/dojo-1.10.10/` or `helm/JSDraw/` sources and no `org.helm.webeditor` / `JSDraw2` / `scil` / `dojo` globals.

Category: **Bioinformatics**.

## The hwe library — `@datagrok-libraries/hwe`

Every editor / parse / render / calculation path routes through the self-contained `@datagrok-libraries/hwe` package. It is an **external npm dependency** (not under this repo's `libraries/`) and injects no global scripts. The package consumes from it:

| Export | Purpose |
|---|---|
| `HelmService` | Parses HELM, renders **synchronously** to a canvas (`renderToCanvas`), and computes molfiles / molecular properties. Reads monomers through `bridgeMonomerLib(...)` |
| `HelmHelperAdapter` | Wraps a `HelmService` to expose Datagrok's `IHelmHelper` plus a legacy-compatible Pistoia `App` shape (`createWebEditorApp`), so existing dialog consumers swap in unchanged |
| `bridgeMonomerLib(lib)` | Adapts Bio's `IMonomerLib` so the editor resolves exactly the same monomers as the rest of the platform |
| `buildHelmMolView(mol)` | Builds a legacy-shaped (`atoms[].p` in model space) laid-out molecule used for hover hit-testing and the highlight overlay |
| `PolymerTypes`; types `IMonomerLibBaseLike`, `MonomersFuncsLike`, `CanvasBounds` | Constants / type shims for the adapter and renderer |
| `@datagrok-libraries/hwe/styles.css` | The editor's `.hw-*` stylesheet. **Must be imported** (hwe ships CSS as a bundler import only — no runtime injection) or the editor renders unstyled. Imported once in `src/package.ts` |

## Architecture

### Entry Point — `src/package.ts`

Imports `@datagrok-libraries/hwe/styles.css`. `PackageFunctions` class registers all platform-visible functions via `@grok.decorators`:
- `initHelm` — `@init`: loads the RDKit module, `SeqHelper`, `MonomerLibHelper`; creates the `HelmHelper` singleton (passing `seqHelper`, logger, `rdKitModule`); then `_package.completeInit(helmHelper, libHelper)`. No Dojo/editor bootstrapping or Pistoia monomer patching.
- `helmCellRenderer` — `@func` (role: `cellRenderer`): returns `HelmGridCellRenderer` for columns with `quality=Macromolecule, units=helm`
- `editMoleculeCell` — `@func` (role: `cellEditor`): opens the full-screen HELM editor app dialog for cell editing (`openWebEditor`)
- `openEditor` (`Edit Helm...`) — `@func` (action): context-menu editor for Macromolecule columns; converts non-HELM notations to HELM first, then writes the edited value back in the column's original notation
- `propertiesWidget` (`Properties`) — `@panel`: molecular formula, weight, extinction coefficient for a Macromolecule value
- `getMolfiles` — `@func`: converts a Macromolecule column to molfile strings via `helmHelper.getMolfiles()`
- `helmInput` — `@func` (role: `valueEditor`): creates the `HelmInput` widget for Macromolecule/HELM semType
- `getHelmHelper` — `@func`: returns the `IHelmHelper` singleton
- `measureCellRenderer` — `@func`: cell-render performance benchmark
- `highlightMonomers` (`Highlight Monomers`) — `@func`: test app that loads peptide data and calls `Bio:toAtomicLevel` with `highlight: true`

> The legacy `getHelmService` function is **gone** — the off-screen render service now lives inside hwe and is owned by the helper's adapter.

**Exports**: `_package` (instance of `HelmPackage`), auto-generated function wrappers from `package.g.ts`

### Detector — `detectors.js`

`HelmPackageDetectors` class:
- `autostart` (meta.role: autostart) — Helm bootstrap; registers the `ui.input.helmAsync` polyfill (delegates to `Helm:helmInput`) so other packages can create HELM inputs before this package fully loads. No semantic-type detector — HELM detection lives in Bio's `detectMacromolecule`.

### Package Singleton — `src/package-utils.ts`

`HelmPackage` class (extends `DG.Package`) — the package singleton managing:
- `helmHelper` / `seqHelper` — accessors for the initialized helpers
- `_libHelper` — monomer library helper; `monomerLib` — accessor for the current `IMonomerLib`
- `completeInit(helmHelper, libHelper)` — stores the helpers and subscribes to `monomerLib.onChanged`
- `monomerLibOnChangedHandler()` — logs the new library summary and shows a "Monomer lib updated" info toast

The legacy Dojo/JSDraw2 loader (`initHELMWebEditor`, `initHelmLoadAndPatchDojo`) and the `org.helm.webeditor.Monomers` patching (`initHelmPatchScilAlert`, `initHelmPatchPistoia`) are **gone** — hwe is self-contained and reads Bio's monomer library directly via `bridgeMonomerLib`.

### HELM Helper — `src/helm-helper.ts`

`HelmHelper` (implements `IHelmHelper`) — singleton. Every editor / parse / calculation path runs through hwe via a lazily-built `HelmHelperAdapter` over a `HelmService` (the `editorAdapter` getter). The service is constructed with `bridgeMonomerLib(monomerLib)`, the RDKit module (for monomer-structure rendering), and an `extraPanes` map providing the editor's package-contributed tabs:

| Pane | Source |
|---|---|
| **FASTA** / **BILN** | Read-only converted text via the seq handler's notation converters (BILN only for peptide chains; FASTA only for a single linear chain) |
| **Molecular Structure** | `Bio:toAtomicLevelPanel` widget |
| **Composition Analysis** | `Bio:compositionAnalysisWidget` widget |

Methods:
- `createHelmInput()` — creates the `HelmInput` widget
- `createHelmWebEditor(host, options)` — view-only editor instance (`IHelmWebEditor`) from the adapter
- `createWebEditorApp(host, helm)` — full editor app (palette + toolbar + canvas + built-in tabs + the extra panes above), wrapped in hwe's `LegacyAppWrapper` so it presents the legacy Pistoia `App` shape (`canvas.getHelm(true)`, `canvas.helm.jsd.m.atoms[i].selected`, …) — existing dialog consumers work unchanged
- `getMolfiles()`, `parse(helm, origin)`, `removeGaps(srcHelm)` — delegate to the adapter
- `getHoveredAtom(x, y, mol, height)` — delegates to `getHoveredMonomerFromEditorMol`
- `originalMonomersFuncs` / `revertOriginalMonomersFuncs()` / `overrideMonomersFuncs()` / `buildMonomersFuncsFromLib()` — thin pass-throughs to the adapter, kept only for the `IHelmHelper` contract (hwe reads `bridgeMonomerLib(...)` directly, so there is no global monomer registry to patch)
- private `_getSingleHelmDf()` / `_getSingleHelmDfSh()` — build a cached single-row HELM dataframe (tagged Macromolecule/helm; sets `.toAtomicWidgetWidth` to ~60% of window width) used by the Molecular-Structure / Composition / FASTA / BILN pane providers

## Utilities (`src/utils/`)

| File | Purpose |
|---|---|
| `helm-grid-cell-renderer.ts` | `HelmGridCellRenderer` + `HelmGridCellRendererBack`: **synchronous** grid cell renderer painting straight into the grid canvas via hwe `HelmService.renderToCanvas` (no async SVG→image / ImageData round-trip). Keeps one hwe service per monomer-library source (`serviceByLib`) and per-cell-value layout bookkeeping (`auxList`) for hover + highlight. Subscribes to monomer-lib changes (invalidate all) and to `MACROMOLECULE_HIGHLIGHT_EVENT_ID`; `drawHighlightOverlay` paints translucent rings around highlighted monomers. `onMouseMove` hit-tests the laid-out mol and shows the monomer tooltip + hover links |
| `get-hovered.ts` | `getHoveredMonomerFromEditorMol()`: nearest-atom hit testing (model coordinates) against a `HelmMol`. `getSeqMonomerFromHelmAtom()`: converts a `HelmAtom` to `ISeqMonomer` with canonical symbol (DNA/RNA sugar + phosphate aware) |
| `index.ts` | HELM string utilities: `parseHelm()` (extracts monomer symbols, bracket/paren aware), `removeGapsFromHelm()` (regex gap removal), `findMonomers()` (resolves via `monomerLib.getMonomer`), `getParts()`, and the `split()` / `detachAnnotation()` HELM section tokenizer |
| `err-info.ts` | `defaultErrorHandler()`: logs via the package logger + `grok.shell.error` |

## Widgets (`src/widgets/`)

| File | Purpose |
|---|---|
| `helm-input.ts` | `HelmInput` (extends `HelmInputBase`): interactive HELM input widget embedding a **view-only hwe editor** (`IHelmWebEditor`). Defaults to a compact 250×250 box; shows monomer tooltips on hover; "Click to edit" opens the full hwe editor app in a full-screen dialog and reads back the edited HELM + monomer selection (via the legacy `App` canvas shape). Re-fits the canvas on resize (`ui.onSizeChanged`). Supports `SeqValueBase`/string get/set and `molValue` |
| `properties-widget.ts` | `getPropertiesWidget()` / `getPropertiesDict()`: **async** property panel — molecular formula, weight, extinction coefficient computed from an off-screen hwe editor (obtained via the async `getHelmHelper()` factory, replacing the old off-screen JSDraw2 editor). Limited to sequences <1000 chars |

## Constants — `src/constants.ts`

- `jsonSdfMonomerLibDict` — mapping between SDF monomer fields and internal property names
- `SMILES`, `RGROUPS`, `MONOMER_SYMBOL`, `RGROUP_*`, `SDF_MONOMER_NAME` — field-name constants for monomer data
- `TAGS.cellRendererRenderError` — tag for tracking render errors

> HELM types (`HelmType`, `HelmMol`, `HelmAtom`, `App`, `Point`, `MonomersFuncs`, `IHelmWebEditor`, …) are imported directly from `@datagrok-libraries/bio/src/helm/types`. The old `src/types/` re-export module (`ScilModule` / `JSDraw2Module` / `OrgHelmModule`) was removed in the migration.

## Tests (`src/tests/`)

Test entry point: `src/package-test.ts` — imports all test files, exports `test()` and `initAutoTests()`.

| File | What it tests |
|---|---|
| `_first-tests.ts` | Package import sanity check |
| `helm-tests.ts` | HELM parsing, `parseHelm` utility |
| `findMonomers-tests.ts` | Missing monomer detection |
| `renderers-tests.ts` | Cell renderer behavior |
| `get-molfiles-tests.ts` | HELM → molfile conversion |
| `properties-widget-tests.ts` | Properties panel (formula, MW, ext. coeff.) |
| `parse-helm-tests.ts` | HELM string parsing edge cases |
| `helm-web-editor-tests.ts` | hwe editor / web-editor widget creation |
| `helm-input-tests.ts` | `HelmInput` widget behavior |
| `helm-helper-tests.ts` | `HelmHelper` methods (parse, removeGaps, getMolfiles) |
| `helm-substructure-filter.ts` | Substructure filtering for HELM |
| `helm-activity-cliffs.ts` | Activity cliffs computation for HELM |
| `to-atomic-level-ui-non-linear.ts` | Non-linear HELM → atomic level conversion |
| `utils.ts` | Shared test helpers (not a `test()` suite) |

## Initialization Flow

1. `detectors.js` `autostart` runs first — registers the `ui.input.helmAsync` polyfill
2. `initHelm()` on package init → loads the RDKit module, `SeqHelper`, `MonomerLibHelper`; creates `HelmHelper(seqHelper, logger, rdKitModule)`
3. `completeInit(helmHelper, libHelper)` → wires the helpers and subscribes to monomer-library changes
4. The hwe `HelmHelperAdapter` / `HelmService` is built **lazily** on first editor / parse / calculation use (after the library is available), reading Bio's monomer library via `bridgeMonomerLib`

## Key Dependencies

- `@datagrok-libraries/hwe` — the standalone Datagrok HELM Web Editor: HELM parsing, synchronous canvas rendering, molfile / property calculators, and the editor app (`HelmService`, `HelmHelperAdapter`, `bridgeMonomerLib`, `buildHelmMolView`, `PolymerTypes`, `styles.css`)
- `@datagrok-libraries/bio` — macromolecule + HELM types (`HelmType`, `HelmAtom`, `HelmMol`, `IHelmHelper`, `ISeqHelper`, `IMonomerLib`), cell-renderer base classes, monomer hover links, macromolecule-highlight events/consts
- `@datagrok-libraries/chem-meta` — RDKit API types (`RDModule`)
- `@datagrok-libraries/utils` — error / console helpers
- `cash-dom` — jQuery-like DOM manipulation
- `wu` — lazy iteration
- `rxjs` — reactive subscriptions

## Quick Lookups

| Looking for... | Check first |
|---|---|
| Function/panel/renderer registration | `src/package.ts` |
| Package singleton, init wiring | `src/package-utils.ts` (`HelmPackage`) |
| HELM helper + hwe adapter (parse, convert, getMolfiles, editor app, panes) | `src/helm-helper.ts` (`HelmHelper.editorAdapter`) |
| The hwe editor / service / adapter | `@datagrok-libraries/hwe` (external) |
| Grid cell renderer (synchronous canvas) | `src/utils/helm-grid-cell-renderer.ts` |
| Monomer highlight overlay (event-driven) | `src/utils/helm-grid-cell-renderer.ts` (`handleHighlightEvent` / `drawHighlightOverlay`) |
| Mouse hover / hit testing | `src/utils/get-hovered.ts` |
| HELM string parsing | `src/utils/index.ts` |
| HELM input widget | `src/widgets/helm-input.ts` |
| Properties panel (MW, formula) | `src/widgets/properties-widget.ts` |
| Constants (field names, tags) | `src/constants.ts` |
| HELM types (`HelmMol`, `HelmAtom`, `App`, …) | `@datagrok-libraries/bio/src/helm/types` |
| Auto-generated function wrappers | `src/package.g.ts` / `src/package-api.ts` |
| Autostart / input polyfill | `detectors.js` |
| Test entry point | `src/package-test.ts` |
| Sample data | `files/samples/`, `files/tests/` |
