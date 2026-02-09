# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

**Helm** (`@datagrok/helm`) is a Datagrok plugin providing full support for **HELM notation** (Hierarchical Editing Language for Macromolecules) — importing, detecting, rendering, editing, and conversion. It wraps the Pistoia HELM Web Editor and Scilligence JSDraw2 libraries, integrating them into the Datagrok platform with custom cell renderers, input widgets, and monomer library management.

Category: **Bioinformatics**.

## Build Commands

```bash
npm install
npm run build              # grok api && grok check --soft && webpack
npm run test               # grok test
npm run lint               # eslint src --ext .ts
npm run lint-fix
npm run build-all          # Builds chem-meta → js-api → utils → bio → this package
```

## Architecture

### Entry Point — `src/package.ts`

`PackageFunctions` class registers all platform-visible functions via `@grok.decorators`:
- `initHelm` — `@init`: loads RdKit, SeqHelper, MonomerLibHelper, initializes HELMWebEditor, creates `HelmHelper` singleton, patches Pistoia monomer system
- `getHelmService` — `@func`: returns `HelmServiceBase` singleton for off-screen HELM rendering
- `helmCellRenderer` — `@func` (role: `cellRenderer`): returns `HelmGridCellRenderer` for columns with `quality=Macromolecule, units=helm`
- `editMoleculeCell` — `@func` (role: `cellEditor`): opens full-screen HELM Web Editor dialog for cell editing
- `openEditor` (`Edit Helm...`) — `@func` (action): context-menu editor for Macromolecule columns, converts non-HELM notations to HELM before editing
- `propertiesWidget` (`Properties`) — `@panel`: shows molecular formula, weight, extinction coefficient for Macromolecule values
- `getMolfiles` — `@func`: converts a Macromolecule column to molfile strings
- `helmInput` — `@func` (role: `valueEditor`): creates `HelmInput` widget for Macromolecule/HELM semtype
- `getHelmHelper` — `@func`: returns `IHelmHelper` singleton
- `measureCellRenderer` — `@func`: performance measurement utility for cell rendering
- `highlightMonomers` — `@func`: test app that loads peptide data and calls `Bio:toAtomicLevel` with highlight

**Exports**: `_package` (instance of `HelmPackage`), auto-generated function wrappers from `package.g.ts`

### Detector — `detectors.js`

`HelmPackageDetectors` class:
- `autostart` — registers `ui.input.helmAsync` polyfill so other packages can create HELM inputs before this package fully loads

### Package Utilities — `src/package-utils.ts`

`HelmPackage` class (extends `DG.Package`) — the package singleton managing:
- `helmHelper` / `seqHelper` — accessors for initialized helpers
- `_libHelper` — monomer library helper
- `monomerLib` — accessor for current `IMonomerLib`
- `initHELMWebEditor()` — loads Dojo toolkit (bundled), JSDraw2 + HELM Web Editor, patches `dojox.gfx.svg.Text.getTextWidth` for hang prevention
- `initHelmPatchScilAlert()` — suppresses Scilligence alert popups for missing monomers
- `initHelmPatchPistoia()` — overrides `org.helm.webeditor.Monomers.getMonomer` to use Datagrok's monomer library instead of built-in dictionaries
- `completeInit()` — ties helmHelper + libHelper together, rewrites Pistoia monomer dictionaries, subscribes to monomer lib changes
- `monomerLibOnChangedHandler()` — resyncs Pistoia dictionaries when monomer library updates

`initHelmLoadAndPatchDojo()` — complex initialization:
1. Configures and loads bundled Dojo 1.10.10 with custom `loaderPatch.injectUrl`
2. Requires Dojo modules (dojo/window, dojox/gfx, dijit/*, etc.)
3. Waits for all modules to be ready
4. Patches `dojox.gfx.svg.Text.prototype.getTextWidth` to prevent infinite loops

### HELM Helper — `src/helm-helper.ts`

`HelmHelper` (implements `IHelmHelper`) — singleton providing:
- `createHelmInput()` — creates `HelmInput` widget
- `createHelmWebEditor()` — creates `HelmWebEditor` viewer instance
- `createWebEditorApp()` — creates full HELM Web Editor app (with monomer explorer, tabs, sequence/notation views) for dialog-based editing. Includes custom "Placeholders" tab
- `overrideMonomersFuncs()` / `revertOriginalMonomersFuncs()` — swap Pistoia's `getMonomer`/`getMonomerSet` with custom implementations
- `buildMonomersFuncsFromLib()` — builds `getMonomer` function that reads from Datagrok's `IMonomerLib`
- `getMolfiles()` — batch converts HELM strings to molfile format using a cached off-screen JSDraw2 editor
- `parse()` — parses HELM string into `HelmMol` object model
- `removeGaps()` — removes gap monomers from HELM, re-links bonds, returns mapping of old→new monomer positions
- `getHoveredAtom()` — hit-tests a point against a `HelmMol` to find the nearest atom

### HELM Web Editor Wrapper — `src/helm-web-editor.ts`

`HelmWebEditor` (implements `IHelmWebEditor`) — lightweight wrapper around `JSDraw2.Editor`:
- Creates a viewonly editor in a host div
- Handles resize events to keep editor synced with container size
- Exposes `editor` (the JSDraw2 `HelmEditor` instance) and `host` div

## Utilities (`src/utils/`)

| File | Purpose |
|---|---|
| `helm-grid-cell-renderer.ts` | `HelmGridCellRenderer` + `HelmGridCellRendererBack`: async cell renderer using SVG→image pipeline. Renders HELM structures in grid cells with LRU caching, handles mouse hover for monomer tooltips and hover links |
| `helm-service.ts` | `HelmService` (extends `HelmServiceBase`): off-screen rendering service. Uses hidden JSDraw2 editors (LRU-cached per monomer lib) to render HELM→SVG→ImageData. Handles scaling/fitting within cell bounds |
| `get-monomer.ts` | `rewriteLibraries()`: syncs Datagrok monomer library → Pistoia `org.helm.webeditor.Monomers` dictionary. `getMonomerOverrideAndLogAlert()`: temporary getMonomer override with alert suppression |
| `get-hovered.ts` | `getHoveredMonomerFromEditorMol()`: nearest-atom hit testing in editor molecule coordinates. `getSeqMonomerFromHelmAtom()`: converts `HelmAtom` to `ISeqMonomer` with canonical symbol |
| `index.ts` | HELM string parsing utilities: `parseHelm()` (extracts monomer symbols), `removeGapsFromHelm()` (regex gap removal), `findMonomers()` (finds missing monomers), `split()` / `detachAnnotation()` (HELM section tokenizer) |
| `err-info.ts` | `defaultErrorHandler()`: logs errors via package logger + `grok.shell.error` |

## Widgets (`src/widgets/`)

| File | Purpose |
|---|---|
| `helm-input.ts` | `HelmInput` (extends `HelmInputBase`): interactive HELM sequence input widget. Embeds a viewonly `HelmWebEditor`, shows monomer tooltips on hover, "Click to edit" hint opens full editor dialog. Supports `SeqValueBase` get/set, string get/set, `molValue` access. Handles mouse events for hover/click/tooltip |
| `properties-widget.ts` | `getPropertiesWidget()`: property panel showing molecular formula, weight, and extinction coefficient. Uses off-screen JSDraw2 editor to compute properties. Limits to sequences <1000 chars |

## Constants — `src/constants.ts`

- `jsonSdfMonomerLibDict` — mapping between SDF monomer fields and internal property names
- `SMILES`, `RGROUPS`, `MONOMER_SYMBOL`, `RGROUP_*` — field name constants for monomer data
- `SDF_MONOMER_NAME` — SDF field name
- `TAGS.cellRendererRenderError` — tag for tracking render errors

## Types — `src/types/index.ts`

Re-exports from `@datagrok-libraries/bio/src/helm/types`:
- `ScilModule` — Scilligence utility module (`scil.*`)
- `JSDraw2Module` — JSDraw2 editor module (`JSDraw2.*`)
- `OrgHelmModule` — Pistoia HELM module (`org.helm.*`)

## Vendor Libraries

### `vendor/dojo-1.10.10/`
Bundled Dojo Toolkit (dijit, dojo, dojox) required by JSDraw2 and HELM Web Editor. Loaded via custom `dojoConfig.loaderPatch.injectUrl` to work within webpack bundle.

### `helm/JSDraw/`
Legacy uncompressed sources of Scilligence JSDraw2 Lite and Pistoia HELM Web Editor. Currently replaced by `@datagrok-libraries/helm-web-editor` npm package.

## Tests (`src/tests/`)

Test entry point: `src/package-test.ts` — imports all test files, exports `test()` and `initAutoTests()`.

| File | What it tests |
|---|---|
| `_first-tests.ts` | Package import sanity check |
| `helm-tests.ts` | HELM parsing, parseHelm utility |
| `findMonomers-tests.ts` | Missing monomer detection |
| `helm-service-tests.ts` | HelmService rendering |
| `renderers-tests.ts` | Cell renderer behavior |
| `get-molfiles-tests.ts` | HELM → molfile conversion |
| `properties-widget-tests.ts` | Properties panel (formula, MW, ext. coeff.) |
| `get-monomer-tests.ts` | Monomer lookup and library integration |
| `parse-helm-tests.ts` | HELM string parsing edge cases |
| `helm-web-editor-tests.ts` | Web editor widget creation |
| `helm-input-tests.ts` | HelmInput widget behavior |
| `helm-helper-tests.ts` | HelmHelper methods (parse, removeGaps, getMolfiles) |
| `helm-substructure-filter.ts` | Substructure filtering for HELM |
| `helm-activity-cliffs.ts` | Activity cliffs computation for HELM |
| `to-atomic-level-ui-non-linear.ts` | Non-linear HELM → atomic level conversion |

## Initialization Flow

1. `detectors.js` `autostart` runs first — registers `ui.input.helmAsync` polyfill
2. `initHelm()` called on package init → loads RdKit, SeqHelper, MonomerLibHelper
3. `initHELMWebEditor()` → loads bundled Dojo, waits for all Dojo modules, patches SVG text width
4. Loads `@datagrok-libraries/helm-web-editor` bundle, waits for `helmWebEditor$.initPromise`
5. Creates `HelmHelper` singleton
6. `completeInit()` → rewrites Pistoia monomer dictionaries from Datagrok lib, patches `getMonomer`, subscribes to lib changes

## Key Dependencies

- `@datagrok-libraries/bio` — macromolecule types, HELM types (`HelmType`, `HelmAtom`, `HelmMol`, `IHelmHelper`, `ISeqHelper`, `IMonomerLib`), cell renderer base classes, monomer hover links
- `@datagrok-libraries/helm-web-editor` — bundled JSDraw2 + Pistoia HELM Web Editor
- `@datagrok-libraries/chem-meta` — chemistry metadata
- `@datagrok-libraries/utils` — SVG utilities (`svgToImage`), console helpers
- `cash-dom` — jQuery-like DOM manipulation
- `lru-cache` — LRU caching for editors and rendered images
- `wu` — lazy iteration
- `rxjs` — reactive subscriptions

## Quick Lookups

| Looking for... | Check first |
|---|---|
| Function/panel/renderer registration | `src/package.ts` |
| Package singleton, Dojo/HWE loading | `src/package-utils.ts` (`HelmPackage`) |
| HELM helper (parse, convert, getMolfiles) | `src/helm-helper.ts` |
| Web editor wrapper | `src/helm-web-editor.ts` |
| Grid cell renderer (async SVG→image) | `src/utils/helm-grid-cell-renderer.ts` |
| Off-screen render service | `src/utils/helm-service.ts` |
| Monomer library sync with Pistoia | `src/utils/get-monomer.ts` |
| Mouse hover / hit testing | `src/utils/get-hovered.ts` |
| HELM string parsing | `src/utils/index.ts` |
| HELM input widget | `src/widgets/helm-input.ts` |
| Properties panel (MW, formula) | `src/widgets/properties-widget.ts` |
| Constants (field names, tags) | `src/constants.ts` |
| Type re-exports (Scil, JSDraw2, Org) | `src/types/index.ts` |
| Auto-generated function wrappers | `src/package.g.ts` / `src/package-api.ts` |
| Autostart / input polyfill | `detectors.js` |
| Test entry point | `src/package-test.ts` |
| Sample data | `files/samples/`, `files/tests/` |
