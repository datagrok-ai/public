# KetcherSketcher

Registers EPAM's [Ketcher](https://lifescience.opensource.epam.com/ketcher/index.html) as a
molecule sketcher option in Datagrok. When the user picks "Ketcher" in any chem sketcher
dropdown (grid cell editor, dialog, filter, etc.), the platform instantiates the widget
defined here and uses it in place of the built-in OpenChemLib sketcher.

This package is small: one real source file (`src/ketcher.tsx`) wraps the third-party
`ketcher-react` editor and adapts it to the platform's `SketcherBase` contract.

## Architecture

| File | Role |
|------|------|
| `src/package.ts` | Registers the single function `Ketcher` with `meta.role: moleculeSketcher`. That's the whole public surface. |
| `src/ketcher.tsx` | `KetcherSketcher extends grok.chem.SketcherBase`. Mounts a React root holding `<Editor>` from `ketcher-react`, subscribes to `'change'` events, and exposes `smiles / molFile / molV3000 / smarts` getters and setters that bridge Ketcher's async API to the platform's synchronous contract. |
| `src/constants.ts` | String constants for molfile version args (`v2000` / `v3000`) passed to `ketcher.getMolfile(...)`. |
| `src/tests/ketcher-utils.ts` | Test helpers: `createKetcher()` opens a dialog containing a fresh `Sketcher`; `_testSetSmiles / _testSetMolfile / _testSetSmarts` verify each notation round-trips. |
| `css/editor.css` | Scoped overrides on Ketcher's generated class names — fixed min-width, hides About/Help, and a workaround rule for the macromolecules editor (see below). |

## Glossary

| Concept | What it maps to |
|---------|-----------------|
| Sketcher | `grok.chem.Sketcher` — the Datagrok host widget (dropdown + container). Delegates to an implementation like `KetcherSketcher` chosen by `grok.chem.currentSketcherType`. |
| Sketcher implementation | A `DG.Widget` that `extends grok.chem.SketcherBase`. Registered via a package function with `meta.role: moleculeSketcher`. The function name (here `Ketcher`) is what shows up in the sketcher-picker dropdown. |
| Notation | One of `smiles`, `molblock` (V2000), `molblockV3000`, `smarts`. `DG.chem.Notation` holds the enum; `DG.chem.convert(src, from, to)` converts between them. |
| `explicitMol` | A `{notation, value}` field on `SketcherBase` that short-circuits the getters: if the caller last set SMILES explicitly, `get smiles` returns that string verbatim instead of reconverting from a cached molblock. |
| Cached formats (`_smiles`, `_molV2000`, `_molV3000`, `_smarts`) | What Ketcher actually produced on the last `'change'` event (or what was pushed in via a setter). Getters for other notations fall through these via `DG.chem.convert`. Setters for one notation clear the others. |
| `host` / `setMoleculeFromHost` | The surrounding `Sketcher` container may already have a molecule set before this implementation initializes. `onInit` calls `setMoleculeFromHost()` to pick it up. |

## Conventions and gotchas

- **Editor instance not ready until `onInit` fires.** The React render is synchronous but `ketcher-react` calls `onInit(ketcher)` asynchronously. Before that, `this._sketcher` is `null` and `isInitialized` is `false`. Tests wait on `sketcher.isInitialized === true`.
- **`getSmiles()` throws on SMARTS.** Ketcher's SMILES export fails when the structure is a query. The `change` handler catches and sets `this._smiles = null`; the `smiles` getter then falls through to `DG.chem.smilesFromSmartsWarning()`.
- **Setters clear sibling caches.** `set smiles` nulls `_molV2000 / _molV3000 / _smarts`, etc. Don't reuse a cached format after calling a setter on a different one.
- **User-settings persistence is intentionally off.** The commented-out `grok.dapi.userDataStorage` calls in `onInit` and `detach` were removed in 2.1.8 because restoring saved options broke SMARTS input. Do not uncomment without re-verifying SMARTS.
- **Popup-truncation workaround (`ketcher.tsx:40-45`).** When the sketcher opens in a popup at the edge of the screen (last grid column), `.Ketcher-root` is forced to `min-width: 0` so it shrinks to fit. Keyed off `host.isInPopupContainer()`.
- **Macromolecules editor `display: none` hack (`css/editor.css:46-52`).** Ketcher mounts the macromolecules editor inside a `display: none` container; `RulerArea` reads `SVGLength.value` on mount and throws when the SVG is hidden. The CSS rule swaps `display: none` for `visibility: hidden` so SVG measurements still work. When the user switches to macro mode the inline style is removed and the rule stops matching.
- **CSS class names are hashed.** `ketcher-react` ships CSS modules, so selectors here use attribute-substring matchers like `[class*="App-module_top__"]`. If styles stop working after a Ketcher upgrade, first check whether the module name changed.
- **`file-loader` for `.sdf` in webpack.** `webpack.config.js` routes `.sdf` through `file-loader`; `templates/fg.sdf` and `templates/library.sdf` are shipped as static assets.
- **No direct `fetch` / `grok.dapi` use.** The sketcher is fully client-side via `StandaloneStructServiceProvider`. There is no server round-trip for molecule conversion.
