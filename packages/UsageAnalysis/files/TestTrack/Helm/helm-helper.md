---
feature: helm
target_layer: apitest
coverage_type: regression
priority: p2
realizes: []
realized_as:
  - helm-helper-api-spec.ts
related_bugs: []
---

# Helm — HelmHelper accessor + factory + hit-test API surface

Exercises the `IHelmHelper` singleton returned by `Helm:getHelmHelper` — the
JS API surface plugins use to parse HELM sequences, build the HELM input and
editor widgets, and hit-test a point against a drawn structure. Confirms the
singleton and its full method set are present, that `createHelmInput` and
`createHelmWebEditor` build working widgets bound to the active monomer
library, and that `getHoveredAtom` returns the correct atom (or `null`) for a
given point. This targets the API after the 2026 rewrite that replaced the
JSDraw2/Dojo editor with a native SVG one — the public method surface is
unchanged, but `createHelmWebEditor().editor` is now a plain adapter object
rather than a JSDraw2 instance; see grok-browser `references/helm.md`.

## Setup

1. Authenticate to Datagrok as the test user.
2. Helm package init complete (`initHelm` `@init` has run): the new
   `@datagrok-libraries/helm-web-editor` SVG editor is loaded, the monomer
   library is wired through Bio, and the RDKit module is loaded from Chem.
   The functions exercised here all depend on a fully initialised
   `HelmHelper` singleton.
3. Datagrok monomer library: default platform monomer library loaded by Bio
   is sufficient. The widget factories pull the lib via
   `_package._libHelper!.getMonomerLib()` so the lib MUST be present or
   `createHelmInput` will throw.
4. Test fixture HELM strings (deterministic, atom-count-checkable):
   - peptide linear:
     `PEPTIDE1{meI.hHis.Aca.N.T.dE.Thr_PO3H2.Aca.D-Tyr_Et}$$$$`
     (9 monomers → 9 atoms, 8 bonds).
   - RNA short: `RNA1{r(A)p.r(C)p.r(G)p.r(U)p}$$$$` (→ 12 atoms, 11 bonds).

## Scenarios

### Scenario 1: Helm:getHelmHelper returns the singleton IHelmHelper

Steps:
1. Call the JS-API entry point:
   `const hh = await grok.functions.call('Helm:getHelmHelper');`.
2. Assert `hh != null` and that it is the same reference returned by a
   second call:
   `const hh2 = await grok.functions.call('Helm:getHelmHelper');
   assert(hh === hh2);`.
3. Inspect the returned object: it MUST carry the methods declared on the
   `IHelmHelper` interface — `parse`, `removeGaps`, `getMolfiles`,
   `getHoveredAtom`, `createHelmInput`, `createHelmWebEditor`,
   `createWebEditorApp`, `overrideMonomersFuncs`,
   `revertOriginalMonomersFuncs`, `buildMonomersFuncsFromLib`.

Expected:
- `Helm:getHelmHelper` resolves without throwing; the returned object is
  non-null.
- The singleton invariant holds — every call returns the same reference.
- Every method listed above is present as a function on the returned object.
  No method is `undefined`. (Verified live: all 10 present.)

### Scenario 2: HelmHelper.createHelmInput builds a HelmInput widget bound to the monomer lib

Steps:
1. With `hh` from Scenario 1, call
   `const input = hh.createHelmInput('macromolecule', {});`.
2. Assert `input != null`; assert it is an instance of `HelmInputBase` (the
   interface re-exported from
   `@datagrok-libraries/bio/src/helm/helm-helper`).
3. Read the widget's host element (`input.root` / `input.viewerHost` or
   equivalent `HelmInputBase` accessor) and confirm it is an `HTMLElement`
   that can be appended to a host div without throwing.
4. Negative path: call `createHelmInput` while the monomer library is
   intentionally unavailable. The wrapper in `helm-helper.ts` catches, logs
   via `errInfo` with prefix `Helm: HelmHelper.createHelmInput()`, then
   re-throws.

Expected:
- The widget construction returns a non-null `HelmInputBase` instance on the
  happy path; no console errors.
- The returned widget's root element is a real `HTMLElement`.
- On the negative path, the call rethrows; the error is surfaced via
  `logger.error`.

### Scenario 3: HelmHelper.createHelmWebEditor builds a view-only SVG editor wrapper

Steps:
1. Build a host div:
   `const host = ui.div([], {style: {width: '400px', height: '300px'}});`.
2. With `hh` from Scenario 1, call
   `const we = hh.createHelmWebEditor(host, {});`.
3. Assert `we != null`; assert it implements `IHelmWebEditor` (host div +
   `editor` accessor).
4. Read `we.host`; confirm it is the same host element passed in (the
   wrapper stores the host without re-parenting it).
5. Read `we.editor`; confirm it is a non-null **adapter object** exposing the
   editor methods (`setHelm`, `getHelm`, `getMolfile`, `getFormula`,
   `getMolWeight`, `resize`, `redraw`, …). **Do NOT assert a JSDraw2
   instance or a `viewonly` flag** — the new editor is not JSDraw2 and does
   not expose `editor.options.viewonly` (it is null).
6. Call the wrapper a second time without a host argument:
   `const we2 = hh.createHelmWebEditor(undefined, {});`. Confirm the no-host
   construction also returns a non-null wrapper whose `host` is a
   freshly-created `HTMLElement`.

Expected:
- `createHelmWebEditor` returns a non-null `IHelmWebEditor` instance.
- `we.host` matches the host div passed in by reference; no clone.
- `we.editor` is a non-null adapter object carrying the editor methods
  listed above.
- The no-host call does not throw; the returned wrapper's `host` is an
  `HTMLElement`.
- No console errors on either call.

### Scenario 4: HelmHelper.getHoveredAtom hit-tests a point against a HelmMol and returns the nearest atom

Steps:
1. With `hh` from Scenario 1, parse a deterministic fixture to a `HelmMol`:
   `const mol = hh.parse('PEPTIDE1{meI.hHis.Aca.N.T.dE.Thr_PO3H2.Aca.D-Tyr_Et}$$$$');`.
2. Assert `mol != null`; record `mol.atoms.length` (expected 9 for the
   linear 9-monomer chain).
3. Pick atom 0 and read its coordinate:
   `const a0 = mol.atoms[0]; const p0 = a0.p;`.
4. Call `hh.getHoveredAtom(p0.x, p0.y, mol, 100)` and confirm the returned
   atom is `a0` (or another atom co-located within the hit-radius of `p0`).
5. Pick a coordinate guaranteed to be far from any atom (e.g.
   `(-10000, -10000)`); call `hh.getHoveredAtom(-10000, -10000, mol, 100)`
   and confirm the result is `null`.
6. Vary the `height` argument across realistic cell heights (40, 80, 200) at
   a known atom coordinate; confirm the method does not throw on any of them.

Expected:
- The happy-path call returns a non-null `HelmAtom` reference; for a query at
  the exact atom coordinate, the returned atom IS `mol.atoms[0]` (or
  co-located).
- The far-away call returns `null`; the method does not invent a match.
- The method runs without throwing across the height range; it is DOM-free
  (delegates to `getHoveredMonomerFromEditorMol`).
- The fixture's atom count (9) is invariant across the three calls; the
  method does not mutate `mol.atoms`.
