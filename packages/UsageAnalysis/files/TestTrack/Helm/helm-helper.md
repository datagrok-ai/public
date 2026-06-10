# Helm — HelmHelper accessor + factory + hit-test API surface

## Setup

1. Authenticate to Datagrok as the test user.
2. Helm package init complete (`initHelm` `@init` has run): the
   bundled Dojo 1.10.10 + JSDraw2 + Pistoia HELM Web Editor stack is
   loaded, the monomer dictionary has been rewritten via
   `_package.completeInit`, and the RDKit module is loaded from
   Chem. The four functions exercised here all depend on a fully
   initialised `HelmHelper` singleton.
3. Datagrok monomer library: default platform monomer library loaded
   by Bio is sufficient. The widget factories pull the lib via
   `_package._libHelper!.getMonomerLib()` so the lib MUST be present
   or `createHelmInput` will throw.
4. Test fixture HELM strings (deterministic, atom-count-checkable):
   - peptide linear:
     `PEPTIDE1{meI.hHis.Aca.N.T.dE.Thr_PO3H2.Aca.D-Tyr_Et}$$$$`
     (9 monomers, matches `parseHelm` test fixture in atlas
     `source_classes[0].examples`).
   - RNA short:
     `RNA1{r(A)p.r(C)p.r(G)p.r(U)p}$$$$`
     (matches `renderers-tests scatterPlotTooltip` fixture in
     atlas `source_classes[0].examples`).

## Scenarios

### Scenario 1: Helm:getHelmHelper returns the singleton IHelmHelper

Steps:
1. Call the JS-API entry point:
   `const hh = await grok.functions.call('Helm:getHelmHelper');`.
2. Assert `hh != null` and that it is the same reference returned by a
   second call:
   `const hh2 = await grok.functions.call('Helm:getHelmHelper');
   assert(hh === hh2);`.
3. Inspect the returned object: it MUST carry the methods declared on
   the `IHelmHelper` interface — `parse`, `removeGaps`, `getMolfiles`,
   `getHoveredAtom`, `createHelmInput`, `createHelmWebEditor`,
   `createWebEditorApp`, `overrideMonomersFuncs`,
   `revertOriginalMonomersFuncs`, `buildMonomersFuncsFromLib`.

Expected:
- `Helm:getHelmHelper` resolves without throwing; the returned object
  is non-null.
- The singleton invariant holds — every call returns the same
  reference (the underlying `HelmHelper` constructor enforces
  `instanceCount > 1 -> throw`, so two distinct instances would have
  surfaced as a startup error long before this scenario runs).
- Every method listed above is present as a function on the returned
  object. No method is `undefined`.

### Scenario 2: HelmHelper.createHelmInput builds a HelmInput widget bound to the monomer lib

Steps:
1. With `hh` from Scenario 1, call
   `const input = hh.createHelmInput('macromolecule', {});`.
2. Assert `input != null`; assert it is an instance of
   `HelmInputBase` (the interface re-exported from
   `@datagrok-libraries/bio/src/helm/helm-helper`).
3. Read the widget's host element (`input.root` or equivalent
   `HelmInputBase` accessor) and confirm it is an `HTMLElement` that
   can be appended to a host div without throwing.
4. Negative path: call `createHelmInput` while the monomer library is
   intentionally unavailable. Simulate by stubbing
   `_package._libHelper.getMonomerLib` to throw, invoke
   `hh.createHelmInput()`, and confirm the wrapper rethrows the
   original error AFTER logging via the package logger.

Expected:
- The widget construction returns a non-null `HelmInputBase` instance
  on the happy path; no console errors.
- The returned widget's root element is a real `HTMLElement` (the
  widget does its own DOM building inside the constructor).
- On the negative path, the call rethrows; the error message is
  surfaced via `logger.error` with prefix
  `Helm: HelmHelper.createHelmInput()` (the wrapper in
  `helm-helper.ts#L54` catches, logs with `errInfo`, then re-throws).
- The original stubbed `getMonomerLib` is restored in the test
  teardown so subsequent scenarios run cleanly.

### Scenario 3: HelmHelper.createHelmWebEditor builds a view-only JSDraw2 editor

Steps:
1. Build a host div:
   `const host = ui.div([], {style: {width: '400px', height: '300px'}});`.
2. With `hh` from Scenario 1, call
   `const we = hh.createHelmWebEditor(host, {});`.
3. Assert `we != null`; assert it implements `IHelmWebEditor` (host
   div + `editor` accessor).
4. Read `we.host`; confirm it is the same host element passed in
   (the wrapper stores the host without re-parenting it).
5. Read `we.editor`; confirm it is a non-null JSDraw2 editor instance
   (the wrapper instantiates `new JSDraw2.Editor(host, { viewonly:
   true })`).
6. Call the wrapper a second time without a host argument:
   `const we2 = hh.createHelmWebEditor(undefined, {});`. The
   `HelmWebEditor` wrapper accepts `host?: HTMLDivElement`; confirm
   the no-host construction also returns a non-null wrapper (the
   internal editor either creates its own host or stays detached
   until one is supplied).

Expected:
- `createHelmWebEditor` returns a non-null `IHelmWebEditor`
  instance.
- `we.host` matches the host div passed in by reference; no clone.
- `we.editor` is a JSDraw2 editor instance whose `viewonly` flag is
  `true` (the wrapper passes `{viewonly: true}` to the JSDraw2
  constructor at `helm-web-editor.ts#L20`).
- The no-host call does not throw; the returned wrapper is still a
  well-shaped `IHelmWebEditor` (host accessor may be null or a
  freshly-created div, depending on the wrapper's internal default).
- No console errors on either call.

### Scenario 4: HelmHelper.getHoveredAtom hit-tests a point against a HelmMol and returns the nearest atom

Steps:
1. With `hh` from Scenario 1, parse a deterministic fixture to a
   `HelmMol`:
   `const mol = hh.parse('PEPTIDE1{meI.hHis.Aca.N.T.dE.Thr_PO3H2.Aca.D-Tyr_Et}$$$$');`.
2. Assert `mol != null`; record `mol.atoms.length` (expected 9 for
   the linear 9-monomer chain).
3. Pick atom 0 (the first monomer) and read its coordinate:
   `const a0 = mol.atoms[0]; const p0 = a0.p;`.
4. Call `hh.getHoveredAtom(p0.x, p0.y, mol, 100)` and confirm the
   returned atom is `a0` (or another atom co-located within the
   hit-radius of `p0`).
5. Pick a coordinate guaranteed to be far from any atom (e.g.
   `(-10000, -10000)`); call `hh.getHoveredAtom(-10000, -10000, mol,
   100)` and confirm the result is `null` (no atom within the hit
   radius).
6. Vary the `height` argument across realistic cell heights (40, 80,
   200) at a known atom coordinate; confirm the method does not
   throw on any of them.

Expected:
- The happy-path call returns a non-null `HelmAtom` reference; for a
  query at the exact atom coordinate, the returned atom IS
  `mol.atoms[0]` (or co-located).
- The far-away call returns `null`; the method does not invent a
  match.
- The method runs without throwing across the height range; it is
  DOM-free (delegates to `getHoveredMonomerFromEditorMol` from
  `utils/get-hovered.ts`) so it does not require a mounted JSDraw2
  editor.
- The fixture's atom count (9) is invariant across the three calls;
  the method does not mutate `mol.atoms`.