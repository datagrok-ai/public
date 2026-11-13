---
feature: helm
target_layer: apitest
coverage_type: regression
priority: p2
realizes_atlas: [helm.int.input-bio-monomer-lib]
realizes: []
realized_as:
  - helm-monomer-funcs-override-api-spec.ts
related_bugs: []
---

# Helm — monomer-funcs path (buildMonomersFuncsFromLib)

Exercises `buildMonomersFuncsFromLib`, which builds a monomer-lookup
implementation backed by the Datagrok monomer library (SequenceTranslator's
PolyTool enumeration uses it). Reflects the 2026 HELM editor rewrite: the editor
reads Bio's monomer library directly, so there is no global `Monomers`
dictionary to patch (that mechanism, `rewriteLibraries` and the
`overrideMonomersFuncs` / `revertOriginalMonomersFuncs` /
`originalMonomersFuncs` swap path are gone), and bracketed monomer symbols like
`[meI]` are no longer unwrapped by this layer; see grok-browser
`references/helm.md` § "Monomer-override contract (CHANGED)".

## Setup

1. Authenticate to Datagrok as the test user.
2. Helm package `@init` (`initHelm`) has completed: the new SVG HELM editor is
   loaded, the monomer library is wired through Bio, and the RDKit module is
   loaded from Chem.
3. Acquire the `IHelmHelper` singleton:
   `const hh = await grok.functions.call('Helm:getHelmHelper');`
4. Acquire a Datagrok `IMonomerLib` reference:
   `const lib = await grok.functions.call('Bio:getMonomerLibHelper')
   .then(h => h.getMonomerLib());`
   The library carries at least one peptide monomer (`A`, `meI`, `hHis`, …)
   when the default library is loaded.

## Scenarios

### Scenario 1: buildMonomersFuncsFromLib reads symbols from the Datagrok IMonomerLib

Steps:
1. Build the `MonomersFuncs` from the Datagrok lib:
   `const funcs = hh.buildMonomersFuncsFromLib(lib);`.
2. Assert `funcs != null`; assert it carries `getMonomer` and `getMonomerSet`
   as functions.
3. Plain-symbol lookup — pass a symbol the default lib resolves:
   `const m = funcs.getMonomer('PEPTIDE', 'A');`. Confirm `m != null`,
   `m.id === 'A'`, and `m` matches `lib.getWebEditorMonomer('PEPTIDE', 'A')`
   (same `id`; e.g. `n === 'Alanine'`).
4. Unknown-symbol lookup — pass a symbol the lib does NOT carry:
   `const mUnknown = funcs.getMonomer('PEPTIDE', 'Xz_NotInLib');`. Confirm
   it returns a `missing` placeholder: `{ id: 'Xz_NotInLib', n: 'missing' }`.
5. Bracketed-symbol lookup — pass a bracketed symbol:
   `const mBracketed = funcs.getMonomer('PEPTIDE', '[meI]');`.
   * **Behaviour change:** the outer brackets are **no longer stripped** at
     this layer — the lookup returns the `missing` placeholder
     (`{ id: '[meI]', n: 'missing' }`), NOT the unbracketed `meI` record.
     (Bracket handling now happens upstream in the notation parser.)
6. `getMonomerSet` smoke — call `funcs.getMonomerSet('PEPTIDE')` and confirm
   it does NOT throw.

Expected:
- Plain-symbol lookups resolve against the Datagrok lib's
  `getWebEditorMonomer` and return a record whose `id`/`n` match the lib's
  record for the same symbol.
- Unknown symbols resolve to a `missing` placeholder (`{id, n:'missing'}`);
  the wrapper does not throw.
- Bracketed symbols return the `missing` placeholder (brackets NOT stripped —
  changed behaviour).
- `getMonomerSet` is callable and does not throw.

### Removed scenarios: the override swap path and rewriteLibraries

The previous scenarios verified the `overrideMonomersFuncs` /
`revertOriginalMonomersFuncs` state machine (tracked via
`hh.originalMonomersFuncs`) and `rewriteLibraries(monomerLib)`, which mutated
the Pistoia dictionary via `org.helm.webeditor.Monomers.getMonomer /
addOneMonomer / clear`. In the new build none of these is part of `IHelmHelper`
and `org.helm.webeditor.Monomers` no longer exists, so there is **no surface**
left to test. The lib→editor sync is exercised by Scenario 1
(`buildMonomersFuncsFromLib`) and by the cross-feature interaction
`helm-input-bio-monomer-lib` (monomer-lib reload re-renders the renderer /
input / properties surfaces).
