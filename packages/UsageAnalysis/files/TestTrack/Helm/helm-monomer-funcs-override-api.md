---
feature: helm
target_layer: apitest
coverage_type: regression
priority: p2
realizes_atlas: [helm-input-bio-monomer-lib]
realizes: []
realized_as:
  - helm-monomer-funcs-override-api-spec.ts
related_bugs: []
---

# Helm — monomer-funcs swap path (overrideMonomersFuncs + buildMonomersFuncsFromLib)

Exercises the API that lets a plugin temporarily swap out how the HELM editor
looks up monomers: `overrideMonomersFuncs` installs a substitute
monomer-lookup implementation and `revertOriginalMonomersFuncs` restores the
original, while `buildMonomersFuncsFromLib` builds a monomer-lookup
implementation backed by the Datagrok monomer library. Confirms the
override/revert state machine (tracked via `hh.originalMonomersFuncs`) behaves
correctly, including that double-applying an override or reverting with
nothing to revert no longer throw. Reflects the 2026 HELM editor rewrite: the
override no longer touches a global `Monomers` dictionary (that mechanism,
along with `rewriteLibraries`, is gone), and bracketed monomer symbols like
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
5. Baseline hygiene: if `hh.originalMonomersFuncs != null` at start, call
   `hh.revertOriginalMonomersFuncs()` once so the swap state is clean.

## Scenarios

### Scenario 1: overrideMonomersFuncs swaps the adapter funcs and revert restores them (state-machine contract)

Steps:
1. Confirm the clean baseline: `hh.originalMonomersFuncs === null`.
2. Build a sentinel `MonomersFuncs` whose handlers are recognisable stubs:
   `const sentinel = { getMonomer: (a, name) => ({id: 'S1'}), getMonomerSet: (a) => null };`.
3. Call `const previous = hh.overrideMonomersFuncs(sentinel);`. Confirm
   `previous != null` and that `hh.originalMonomersFuncs` is now **non-null**
   (the override is recorded).
4. Double-apply: call `hh.overrideMonomersFuncs(anotherSentinel)` while the
   first override is still in effect.
   * **Behaviour change:** this **no longer throws** (the old
     `originalGetMonomer is overridden already` guard was removed). Record
     that it returns without throwing.
5. Restore: `const outgoing = hh.revertOriginalMonomersFuncs();`. Confirm
   `outgoing != null` and that `hh.originalMonomersFuncs` is back to `null`.
6. Pre-revert (no override in effect): call `hh.revertOriginalMonomersFuncs()`
   again.
   * **Behaviour change:** this **no longer throws** (the old `Unable to
     revert original getMonomer` guard was removed). Record that it returns
     without throwing.

Expected:
- `overrideMonomersFuncs` returns the previous funcs (non-null) and flips
  `hh.originalMonomersFuncs` from `null` to non-null.
- Double-apply does NOT throw (changed behaviour).
- `revertOriginalMonomersFuncs` returns the outgoing funcs (non-null) and
  flips `hh.originalMonomersFuncs` back to `null`.
- Pre-revert does NOT throw (changed behaviour).
- No error balloon; the teardown leaves `hh.originalMonomersFuncs === null`
  so subsequent scenarios run cleanly.

> NOTE: there is no `window.org.helm.webeditor.Monomers` to assert against in
> the new build. The override's effect on monomer lookups is internal to the
> `editorAdapter`; the durable, observable contract is the
> `originalMonomersFuncs` null↔non-null state machine plus the non-null
> return values above.

### Scenario 2: buildMonomersFuncsFromLib reads symbols from the Datagrok IMonomerLib

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

### Scenario 3 (REMOVED in the new build): rewriteLibraries

The previous Scenario 3 verified `rewriteLibraries(monomerLib)` by reading and
mutating the Pistoia dictionary via
`org.helm.webeditor.Monomers.getMonomer / addOneMonomer / clear`. In the new
build:

- `rewriteLibraries` is not a method on `IHelmHelper` and the bundled module is
  not import()-able from the apitest layer.
- `org.helm.webeditor.Monomers` (and `addOneMonomer` / `clear`) no longer exist.

There is therefore **no automatable surface** for a direct rewrite-libraries
test. The lib→editor sync is exercised indirectly by Scenario 2
(`buildMonomersFuncsFromLib`) and by the cross-feature interaction
`helm-input-bio-monomer-lib` (monomer-lib reload re-renders the renderer /
input / properties surfaces). This scenario is intentionally a no-op placeholder
documenting the removal.
