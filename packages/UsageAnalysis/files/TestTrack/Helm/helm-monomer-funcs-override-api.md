# Helm — Pistoia monomer-dict swap path (overrideMonomersFuncs + buildMonomersFuncsFromLib + rewriteLibraries)

## Setup

1. Authenticate to Datagrok as the test user.
2. Helm package `@init` (`initHelm`) has completed: bundled
   Dojo 1.10.10 + JSDraw2 + Pistoia HELM Web Editor are loaded,
   `_package.completeInit` has rewritten the Pistoia monomer
   dictionary from the platform monomer library, and the RDKit
   module is loaded from Chem. After `completeInit` the swap
   path's "original" reference state is the Datagrok-rewritten
   dictionary, NOT the unmodified Pistoia defaults.
3. Acquire the `IHelmHelper` singleton:
   `const hh = await grok.functions.call('Helm:getHelmHelper');`
   The singleton invariant (every call returns the same reference)
   is exercised by sibling scenario `helm-api-helm-helper.md`; this
   scenario consumes the singleton as a precondition.
4. Acquire a Datagrok `IMonomerLib` reference:
   `const lib = grok.functions.call('Bio:getMonomerLibHelper').then(h => h.getMonomerLib())`
   resolves to the active platform monomer library (the same one
   `completeInit` already loaded into Pistoia's dictionary). The
   library carries at least one peptide monomer (`A`, `meI`,
   `hHis`, ...) and at least one RNA monomer (`r(A)`, `r(C)`, ...)
   when the default library is loaded.
5. Test fixture HELM strings for the symbol-lookup assertions:
   - `PEPTIDE1{A.G.S}$$$$` — three canonical peptide monomers
     present in the default lib.
   - `PEPTIDE1{[meI].[hHis]}$$$$` — bracketed multi-char monomer
     symbols (the bracket-stripping path of
     `buildMonomersFuncsFromLib`).

## Scenarios

### Scenario 1: overrideMonomersFuncs swaps the Pistoia dictionary and revert restores it; double-apply and pre-revert throw

Steps:
1. Capture the pre-swap dictionary references:
   `const preMonomers = (window as any).org.helm.webeditor.Monomers;`
   `const preGetMonomer = preMonomers.getMonomer;`
   `const preGetMonomerSet = preMonomers.getMonomerSet;`.
2. Build a sentinel `MonomersFuncs` whose handlers are recognisable
   stubs (return distinct sentinel values for any input). This
   isolates the swap from the production monomer-lib lookup.
3. Call `const original = hh.overrideMonomersFuncs(sentinelFuncs);`.
   Confirm `original` is non-null and carries `getMonomer` /
   `getMonomerSet` properties whose values are
   `=== preGetMonomer` / `=== preGetMonomerSet` (i.e. the pre-
   swap dictionary handlers, captured for restoration).
4. Confirm the live dictionary has been replaced:
   `(window as any).org.helm.webeditor.Monomers.getMonomer === sentinelFuncs.getMonomer`
   and the same for `getMonomerSet`.
5. Negative path — double-apply: call
   `hh.overrideMonomersFuncs(anotherSentinelFuncs)` while the
   first override is still in effect. Confirm the call throws
   with message containing `originalGetMonomer is overridden
   already` (per `helm-helper.ts#L249`). Confirm the dictionary
   is unchanged — `getMonomer` still references the first
   sentinel.
6. Restore: `const overridden = hh.revertOriginalMonomersFuncs();`.
   Confirm `overridden` is non-null and its `getMonomer` /
   `getMonomerSet` are the first sentinel handlers (the return is
   the SECOND override candidate of the swap — what was just
   removed, not what was restored).
7. Confirm the live dictionary is restored:
   `(window as any).org.helm.webeditor.Monomers.getMonomer === preGetMonomer`
   and the same for `getMonomerSet`.
8. Negative path — pre-revert: call
   `hh.revertOriginalMonomersFuncs()` again. Confirm the call
   throws with message containing `Unable to revert original
   getMonomer` (per `helm-helper.ts#L230`).

Expected:
- `overrideMonomersFuncs` returns the previous dictionary handlers
  by reference (not a clone); the live
  `org.helm.webeditor.Monomers` slots are replaced by the supplied
  `monomersFuncs.getMonomer` / `getMonomerSet`.
- Double-apply throws synchronously; the dictionary is NOT
  mutated by the failed call.
- `revertOriginalMonomersFuncs` restores the live dictionary slots
  to the references captured at override time and returns the
  outgoing (sentinel) handlers.
- Calling revert with no override in effect throws synchronously.
- The teardown (revert called once) leaves the platform dictionary
  in the same state subsequent scenarios expect — the renderer /
  helm input / properties panel will keep looking up monomers
  against the production lib.

### Scenario 2: buildMonomersFuncsFromLib reads symbols from the Datagrok IMonomerLib and strips outer square brackets

Steps:
1. Build the `MonomersFuncs` from the Datagrok lib:
   `const funcs = hh.buildMonomersFuncsFromLib(lib);`.
2. Assert `funcs != null`; assert it carries `getMonomer` and
   `getMonomerSet` as functions.
3. Plain-symbol lookup — pass a symbol that the default Datagrok
   lib resolves (e.g. peptide `A`):
   `const m = funcs.getMonomer(null, 'A');`. Confirm `m != null`
   and that the returned object matches the shape
   `lib.getWebEditorMonomer('PEPTIDE', 'A')` returns directly
   (same `id` / `symbol` / molfile-like fields, depending on the
   `IMonomerLib` shape — assert at least one stable identifier
   field is present and equal).
4. Bracketed-symbol lookup — pass a bracketed symbol:
   `const mBracketed = funcs.getMonomer(null, '[meI]');`. Confirm
   the call returns the same monomer record as
   `funcs.getMonomer(null, 'meI')` (i.e. the outer brackets are
   stripped before the lib lookup, per `helm-helper.ts#L269`).
5. Unknown-symbol lookup — pass a symbol the lib does NOT carry:
   `const mUnknown = funcs.getMonomer(null, 'Xz_NotInLib');`.
   Confirm the result is null / undefined (the lib's
   `getWebEditorMonomer` returns null for unknown symbols and
   the wrapper passes that through unchanged).
6. `getMonomerSet` fall-through — call
   `funcs.getMonomerSet(null, 'PEPTIDE')` and confirm it does
   NOT throw. The atlas notes `getMonomerSet` falls through to
   the original Pistoia implementation; the test asserts the
   call shape works regardless of the production lib's
   contents.

Expected:
- Plain-symbol lookups resolve against the Datagrok lib's
  `getWebEditorMonomer` and return a monomer record whose stable
  identifier (id/symbol) matches the lib's record for the same
  symbol.
- Bracketed symbols (e.g. `[meI]`) resolve identically to their
  unbracketed equivalents (`meI`) — the wrapper strips a single
  outer pair of `[]` before lookup.
- Unknown symbols resolve to null / undefined; the wrapper does
  not throw.
- `getMonomerSet` is callable and does not throw on a polymer-
  type argument the lib doesn't enumerate (it delegates to the
  original Pistoia implementation, which is the production
  fall-through path).

### Scenario 3: rewriteLibraries syncs the Datagrok IMonomerLib into the Pistoia dictionary

Steps:
1. Import the `rewriteLibraries` export from the package:
   `const { rewriteLibraries } = await import('@datagrok/helm/utils/get-monomer');`
   (or invoke via the package's exposed surface; the atlas
   `source:` for `helm.utils.rewrite-libraries` is
   `public/packages/Helm/src/utils/get-monomer.ts#L1` and the
   function is public per `CLAUDE.md#Package Utilities`).
2. Capture a "before" snapshot of the Pistoia dictionary
   contents: for each fixture monomer symbol (`A`, `meI`,
   `r(A)`), read
   `(window as any).org.helm.webeditor.Monomers.getMonomer(null, sym)`
   and record presence + identifier.
3. Build a modified `IMonomerLib` view that removes one known
   monomer symbol (e.g. omit `meI`) and adds a new sentinel
   symbol (e.g. `TestSentinelMonomer`) with a stable molfile.
   Use the platform's monomer-lib helper to create the modified
   view; do NOT mutate the global lib reference.
4. Call `rewriteLibraries(modifiedLib);` and wait for any
   library-rewrite promise the function returns to settle.
5. Re-snapshot the Pistoia dictionary: read
   `Monomers.getMonomer(null, 'meI')`,
   `Monomers.getMonomer(null, 'TestSentinelMonomer')`, and
   `Monomers.getMonomer(null, 'A')`.
6. Teardown — call `rewriteLibraries(lib)` (the original library
   from Setup) so subsequent scenarios run against the canonical
   monomer dictionary again.

Expected:
- After `rewriteLibraries(modifiedLib)`, `meI` is no longer
  resolvable from the Pistoia dictionary (the post-rewrite
  lookup returns null / undefined or yields the
  Pistoia-default placeholder, depending on Pistoia's
  unknown-symbol behaviour — the assertion is that the previous
  Datagrok-rewritten `meI` record is gone).
- The new `TestSentinelMonomer` IS resolvable post-rewrite; its
  molfile field matches the sentinel value supplied in the
  modified lib.
- Unaffected symbols (e.g. `A`) remain resolvable across the
  rewrite — `rewriteLibraries` replaces the dictionary contents
  in-place against the new lib, not against a diff.
- Teardown restores the original lib's reach into the Pistoia
  dictionary; the renderer / helm input / properties panel
  continue to resolve their fixture monomers in later
  scenarios.