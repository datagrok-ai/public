---
feature: dendrogram
target_layer: apitest
coverage_type: smoke
priority: p1
realizes_atlas: [dendrogram.cross.tree-helper-cross-package]
realizes: [dendrogram]
produced_from: atlas-driven
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
realized_as:
  - dendrogram-tree-helper-cross-package-api-spec.ts
gate_verdicts:
  e:
    verdict: PASS
    cycle_id: 2026-06-03-dendrogram-automate-02
    timestamp: 2026-06-03T18:30:00Z
    failure_keys: []
  b:
    verdict: PASS
    cycle_id: 2026-06-03-dendrogram-automate-02
    timestamp: 2026-06-03T20:30:00Z
    spec_runs:
      - spec: dendrogram-tree-helper-cross-package-api-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 16
        failure_keys: []
---

# Dendrogram — Cross-package TreeHelper / DendrogramService consumption

Three public package functions are exercised end-to-end as the
cross-package consumption contract that Bio, BiostructureViewer,
HitTriage, and Peptides rely on:

- `getTreeHelper` — returns an `ITreeHelper` instance (the package's
  `TreeHelper`), invoked through
  `grok.functions.call('Dendrogram:getTreeHelper')`.
- `getDendrogramService` — returns an `IDendrogramService`, invoked
  through `grok.functions.call('Dendrogram:getDendrogramService')`;
  the implementation caches the instance under
  `window.$dendrogramService` so the same singleton is reused across
  consumer plugins.
- `DendrogramService.injectTreeForGrid` — the facade that consumer
  plugins call to attach a dendrogram tree as a grid neighbor; the
  attached neighbor surfaces the `.dendrogram-assign-clusters-bttn`
  magic-wand icon.

Scenario 1 is the smoke happy-path through all three functions on a
small synthetic DataFrame + 4-leaf tree, ending in a grid-neighbor
injection that surfaces the magic-wand assign-clusters affordance.
Scenario 2 is the singleton invariant — repeated `getDendrogramService`
calls return the same instance and the `window.$dendrogramService`
cache is populated as a side effect.

## Setup

1. Build a 4-row synthetic `DG.DataFrame` with a single string column
   `leaf` whose values are `[A, B, C, D]` — these will line up with the
   tree's leaf names so `injectTreeForGrid` can sync rows to leaves.
2. Build a 4-leaf `NodeType` tree by calling
   `TreeHelper.newickToDf('((A:1,B:1):1,(C:1,D:1):1):0;')` to obtain a
   tagged DataFrame, then retrieving the resulting root through the
   helper API (`TreeHelper.getNodeList` chain) — or, equivalently,
   construct a `NodeType` literal whose leaf `name` values are
   `[A, B, C, D]`. The Setup block leaves the choice to the test
   author; either path produces the same 4-leaf root used downstream.
3. Open a `TableView` on the synthetic DataFrame so a `DG.Grid` is
   available for `injectTreeForGrid` to attach onto.

## Scenarios

### Scenario 1: getTreeHelper + getDendrogramService + injectTreeForGrid happy path

Steps:
1. Call
   `await grok.functions.call('Dendrogram:getTreeHelper')` and
   capture the returned `treeHelper`.
2. Assert `treeHelper` is truthy and exposes the expected
   `ITreeHelper` surface — at minimum the methods `newickToDf`,
   `toNewick`, `getLeafList`, `getNodeList` exist as callable
   properties.
3. Call
   `await grok.functions.call('Dendrogram:getDendrogramService')` and
   capture the returned `dendrogramService`.
4. Assert `dendrogramService` is truthy and exposes
   `injectTreeForGrid` as a callable property.
5. Use `treeHelper` to materialize the 4-leaf root prepared in Setup
   (either via `newickToDf` parse-and-walk or by direct `NodeType`
   construction, per Setup step 2).
6. Invoke
   `dendrogramService.injectTreeForGrid(grid, treeRoot, 'leaf')` on
   the open `TableView`'s grid.
7. Await DOM settle and query the grid neighbor area for the
   `.dendrogram-assign-clusters-bttn` magic-wand icon.

Expected:
- `getTreeHelper` resolves to a non-null `ITreeHelper` whose method
  surface matches the contract documented at
  `public/packages/Dendrogram/src/utils/tree-helper.ts`.
- `getDendrogramService` resolves to a non-null `IDendrogramService`
  whose `injectTreeForGrid` method is callable.
- `injectTreeForGrid` returns without throwing; the call leaves the
  open `TableView` with a grid neighbor populated by a dendrogram
  canvas.
- The `.dendrogram-assign-clusters-bttn` magic-wand icon is present
  in the grid neighbor (per
  `public/packages/Dendrogram/src/viewers/inject-tree-for-grid2.ts#L74`).
- No console errors are emitted during the call sequence.

### Scenario 2: getDendrogramService singleton invariant

Steps:
1. Clear or capture the initial value of `window.$dendrogramService`.
2. Call
   `await grok.functions.call('Dendrogram:getDendrogramService')`,
   capture the returned instance as `svc1`.
3. Assert `window.$dendrogramService` is now set and is referentially
   equal to `svc1`.
4. Call
   `await grok.functions.call('Dendrogram:getDendrogramService')` a
   second time, capture as `svc2`.
5. Assert `svc2 === svc1` (referential equality — the singleton is
   reused, not re-created).
6. Assert `window.$dendrogramService === svc1` continues to hold
   after the second call.

Expected:
- The first `getDendrogramService` call populates the
  `window.$dendrogramService` window-scoped cache.
- Subsequent calls return the same instance — confirming the
  singleton invariant the atlas derives from
  `public/packages/Dendrogram/src/package.ts#L83-L84` and the
  cross-package consumption contract relies on (Bio,
  BiostructureViewer, HitTriage, Peptides all import and reuse this
  one instance).

## Notes

- The one DOM observable used in Scenario 1 (the magic-wand icon) is
  asserted as a single-anchor check, not a UI flow. The UI flow that
  follows from `injectTreeForGrid` — the Assign Clusters dialog,
  Ctrl+wheel zoom, double-click reset — is already owned by
  `assign-clusters.md` and isn't re-asserted here.
- See: `public/packages/Dendrogram/src/package.ts#L69` (getTreeHelper).
- See: `public/packages/Dendrogram/src/package.ts#L83` (getDendrogramService).
- See: `public/packages/Dendrogram/src/utils/dendrogram-service.ts#L10`
  (DendrogramService.injectTreeForGrid).
- See: `.claude/skills/grok-browser/references/dendrogram.md`
  `## tree-helper-cross-package`.
