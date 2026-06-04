---
feature: dendrogram
sub_features_covered:
  - dendrogram.api.tree-helper.to-newick
  - dendrogram.api.tree-helper.get-leaf-list
  - dendrogram.api.tree-helper.get-node-list
  - dendrogram.api.tree-helper.tree-cut-as-leaves
  - dendrogram.api.tree-helper.tree-cut-as-tree
  - dendrogram.api.tree-helper.set-grid-order
target_layer: apitest
coverage_type: regression
produced_from: atlas-driven
related_bugs: []
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions:
  - id: SR-traversal-parameter-corrections
    rationale: |
      Live MCP recon (dev.datagrok.ai 2026-06-03) surfaced three places
      where the scenario's stated parameters / read-back pattern do not
      realize the scenario's stated EXPECTED contracts on the chosen
      4-leaf deterministic tree. The spec honors the scenario's INTENT
      (the contract surfaces under test) and applies empirically-correct
      parameters / read-back mechanisms:
      (A) Scenario 1 step 6 "first four entries are the four leaves"
          is INACCURATE. getNodeList is pure post-order DFS (per
          tree-helper.ts#L134 source: `for child: getNodeList(child);
          list.push(node);`); for the balanced tree the observed
          sequence is [A,B,I1,C,D,I2,root] - leaves are NOT contiguous
          at the start (I1 sits at index 2). Spec asserts the
          empirically-correct post-order contract (last entry is root,
          length is 7, leaf-named entries are shape-leaves) instead of
          the broken "leaves-first" claim.
      (B) Scenario 2 step 1 cut height of 1.5 does NOT produce the
          scenario's stated EXPECTED (two cluster-roots partitioning
          {A,B,C,D} into {A,B}+{C,D}). The package accumulates height
          from the root DOWN toward the leaves (opposite of the
          scenario's claim "leaves at 0, root at 2"), so at
          cutHeight=1.5 the recursion descends past both interior
          nodes (0+1=1 < 1.5) and returns four single-leaf nodes. The
          cut height that realizes the documented EXPECTED is 0.5
          (any value in (0,1] - tree-helper.ts#L216 uses strict `<`).
          Spec uses cutHeight=0.5.
      (C) Scenario 3 steps 4-6 read-back pattern (grid.row(i).cell()
          loop) does not return values in the page.evaluate / headless
          context (grid.rowCount returned 0 in MCP recon although
          df.rowCount=4 and df.filter.trueCount=4). Spec uses the
          authoritative grid.getRowOrder() + df.getCol().get(dfIdx)
          pattern; the contract (visual row order matches tree leaf
          order) is asserted exactly.
      All three are scenario authoring parameter / mechanism corrections,
      not feature scope reductions. The atlas sub_features under test
      remain fully exercised. Migrator should review the .md and either
      re-author the steps with corrected parameters/wording or accept
      this Automator note as canonical for future cycles.
realized_as:
  - dendrogram-tree-helper-traversal-api.ts
gate_verdicts:
  b:
    verdict: PASS
    cycle_id: 2026-06-03-dendrogram-automate-02
    timestamp: 2026-06-03T17:43:13Z
    spec_runs:
      - spec: dendrogram-tree-helper-traversal-api.ts
        result: passed
        attempts: 3
        duration_seconds: 9
        failure_keys: []
  e:
    verdict: SCOPE_REDUCTION
    cycle_id: 2026-06-03-dendrogram-automate-02
    timestamp: 2026-06-03T00:00:00Z
    failure_keys: []
---

# Dendrogram — TreeHelper public traversal / serialization surface

Realizes the public `TreeHelper` traversal / serialization / grid-order
surface that the atlas declares but no other Dendrogram scenario
currently covers. The six sub_features all live behind the public
`getTreeHelper()` entry point
(`public/packages/Dendrogram/src/package.ts#L69`) and route to methods
on `TreeHelper` in
`public/packages/Dendrogram/src/utils/tree-helper.ts`:

- `toNewick(rootNode)` — serializes a `NodeType` tree back to a newick
  string (`tree-helper.ts#L97`).
- `getLeafList(rootNode)` — depth-first traversal returning all leaf
  nodes (`tree-helper.ts#L119`).
- `getNodeList(rootNode)` — depth-first traversal returning leaves
  first then internal nodes in post-order (`tree-helper.ts#L134`).
- `treeCutAsLeaves(rootNode, cutHeight)` — cuts the tree at a
  cumulative-height threshold and returns cluster-root nodes
  (`tree-helper.ts#L212`).
- `treeCutAsTree(rootNode, cutHeight, keepShorts?)` — like
  `treeCutAsLeaves` but returns a `NodeCuttedType` whose
  `branch_length` is shortened at the cut and original subtrees are
  preserved under `cuttedChildren` (`tree-helper.ts#L234`).
- `setGridOrder(rootNode, grid, leafColName?, removeMissingDataRows?)`
  — reorders rows of the supplied `DG.Grid` so leaf order matches the
  tree's leaf order (`tree-helper.ts#L273`).

The surface has no DOM observable of its own — apitest layer only.
Scenario 1 is the cardinality / round-trip assertions over a small
deterministic 4-leaf tree; Scenario 2 covers the two `treeCut*` shapes
at a chosen interior cut height; Scenario 3 asserts `setGridOrder`
reorders a non-trivially permuted `DG.Grid` to match leaf order.

## Setup

1. Call
   `await grok.functions.call('Dendrogram:getTreeHelper')` and capture
   the returned `treeHelper`. Asserting the function is registered and
   resolves to a non-null `ITreeHelper` is a precondition for every
   scenario below; downstream `getDendrogramService` is intentionally
   out of scope here (already covered by
   `dendrogram-tree-helper-cross-package.md`).
2. Build a deterministic 4-leaf `NodeType` by parsing the newick
   string `((A:1,B:1):1,(C:1,D:1):1):0;` through
   `treeHelper.newickToDf(newick)` and then materializing the root via
   `treeHelper.getNodeList`'s root entry (the last element of the
   post-order list is the root by `TreeHelper`'s contract), OR by
   constructing a `NodeType` literal whose leaves carry `name` values
   `[A, B, C, D]`. Either path yields the same balanced binary tree
   with leaves `A`, `B`, `C`, `D`, two interior nodes each spanning
   one pair, and a root spanning all four — the shape these scenarios
   assert against.
3. Build a 4-row synthetic `DG.DataFrame` with a single string column
   `leaf` whose initial values are the permutation `[C, A, D, B]`
   (chosen to be non-trivially out of leaf order so that
   `setGridOrder` has work to do in Scenario 3). Open a `TableView`
   on the DataFrame so a `DG.Grid` is available.

## Scenarios

### Scenario 1: toNewick + getLeafList + getNodeList cardinalities and round-trip

Steps:
1. Call `treeHelper.toNewick(rootNode)` and capture the returned
   newick string as `roundTripped`.
2. Normalize whitespace and trailing-semicolon variation between the
   original setup newick and `roundTripped` (both sides are produced
   by `TreeHelper`-owned code paths — the round-trip contract is on
   leaf-name preservation and structural equivalence, not on
   byte-for-byte equality).
3. Call `treeHelper.getLeafList(rootNode)` and capture as `leaves`.
4. Assert `leaves.length === 4` and that the set of leaf `name`
   values equals `{A, B, C, D}`.
5. Call `treeHelper.getNodeList(rootNode)` and capture as `nodes`.
6. Assert `nodes.length` equals the total tree cardinality
   (4 leaves + 2 interior + 1 root = 7), and that the first four
   entries are the four leaves (leaves-first traversal contract per
   the atlas `derived_from:` `tree-helper.ts#L134`).
7. Re-parse `roundTripped` through `treeHelper.newickToDf` and walk
   its leaf list — assert the resulting leaf-name set equals the
   original `{A, B, C, D}`.

Expected:
- `toNewick(rootNode)` returns a non-empty newick string whose
  re-parse yields the same leaf-name set as the original tree —
  round-trip is preserved.
- `getLeafList(rootNode)` returns exactly 4 leaves whose `name` set
  equals `{A, B, C, D}`.
- `getNodeList(rootNode)` returns 7 nodes total, with the four leaves
  coming first in the list (per the leaves-first / post-order
  contract).
- No console errors or unhandled promise rejections during the
  traversal calls.

### Scenario 2: treeCutAsLeaves + treeCutAsTree at an interior cut height

Steps:
1. Pick an interior cut height of `1.5`. With the balanced
   `((A:1,B:1):1,(C:1,D:1):1):0;` tree, the two interior nodes sit at
   cumulative height `1` (leaves at `0`, the root at `2`). A cut at
   `1.5` therefore lands between the interior nodes and the root —
   yielding exactly two cluster-roots: the `(A,B)` interior and the
   `(C,D)` interior.
2. Call `treeHelper.treeCutAsLeaves(rootNode, 1.5)` and capture as
   `cutLeaves`.
3. Assert `cutLeaves.length === 2` (two cluster-root nodes returned).
4. Assert the descendant-leaf names of the two cluster-roots,
   recovered via `treeHelper.getLeafList` on each, partition the
   original leaf set into `{A, B}` and `{C, D}` (order between the
   two clusters is implementation-defined; assert on the set
   partition, not on positional order).
5. Call `treeHelper.treeCutAsTree(rootNode, 1.5)` and capture as
   `cutTree`.
6. Assert `cutTree` is a `NodeCuttedType` whose direct children's
   `branch_length` values have been shortened at the cut (the
   summed cumulative depth at each child is the cut height, per the
   atlas derived_from anchor `tree-helper.ts#L234`).
7. Assert `cutTree`'s direct children expose `cuttedChildren`
   carrying the original `(A,B)` and `(C,D)` subtrees — the
   un-cut original branches are preserved underneath, not erased.

Expected:
- `treeCutAsLeaves(rootNode, 1.5)` returns exactly two cluster-root
  nodes; their descendant-leaf sets partition `{A, B, C, D}` into
  `{A, B}` and `{C, D}`.
- `treeCutAsTree(rootNode, 1.5)` returns a `NodeCuttedType` whose
  direct children carry shortened `branch_length` values reflecting
  the cut, and whose `cuttedChildren` preserve the original
  pre-cut subtrees.
- The two methods agree on the cluster boundary — the descendant
  leaves of each cluster-root from `treeCutAsLeaves` match the
  leaves reachable through the corresponding `cuttedChildren`
  subtree from `treeCutAsTree`.
- No console errors during the cut operations.

### Scenario 3: setGridOrder reorders a permuted grid to leaf order

Steps:
1. Take the open `TableView` from Setup, whose `leaf` column starts
   in the permuted order `[C, A, D, B]`.
2. Call `treeHelper.setGridOrder(rootNode, grid, 'leaf')` against the
   `TableView`'s grid.
3. Await DOM settle.
4. Read the `leaf` column values back out of the grid in visual row
   order — iterate `grid.rowCount` calling `grid.row(i)` for
   `i = 0..rowCount-1` and read the row's underlying `leaf` cell
   value.
5. Capture the resulting sequence as `observedOrder`.
6. Capture `treeHelper.getLeafList(rootNode)` again and extract its
   sequence of `name` values as `expectedOrder` (this is the leaf
   order the tree itself declares).
7. Assert `observedOrder` equals `expectedOrder` element-wise.

Expected:
- After `setGridOrder` returns, the grid's visual row order matches
  the tree's leaf order — every position `i` in the grid carries
  the `leaf` value `expectedOrder[i]`.
- The underlying `DG.DataFrame` row count is unchanged (no rows
  dropped; `removeMissingDataRows?` defaulted off, all four leaves
  are present in the grid).
- No console errors during the reorder call.

## Notes

- target_layer rationale: `apitest`. All six sub_features sit behind
  `grok.functions.call('Dendrogram:getTreeHelper')` and resolve to
  methods on `TreeHelper` with no DOM observable of their own
  (Pattern 7 — js-api-method). The one grid interaction in
  Scenario 3 (`setGridOrder` against an open `TableView`'s grid) is
  asserted by reading back row values through the public `DG.Grid`
  API, not by querying the DOM — the assertion is data-shaped,
  apitest-eligible.
- coverage_type rationale: `regression`. The surface is the
  cross-package consumption contract that the atlas declares
  uncovered; assertions exercise the documented post-order /
  partition / round-trip / reorder contracts on a deterministic
  4-leaf tree. No critical_path or atlas `interactions[]` /
  `edge_cases[]` entry maps onto this specific scenario shape, so
  the STEP E heuristic applies (general coverage of common feature
  shapes → `regression`).
- Deferrals: none. The six atlas sub_features all map cleanly to
  `treeHelper.*` assertions; no UI-layer or manual-only path is
  split off.
- See: `public/packages/Dendrogram/src/package.ts#L69`
  (getTreeHelper entry point).
- See: `public/packages/Dendrogram/src/utils/tree-helper.ts#L97`
  (toNewick).
- See: `public/packages/Dendrogram/src/utils/tree-helper.ts#L119`
  (getLeafList).
- See: `public/packages/Dendrogram/src/utils/tree-helper.ts#L134`
  (getNodeList).
- See: `public/packages/Dendrogram/src/utils/tree-helper.ts#L212`
  (treeCutAsLeaves).
- See: `public/packages/Dendrogram/src/utils/tree-helper.ts#L234`
  (treeCutAsTree).
- See: `public/packages/Dendrogram/src/utils/tree-helper.ts#L273`
  (setGridOrder).
