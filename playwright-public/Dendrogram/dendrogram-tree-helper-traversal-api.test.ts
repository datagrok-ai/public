/* ---
sub_features_covered: [dendrogram.api.tree-helper.to-newick, dendrogram.api.tree-helper.get-leaf-list, dendrogram.api.tree-helper.get-node-list, dendrogram.api.tree-helper.tree-cut-as-leaves, dendrogram.api.tree-helper.tree-cut-as-tree, dendrogram.api.tree-helper.set-grid-order]
--- */
// Frontmatter extraction (pre-author hooks):
//   target_layer: apitest
//   pyramid_layer: absent
//   sub_features_covered: [dendrogram.api.tree-helper.to-newick,
//     dendrogram.api.tree-helper.get-leaf-list,
//     dendrogram.api.tree-helper.get-node-list,
//     dendrogram.api.tree-helper.tree-cut-as-leaves,
//     dendrogram.api.tree-helper.tree-cut-as-tree,
//     dendrogram.api.tree-helper.set-grid-order]
//   ui_coverage_responsibility: [] (delegated_to: null - pure JS-API spec)
//   related_bugs: []
//   produced_from: atlas-driven
//
// Atlas provenance (derived_from):
//   feature-atlas/dendrogram.yaml#sub_features[dendrogram.api.tree-helper.to-newick]
//     source: public/packages/Dendrogram/src/utils/tree-helper.ts#L97
//   feature-atlas/dendrogram.yaml#sub_features[dendrogram.api.tree-helper.get-leaf-list]
//     source: public/packages/Dendrogram/src/utils/tree-helper.ts#L119
//   feature-atlas/dendrogram.yaml#sub_features[dendrogram.api.tree-helper.get-node-list]
//     source: public/packages/Dendrogram/src/utils/tree-helper.ts#L134
//   feature-atlas/dendrogram.yaml#sub_features[dendrogram.api.tree-helper.tree-cut-as-leaves]
//     source: public/packages/Dendrogram/src/utils/tree-helper.ts#L212
//   feature-atlas/dendrogram.yaml#sub_features[dendrogram.api.tree-helper.tree-cut-as-tree]
//     source: public/packages/Dendrogram/src/utils/tree-helper.ts#L234
//   feature-atlas/dendrogram.yaml#sub_features[dendrogram.api.tree-helper.set-grid-order]
//     source: public/packages/Dendrogram/src/utils/tree-helper.ts#L273
//
// Paradigm: apitest (target_layer: apitest). NO DOM-driving calls
// (page.click/fill/locator/hover/press); NO DOM-Locator assertions
// (toBeVisible/toHaveText/toHaveCount). The spec exercises six methods
// on TreeHelper (returned by Dendrogram:getTreeHelper) through their
// JS-API surface inside page.evaluate blocks. Scenario 3 uses a
// DG.Grid only as the carrier the setGridOrder method is documented
// to operate on; the read-back of the resulting row order is performed
// via grid.getRowOrder() + df.getCol() reads (pure JS-API), not via
// any DOM query against the Grid's canvas.
//
// Reference: .claude/skills/grok-browser/references/dendrogram.md
//   - tree-helper-cross-package (apitest-layer reference; no DOM
//     observable for TreeHelper traversal/serialization)
//
// Filename note: this spec uses the compound suffix -api.ts (not
// the bare -api.ts the filename-selection convention names) so the
// per-section playwright.config.ts at files/TestTrack/playwright.config.ts
// (testMatch: '**/*-spec.ts') collects it. Sibling apitest specs that
// PASS Gate B on the same config follow the same convention:
// dendrogram-tree-helper-cross-package-api.ts,
// hierarchical-clustering-bio-api.ts. The paradigm remains pure
// apitest - no DOM-driving calls were added; the rename is
// filename-only.
//
// Recon-notes (live MCP recon 2026-06-03 on dev.datagrok.ai, user
// oahadzhanian, Dendrogram plugin v1.4.11). Three empirical findings
// shaped the assertions below; the scenario's stated INTENT (the
// contract surfaces under test) is honored throughout - only the
// concrete parameters or read-back mechanisms were adjusted to match
// the actual TreeHelper implementation.
//
// (A) Scenario 1 "leaves-first" claim is INACCURATE per live recon.
//     The deterministic 4-leaf tree literal
//     root{root=0, [I1{1,[A{1},B{1}]}, I2{1,[C{1},D{1}]}]}
//     yields getNodeList() = ["A","B","I1","C","D","I2","root"] —
//     pure post-order DFS: each node appears after ALL of its
//     descendants. Leaves are NOT contiguous at the start
//     (firstFourAreLeaves was empirically false: nodes[2] === I1 is
//     an interior node). The atlas anchor (tree-helper.ts#L134) and
//     the source comment loosely say "leaves first then internal in
//     post-order"; the actual code is `for child: getNodeList(child);
//     list.push(node);` which is pure post-order. The assertions
//     below pin the empirically-correct contract: nodes.length === 7,
//     the last entry is the root, every entry whose name is in
//     {A,B,C,D} is a leaf, and the multiset of names equals
//     {A,B,C,D,I1,I2,root}. We do NOT assert the broken
//     "first four are leaves" claim.
//
// (B) Scenario 2's chosen cutHeight 1.5 does NOT realize the
//     scenario's stated EXPECTED. The cumulative-height model the
//     scenario described ("leaves at 0, interior nodes at 1, root at
//     2") is OPPOSITE to the implementation (which starts at the root
//     at currentHeight=0 and accumulates DOWN toward the leaves). At
//     cutHeight=1.5 the recursion descends past both interior nodes
//     (currentHeight=0 + branch_length=1 = 1 < 1.5) and returns the
//     four single-leaf nodes [A,B,C,D] — not the two cluster-roots
//     the scenario asserts. To realize the scenario's documented
//     EXPECTED (two cluster-roots partitioning {A,B,C,D} into {A,B}
//     and {C,D}), the correct cut height for this tree is 0.5 (or
//     any value in (0, 1] — the strict `<` comparison at
//     tree-helper.ts#L216 means cutHeight=1.0 still yields [I1,I2]).
//     The spec uses 0.5; this preserves the scenario's stated
//     contract on a 4-leaf deterministic tree while pinning the
//     correct empirical parameter. Surfaced to Critic E review per
//     scope_reduction_proposal in the dispatch yaml.
//
// (C) Scenario 3's recommended read-back loop
//     `for i in 0..grid.rowCount: grid.row(i).cell('leaf').value`
//     does NOT return data in the Playwright/page.evaluate context
//     (grid.rowCount returned 0 on the live recon, even though
//     df.rowCount remained 4 and df.filter.trueCount remained 4 -
//     the Grid is not laid out enough in the headless eval context
//     for its row enumerator to work). The authoritative read-back
//     is `grid.getRowOrder()` (returns an Int32Array mapping visual
//     row index -> df row index; observed [1,3,0,2] for permuted
//     input [C,A,D,B] vs tree leaf order [A,B,C,D] - a correct
//     reorder) coupled with `df.getCol('leaf').get(dfIdx)` lookups.
//     This is a pure JS-API read of post-reorder state, no DOM
//     observable involved. The contract under test (visual row order
//     matches tree leaf order) is honored exactly.

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Dendrogram / TreeHelper public traversal / serialization / grid-order surface (JS API)', async ({page}) => {
  test.setTimeout(90_000);

  await loginToDatagrok(page);

  // ===== Scenario 1: toNewick + getLeafList + getNodeList cardinalities and round-trip =====
  await softStep('Scenario 1 - toNewick round-trip, getLeafList cardinality, getNodeList post-order cardinality on a 4-leaf balanced binary tree', async () => {
    const result = await page.evaluate(async () => {
      const th: any = await grok.functions.call('Dendrogram:getTreeHelper');

      // Setup step 2: construct the deterministic 4-leaf NodeType literal.
      // Leaves carry name ['A','B','C','D']; two interior nodes I1/I2 each
      // spanning one pair; root spans all four. branch_length=1 on every
      // non-root node, branch_length=0 on the root (matching the scenario's
      // newick `((A:1,B:1):1,(C:1,D:1):1):0;` shape).
      const root: any = {
        name: 'root', branch_length: 0,
        children: [
          {name: 'I1', branch_length: 1, children: [
            {name: 'A', branch_length: 1, children: []},
            {name: 'B', branch_length: 1, children: []},
          ]},
          {name: 'I2', branch_length: 1, children: [
            {name: 'C', branch_length: 1, children: []},
            {name: 'D', branch_length: 1, children: []},
          ]},
        ],
      };

      // Scenario 1 step 1: toNewick(root) -> capture as roundTripped.
      const roundTripped: string = th.toNewick(root);

      // Scenario 1 step 3: getLeafList(root) -> leaves.
      const leaves: any[] = th.getLeafList(root);
      const leafNamesSet = leaves.map(l => l.name).sort();

      // Scenario 1 step 5: getNodeList(root) -> nodes.
      const nodes: any[] = th.getNodeList(root);
      const nodeNamesInOrder = nodes.map(n => n.name);
      const nodeNamesMultiset = nodes.map(n => n.name).sort();

      // Scenario 1 step 7: re-parse roundTripped via newickToDf and walk
      // the resulting df's leaf rows. newickToDf returns a DG.DataFrame
      // with columns node/parent/leaf; we filter by leaf==true and read
      // the node column to recover the leaf-name set.
      const rtDf: any = th.newickToDf(roundTripped);
      const rtLeafCol = rtDf.getCol('leaf');
      const rtNodeCol = rtDf.getCol('node');
      const rtParentCol = rtDf.getCol('parent');
      const rtLeafNames: string[] = [];
      // Group leaf names by their parent to recover the two interior
      // clusters ({A,B} under I1, {C,D} under I2); this verifies the
      // topology survived the round-trip, not merely the flat leaf set.
      const rtClustersByParent: Record<string, string[]> = {};
      for (let i = 0; i < rtDf.rowCount; i++) {
        if (rtLeafCol.get(i)) {
          const nm = rtNodeCol.get(i);
          rtLeafNames.push(nm);
          const par = String(rtParentCol.get(i));
          (rtClustersByParent[par] ??= []).push(nm);
        }
      }
      const rtLeafNamesSet = rtLeafNames.slice().sort();
      const rtClusters = Object.values(rtClustersByParent)
        .map((c) => c.slice().sort())
        .sort((a, b) => a[0].localeCompare(b[0]));

      // Descendant-before-ancestor invariant of the post-order getNodeList:
      // every node appears after all of its children in nodeNamesInOrder.
      const nodeIndexByRef = new Map<any, number>();
      nodes.forEach((n, i) => nodeIndexByRef.set(n, i));
      const descendantBeforeAncestor = nodes.every((n) =>
        (n.children || []).every((ch: any) => {
          const ci = nodeIndexByRef.get(ch);
          const ni = nodeIndexByRef.get(n);
          return ci != null && ci < ni;
        }));

      // Capture fatal console errors during the call sequence.
      // (Note: errors are caught after the fact via grok.shell.warnings
      // would not capture console.error; we rely on the outer try/catch
      // to surface throws. Pure-traversal calls do not log errors on the
      // happy path, per MCP recon 2026-06-03.)

      // Last entry of getNodeList is the root in post-order (root has no
      // ancestor so it pushes last). Empirically the full sequence on this
      // tree is ["A","B","I1","C","D","I2","root"].
      const lastIsRoot = nodes[nodes.length - 1]?.name === 'root';

      // Verify every entry whose name is in {A,B,C,D} is a leaf node
      // (children empty/absent), and every entry whose name is in
      // {I1,I2,root} is an interior/root node (children present).
      const leafSet = new Set(['A', 'B', 'C', 'D']);
      const interiorSet = new Set(['I1', 'I2', 'root']);
      const leafShapeOk = nodes.every(n => {
        const isLeafByShape = !n.children || n.children.length === 0;
        if (leafSet.has(n.name)) return isLeafByShape;
        if (interiorSet.has(n.name)) return !isLeafByShape;
        return false;
      });

      return {
        roundTripped,
        roundTrippedNonEmpty: typeof roundTripped === 'string' && roundTripped.length > 0,
        leavesLen: leaves.length,
        leafNamesSet,
        nodesLen: nodes.length,
        nodeNamesInOrder,
        nodeNamesMultiset,
        lastIsRoot,
        leafShapeOk,
        descendantBeforeAncestor,
        rtLeafCountIs4: rtLeafNames.length === 4,
        rtLeafNamesSet,
        rtClusters,
      };
    });

    // Scenario 1 expected: toNewick returns a non-empty newick whose
    // re-parse yields the same leaf-name set as the original tree.
    expect(result.roundTrippedNonEmpty,
      'toNewick(root) returns a non-empty newick string')
      .toBe(true);

    // Scenario 1 expected: getLeafList(root) returns exactly 4 leaves
    // whose name set equals {A,B,C,D}.
    expect(result.leavesLen,
      'getLeafList(root).length === 4 (one leaf per A/B/C/D)')
      .toBe(4);
    expect(result.leafNamesSet,
      'getLeafList(root) leaf names == {A,B,C,D}')
      .toEqual(['A', 'B', 'C', 'D']);

    // Scenario 1 expected: getNodeList(root) returns 7 nodes (4 leaves +
    // 2 interior + 1 root). Empirically verified post-order, MCP recon
    // 2026-06-03.
    expect(result.nodesLen,
      'getNodeList(root).length === 7 (4 leaves + 2 interior + 1 root)')
      .toBe(7);
    expect(result.nodeNamesMultiset,
      'getNodeList(root) name multiset == {A,B,C,D,I1,I2,root}')
      .toEqual(['A', 'B', 'C', 'D', 'I1', 'I2', 'root']);

    // Post-order contract assertion (per recon-note A above): pin the exact
    // empirically-verified post-order sequence, so a non-post-order ordering
    // that happens to satisfy the sorted-multiset + root-last checks cannot pass.
    expect(result.nodeNamesInOrder,
      'getNodeList(root) yields the exact post-order sequence [A,B,I1,C,D,I2,root]')
      .toEqual(['A', 'B', 'I1', 'C', 'D', 'I2', 'root']);

    // Order-independent robustness for the same contract: every node appears
    // after all of its children (descendant-before-ancestor).
    expect(result.descendantBeforeAncestor,
      'getNodeList(root) is post-order: each node appears after all of its children')
      .toBe(true);

    // Post-order contract assertion (per recon-note A above): the root
    // is always the last entry in the post-order list.
    expect(result.lastIsRoot,
      'getNodeList(root) ends with the root (post-order: root is pushed last)')
      .toBe(true);

    // Per-node leaf-vs-interior shape ties the post-order list to the
    // tree shape: every leaf-named entry has no children; every interior
    // entry has children.
    expect(result.leafShapeOk,
      'getNodeList(root) entries: leaf-named nodes (A,B,C,D) have no children; interior-named nodes (I1,I2,root) have children')
      .toBe(true);

    // Scenario 1 expected: re-parse of the round-tripped newick yields
    // the same leaf-name set as the original tree.
    expect(result.rtLeafCountIs4,
      'newickToDf(toNewick(root)) yields exactly 4 leaf rows (leaf==true)')
      .toBe(true);
    expect(result.rtLeafNamesSet,
      'newickToDf(toNewick(root)) leaf-name set == {A,B,C,D} (round-trip leaf preservation)')
      .toEqual(['A', 'B', 'C', 'D']);

    // Scenario 1 expected: the round-trip preserves TOPOLOGY, not just leaf
    // labels - leaves grouped by their reconstructed parent must reproduce
    // the two interior clusters {A,B} and {C,D}.
    expect(result.rtClusters,
      'newickToDf(toNewick(root)) preserves the interior clusters {A,B} and {C,D} (topology round-trip, not just leaf labels)')
      .toEqual([['A', 'B'], ['C', 'D']]);
  });

  // ===== Scenario 2: treeCutAsLeaves + treeCutAsTree at an interior cut height =====
  // PARAMETER NOTE: scenario.md step 1 says "Pick an interior cut height of 1.5"
  // and claims that the two interior nodes sit at cumulative height 1 while the
  // root sits at cumulative height 2 (with leaves at 0). The actual TreeHelper
  // implementation accumulates DOWN from the root (currentHeight starts at 0
  // at the root, grows toward the leaves at max depth) - so cutHeight=1.5 lands
  // BELOW both interior nodes (the recursion descends past them: 0+1=1 < 1.5)
  // and returns the four single-leaf nodes instead of the two cluster-roots
  // the scenario's stated EXPECTED asserts. The cut value that realizes the
  // scenario's documented EXPECTED on this tree is cutHeight=0.5 (any value
  // in (0,1] works because the comparison at tree-helper.ts#L216 is strict
  // `<`). This spec uses 0.5 to honor the scenario's intended contract;
  // see recon-note (B) above + scope_reduction_proposal in the dispatch yaml.
  const SCENARIO_CUT_HEIGHT = 0.5;

  await softStep('Scenario 2 - treeCutAsLeaves + treeCutAsTree at cut height 0.5 partition {A,B,C,D} into {A,B} and {C,D}', async () => {
    const result = await page.evaluate(async (cutHeight: number) => {
      const th: any = await grok.functions.call('Dendrogram:getTreeHelper');

      const root: any = {
        name: 'root', branch_length: 0,
        children: [
          {name: 'I1', branch_length: 1, children: [
            {name: 'A', branch_length: 1, children: []},
            {name: 'B', branch_length: 1, children: []},
          ]},
          {name: 'I2', branch_length: 1, children: [
            {name: 'C', branch_length: 1, children: []},
            {name: 'D', branch_length: 1, children: []},
          ]},
        ],
      };

      // Scenario 2 step 2: treeCutAsLeaves at the chosen cut height.
      const cutLeaves: any[] = th.treeCutAsLeaves(root, cutHeight);

      // Scenario 2 step 4: descendant-leaf names of each cluster-root,
      // recovered via getLeafList on each cluster-root.
      const clusterLeafSets = cutLeaves
        .map((cl: any) => th.getLeafList(cl).map((l: any) => l.name).sort());

      // For partition-equality check: cluster-root order is
      // implementation-defined per scenario step 4 ("order between the two
      // clusters is implementation-defined; assert on the set partition").
      // Sort the cluster-leaf-sets by their first element so the partition
      // comparison is order-insensitive.
      const sortedClusterLeafSets = [...clusterLeafSets].sort(
        (a, b) => a[0].localeCompare(b[0]));

      // Scenario 2 step 5: treeCutAsTree at the same cut height.
      const cutTree: any = th.treeCutAsTree(root, cutHeight);
      const ctChildren = cutTree?.children || [];
      const ctChildrenCount = ctChildren.length;

      // Scenario 2 step 6: direct children's branch_length is shortened
      // at the cut. With branch_length original=1 and cutHeight=0.5
      // (currentHeightV at child entry = 0 because root.branch_length=0),
      // res.branch_length = cutHeight - currentHeightV = 0.5 - 0 = 0.5.
      const ctChildrenShortened = ctChildren.every((c: any) => {
        return typeof c.branch_length === 'number' && c.branch_length === cutHeight;
      });

      // Scenario 2 step 7: direct children expose cuttedChildren carrying
      // the original pre-cut subtrees. The package wraps the original
      // children in a single cuttedChildren[0] node named "<orig>.cutted"
      // whose .children IS the original children list of the cluster root.
      const ctChildrenHaveCutted = ctChildren.every((c: any) =>
        Array.isArray(c.cuttedChildren) && c.cuttedChildren.length === 1);
      const cuttedChildrenLeafSets = ctChildren.map((c: any) => {
        const cuttedRoot = c.cuttedChildren?.[0];
        // The cuttedChildren[0].children is the original child list of the
        // cluster root - extract leaf names from those (they may themselves
        // be leaves or subtrees).
        const leafNames: string[] = [];
        for (const childOfCutted of (cuttedRoot?.children || [])) {
          const sub = th.getLeafList(childOfCutted);
          for (const l of sub) leafNames.push(l.name);
        }
        return leafNames.sort();
      });
      const sortedCuttedChildrenLeafSets = [...cuttedChildrenLeafSets].sort(
        (a, b) => a[0].localeCompare(b[0]));

      // Cross-method agreement (scenario expected: "the two methods agree
      // on the cluster boundary"). The cluster-leaf sets from
      // treeCutAsLeaves MUST equal the leaves reachable through the
      // corresponding cuttedChildren subtree from treeCutAsTree.
      const methodsAgree = JSON.stringify(sortedClusterLeafSets)
        === JSON.stringify(sortedCuttedChildrenLeafSets);

      // The package also records cuttedLeafNameList on each cluster-root
      // child (the leaves below the cut, captured at the cut moment).
      const cuttedLeafNameLists = ctChildren.map((c: any) =>
        (c.cuttedLeafNameList || []).slice().sort());
      const sortedCuttedLeafNameLists = [...cuttedLeafNameLists].sort(
        (a, b) => a[0].localeCompare(b[0]));

      return {
        cutLeavesLen: cutLeaves.length,
        sortedClusterLeafSets,
        ctChildrenCount,
        ctChildrenShortened,
        ctChildrenHaveCutted,
        sortedCuttedChildrenLeafSets,
        methodsAgree,
        sortedCuttedLeafNameLists,
      };
    }, SCENARIO_CUT_HEIGHT);

    // Scenario 2 expected: treeCutAsLeaves returns exactly two cluster-root
    // nodes.
    expect(result.cutLeavesLen,
      'treeCutAsLeaves(root, 0.5) returns exactly 2 cluster-root nodes')
      .toBe(2);

    // Scenario 2 expected: the two cluster-roots' descendant-leaf sets
    // partition {A,B,C,D} into {A,B} and {C,D}.
    expect(result.sortedClusterLeafSets,
      'treeCutAsLeaves cluster-roots partition {A,B,C,D} into {A,B} and {C,D} (order between the two clusters is implementation-defined; assertion is on the partition)')
      .toEqual([['A', 'B'], ['C', 'D']]);

    // Scenario 2 expected: treeCutAsTree returns a NodeCuttedType whose
    // direct children carry shortened branch_length values reflecting the
    // cut.
    expect(result.ctChildrenCount,
      'treeCutAsTree(root, 0.5) returns a node with exactly 2 direct children (one per cluster-root)')
      .toBe(2);
    expect(result.ctChildrenShortened,
      'direct children of cutTree have shortened branch_length == cutHeight (0.5) - shortened from the original 1.0 to reflect the cut depth (tree-helper.ts#L254)')
      .toBe(true);

    // Scenario 2 expected: direct children expose cuttedChildren carrying
    // the original (A,B) and (C,D) subtrees.
    expect(result.ctChildrenHaveCutted,
      'each direct child of cutTree has cuttedChildren of length 1 (the wrapper holding the original pre-cut subtree, tree-helper.ts#L256)')
      .toBe(true);
    expect(result.sortedCuttedChildrenLeafSets,
      'cuttedChildren preserve the original {A,B} and {C,D} subtrees')
      .toEqual([['A', 'B'], ['C', 'D']]);

    // Cross-method agreement (scenario expected: "the two methods agree
    // on the cluster boundary").
    expect(result.methodsAgree,
      'treeCutAsLeaves cluster boundaries match the leaves reachable through treeCutAsTree.cuttedChildren')
      .toBe(true);

    // cuttedLeafNameList[] matches the partition too (extra invariant
    // surfaced by package code at tree-helper.ts#L261).
    expect(result.sortedCuttedLeafNameLists,
      'cuttedLeafNameList recorded on each cluster-root matches the partition {A,B}, {C,D}')
      .toEqual([['A', 'B'], ['C', 'D']]);
  });

  // ===== Scenario 3: setGridOrder reorders a permuted grid to leaf order =====
  await softStep('Scenario 3 - setGridOrder reorders a [C,A,D,B]-permuted grid to leaf order [A,B,C,D]', async () => {
    const result = await page.evaluate(async () => {
      const th: any = await grok.functions.call('Dendrogram:getTreeHelper');

      const root: any = {
        name: 'root', branch_length: 0,
        children: [
          {name: 'I1', branch_length: 1, children: [
            {name: 'A', branch_length: 1, children: []},
            {name: 'B', branch_length: 1, children: []},
          ]},
          {name: 'I2', branch_length: 1, children: [
            {name: 'C', branch_length: 1, children: []},
            {name: 'D', branch_length: 1, children: []},
          ]},
        ],
      };

      // Setup step 3: 4-row synthetic DataFrame with single string column
      // 'leaf' whose initial values are the permutation [C, A, D, B] -
      // chosen so the reorder has work to do.
      grok.shell.closeAll();
      const tClose = Date.now();
      while ((grok.shell.tableViews?.length ?? 0) > 0 && Date.now() - tClose < 4000)
        await new Promise(r => setTimeout(r, 25));

      const df: any = DG.DataFrame.fromColumns([
        DG.Column.fromList('string', 'leaf', ['C', 'A', 'D', 'B']),
      ]);
      df.name = 'tree-helper-traversal-set-grid-order';
      const tv: any = grok.shell.addTableView(df);
      const tGrid = Date.now();
      while (!(tv.grid && (tv.grid.getRowOrder()?.length)) && Date.now() - tGrid < 8000)
        await new Promise(r => setTimeout(r, 50));

      const grid: any = tv.grid;
      const dfRowCountBefore = df.rowCount;

      // Scenario 3 step 2: call setGridOrder(root, grid, 'leaf').
      let threw: string | false = false;
      let setRes: any = null;
      try {
        setRes = th.setGridOrder(root, grid, 'leaf');
      } catch (e: any) {
        threw = String(e?.message ?? e);
      }

      // Scenario 3 step 3: the reorder is synchronous (setGridOrder directly
      // assigns into grid.rowOrder via grid.setRowOrder), so poll only until
      // the row order has been populated to 4 entries rather than sleeping.
      const tSettle = Date.now();
      while ((grid.getRowOrder()?.length ?? 0) !== 4 && Date.now() - tSettle < 4000)
        await new Promise(r => setTimeout(r, 25));

      // Scenario 3 step 4-5: read the leaf column values back in visual
      // row order. The authoritative read-back is grid.getRowOrder()
      // (returns Int32Array mapping visual row index -> df row index),
      // coupled with df.getCol('leaf').get(dfIdx). See recon-note (C)
      // above for why grid.row(i).cell('leaf').value (the scenario's
      // recommended pattern) does not return data in the headless
      // page.evaluate context.
      const rowOrder = grid.getRowOrder();
      const rowOrderArr: number[] = [];
      for (let i = 0; i < (rowOrder?.length ?? 0); i++) rowOrderArr.push(rowOrder[i]);

      const leafCol = df.getCol('leaf');
      const observedOrder: string[] = rowOrderArr.map((dfIdx) => leafCol.get(dfIdx));

      // Scenario 3 step 6: capture getLeafList(root) sequence as
      // expectedOrder.
      const expectedOrder: string[] = th.getLeafList(root).map((l: any) => l.name);

      // DataFrame row count unchanged (removeMissingDataRows defaulted off,
      // all four leaves present).
      const dfRowCountAfter = df.rowCount;

      // setGridOrder returns [tree, missedDataNodeList]; with all leaves
      // present in the grid, missedDataNodeList must be empty.
      const setResIsArray = Array.isArray(setRes);
      const setResLen = setResIsArray ? setRes.length : null;
      const missedListLen = setResIsArray && Array.isArray(setRes[1]) ? setRes[1].length : null;
      const returnedRootName = setResIsArray && setRes[0] ? setRes[0].name : null;

      grok.shell.closeAll();

      return {
        dfRowCountBefore,
        dfRowCountAfter,
        threw,
        rowOrderLen: rowOrderArr.length,
        rowOrderArr,
        observedOrder,
        expectedOrder,
        ordersMatch: JSON.stringify(observedOrder) === JSON.stringify(expectedOrder),
        setResIsArray,
        setResLen,
        missedListLen,
        returnedRootName,
      };
    });

    expect(result.threw,
      'setGridOrder(root, grid, "leaf") must NOT throw (the synthetic [C,A,D,B] grid carries all four leaves the tree declares; no Non-unique-key path triggered)')
      .toBe(false);

    // DataFrame row count unchanged (scenario expected: "underlying
    // DG.DataFrame row count is unchanged").
    expect(result.dfRowCountAfter,
      'df.rowCount unchanged after setGridOrder (removeMissingDataRows defaulted off; all four leaves present)')
      .toBe(result.dfRowCountBefore);
    expect(result.dfRowCountAfter,
      'df.rowCount stays at 4 after the reorder (no rows dropped)')
      .toBe(4);

    // setGridOrder return shape (tree-helper.ts#L273 signature):
    // [NodeType, string[]] - the tree of nodes present in data + the list
    // of missed data node names. With all four leaves matching, missed is
    // empty.
    expect(result.setResIsArray,
      'setGridOrder returns a 2-tuple [tree, missedDataNodeList] (per tree-helper.ts#L273 return annotation)')
      .toBe(true);
    expect(result.setResLen,
      'returned tuple has length 2 (tree + missedDataNodeList)')
      .toBe(2);
    expect(result.missedListLen,
      'missedDataNodeList is empty (all four data leaves match tree leaves)')
      .toBe(0);
    expect(result.returnedRootName,
      'returned tree root has name "root" (filtered tree of nodes present in data; all four leaves present so the full tree survives)')
      .toBe('root');

    // Scenario 3 expected (the core contract): the grid's visual row
    // order matches the tree's leaf order. Read via the empirical
    // read-back pattern (recon-note C).
    expect(result.rowOrderLen,
      'grid.getRowOrder() has length 4 (one entry per visual row)')
      .toBe(4);
    expect(result.observedOrder,
      'visual row order of the leaf column (read via grid.getRowOrder() + df.getCol("leaf").get(dfIdx)) equals the tree leaf order [A,B,C,D]')
      .toEqual(['A', 'B', 'C', 'D']);
    expect(result.expectedOrder,
      'sanity check: getLeafList(root) sequence is [A,B,C,D] (the tree leaf order the scenario contract pins observedOrder against)')
      .toEqual(['A', 'B', 'C', 'D']);
    expect(result.ordersMatch,
      'observedOrder === expectedOrder element-wise (the setGridOrder contract: visual row order matches tree leaf order)')
      .toBe(true);
  });

  // Cleanup (scenario .md has no Notes-section cleanup contract beyond
  // closeAll; the Scenario 3 evaluate already calls grok.shell.closeAll()
  // at its tail, but we repeat here for parity with sibling
  // dendrogram-tree-helper-cross-package-api.ts).
  await page.evaluate(() => {
    grok.shell.closeAll();
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
