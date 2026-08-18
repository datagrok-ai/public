import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../spec-login';

test.use(specTestOptions);

test('Dendrogram: TreeHelper public traversal / serialization / grid-order surface (JS API)', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  await softStep('Scenario 1 - toNewick round-trip, getLeafList cardinality, getNodeList post-order cardinality on a 4-leaf balanced binary tree', async () => {
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

      const roundTripped: string = th.toNewick(root);

      const leaves: any[] = th.getLeafList(root);
      const leafNamesSet = leaves.map(l => l.name).sort();

      const nodes: any[] = th.getNodeList(root);
      const nodeNamesInOrder = nodes.map(n => n.name);
      const nodeNamesMultiset = nodes.map(n => n.name).sort();

      const rtDf: any = th.newickToDf(roundTripped);
      const rtLeafCol = rtDf.getCol('leaf');
      const rtNodeCol = rtDf.getCol('node');
      const rtLeafNames: string[] = [];
      for (let i = 0; i < rtDf.rowCount; i++) {
        if (rtLeafCol.get(i)) rtLeafNames.push(rtNodeCol.get(i));
      }
      const rtLeafNamesSet = rtLeafNames.sort();

      const lastIsRoot = nodes[nodes.length - 1]?.name === 'root';

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
        rtLeafCountIs4: rtLeafNames.length === 4,
        rtLeafNamesSet,
      };
    });

    expect(result.roundTrippedNonEmpty,
      'toNewick(root) returns a non-empty newick string')
      .toBe(true);

    expect(result.leavesLen,
      'getLeafList(root).length === 4 (one leaf per A/B/C/D)')
      .toBe(4);
    expect(result.leafNamesSet,
      'getLeafList(root) leaf names == {A,B,C,D}')
      .toEqual(['A', 'B', 'C', 'D']);

    expect(result.nodesLen,
      'getNodeList(root).length === 7 (4 leaves + 2 interior + 1 root)')
      .toBe(7);
    expect(result.nodeNamesMultiset,
      'getNodeList(root) name multiset == {A,B,C,D,I1,I2,root}')
      .toEqual(['A', 'B', 'C', 'D', 'I1', 'I2', 'root']);

    expect(result.lastIsRoot,
      'getNodeList(root) ends with the root (post-order: root is pushed last)')
      .toBe(true);

    expect(result.leafShapeOk,
      'getNodeList(root) entries: leaf-named nodes (A,B,C,D) have no children; interior-named nodes (I1,I2,root) have children')
      .toBe(true);

    expect(result.rtLeafCountIs4,
      'newickToDf(toNewick(root)) yields exactly 4 leaf rows (leaf==true)')
      .toBe(true);
    expect(result.rtLeafNamesSet,
      'newickToDf(toNewick(root)) leaf-name set == {A,B,C,D} (round-trip leaf preservation)')
      .toEqual(['A', 'B', 'C', 'D']);
  });

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

      const cutLeaves: any[] = th.treeCutAsLeaves(root, cutHeight);

      const clusterLeafSets = cutLeaves
        .map((cl: any) => th.getLeafList(cl).map((l: any) => l.name).sort());

      const sortedClusterLeafSets = [...clusterLeafSets].sort(
        (a, b) => a[0].localeCompare(b[0]));

      const cutTree: any = th.treeCutAsTree(root, cutHeight);
      const ctChildren = cutTree?.children || [];
      const ctChildrenCount = ctChildren.length;

      const ctChildrenShortened = ctChildren.every((c: any) => {
        return typeof c.branch_length === 'number' && c.branch_length === cutHeight;
      });

      const ctChildrenHaveCutted = ctChildren.every((c: any) =>
        Array.isArray(c.cuttedChildren) && c.cuttedChildren.length === 1);
      const cuttedChildrenLeafSets = ctChildren.map((c: any) => {
        const cuttedRoot = c.cuttedChildren?.[0];

        const leafNames: string[] = [];
        for (const childOfCutted of (cuttedRoot?.children || [])) {
          const sub = th.getLeafList(childOfCutted);
          for (const l of sub) leafNames.push(l.name);
        }
        return leafNames.sort();
      });
      const sortedCuttedChildrenLeafSets = [...cuttedChildrenLeafSets].sort(
        (a, b) => a[0].localeCompare(b[0]));

      const methodsAgree = JSON.stringify(sortedClusterLeafSets)
        === JSON.stringify(sortedCuttedChildrenLeafSets);

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

    expect(result.cutLeavesLen,
      'treeCutAsLeaves(root, 0.5) returns exactly 2 cluster-root nodes')
      .toBe(2);

    expect(result.sortedClusterLeafSets,
      'treeCutAsLeaves cluster-roots partition {A,B,C,D} into {A,B} and {C,D} (order between the two clusters is implementation-defined; assertion is on the partition)')
      .toEqual([['A', 'B'], ['C', 'D']]);

    expect(result.ctChildrenCount,
      'treeCutAsTree(root, 0.5) returns a node with exactly 2 direct children (one per cluster-root)')
      .toBe(2);
    expect(result.ctChildrenShortened,
      'direct children of cutTree have shortened branch_length == cutHeight (0.5) - shortened from the original 1.0 to reflect the cut depth (tree-helper.ts#L254)')
      .toBe(true);

    expect(result.ctChildrenHaveCutted,
      'each direct child of cutTree has cuttedChildren of length 1 (the wrapper holding the original pre-cut subtree, tree-helper.ts#L256)')
      .toBe(true);
    expect(result.sortedCuttedChildrenLeafSets,
      'cuttedChildren preserve the original {A,B} and {C,D} subtrees')
      .toEqual([['A', 'B'], ['C', 'D']]);

    expect(result.methodsAgree,
      'treeCutAsLeaves cluster boundaries match the leaves reachable through treeCutAsTree.cuttedChildren')
      .toBe(true);

    expect(result.sortedCuttedLeafNameLists,
      'cuttedLeafNameList recorded on each cluster-root matches the partition {A,B}, {C,D}')
      .toEqual([['A', 'B'], ['C', 'D']]);
  });

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

      grok.shell.closeAll();
      await new Promise(r => setTimeout(r, 400));

      const df: any = DG.DataFrame.fromColumns([
        DG.Column.fromList('string', 'leaf', ['C', 'A', 'D', 'B']),
      ]);
      df.name = 'tree-helper-traversal-set-grid-order';
      const tv: any = grok.shell.addTableView(df);
      await new Promise(r => setTimeout(r, 1200));

      const grid: any = tv.grid;
      const dfRowCountBefore = df.rowCount;

      let threw: string | false = false;
      let setRes: any = null;
      try {
        setRes = th.setGridOrder(root, grid, 'leaf');
      } catch (e: any) {
        threw = String(e?.message ?? e);
      }

      await new Promise(r => setTimeout(r, 500));

      const rowOrder = grid.getRowOrder();
      const rowOrderArr: number[] = [];
      for (let i = 0; i < (rowOrder?.length ?? 0); i++) rowOrderArr.push(rowOrder[i]);

      const leafCol = df.getCol('leaf');
      const observedOrder: string[] = rowOrderArr.map((dfIdx) => leafCol.get(dfIdx));

      const expectedOrder: string[] = th.getLeafList(root).map((l: any) => l.name);

      const dfRowCountAfter = df.rowCount;

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

    expect(result.dfRowCountAfter,
      'df.rowCount unchanged after setGridOrder (removeMissingDataRows defaulted off; all four leaves present)')
      .toBe(result.dfRowCountBefore);
    expect(result.dfRowCountAfter,
      'df.rowCount stays at 4 after the reorder (no rows dropped)')
      .toBe(4);

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

  await page.evaluate(() => {
    grok.shell.closeAll();
  });

  if (stepErrors.length > 0) {
    const summary = stepErrors.map(e => `  - ${e.step}: ${e.error}`).join('\n');
    throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
  }
});
