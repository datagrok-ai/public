import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test, expect, expectArray, expectObject} from '@datagrok-libraries/utils/src/test';
import {getLeafList, treeCutAsLeafs, treeCutAsTree} from '../utils/tree-helper';
import {NodeType} from '@datagrok-libraries/bio/src/types';

/*
R code to plot picture explaining tests

nwk <- '((l1:1,(l2:1.2,l3:1.4)node-l2-l3:0.2)node-l1-l2-l3:1,((l4:0.3,l5:0.7)node-l4-l5:0.3,l6:0.3)node-l4-l5-l6:1.6)root:0.1;';
tree <- ape::read.tree(text=nwk);
ape::plot.phylo(tree, show.node.label = TRUE, root.edge = TRUE)
axis(1);
abline(v=0.05, col='red')
abline(v=1.2,  col='green');
abline(v=1.5,  col='blue');

abline(v=1.8,  col='orange');
abline(v=2.05, col='violet');
abline(v=2.6,  col='brown');

*/
category('treeCut', () => {
  const enum Tests {
    cut1 = 'cut1',
  }

  const data: {
    [test: string]: { tree: NodeType, tgt: { cutLevel: number, leafNameList: string[][], tree?: NodeType | null }[] }
  } = {
    [Tests.cut1]: {
      tree: {
        name: 'root',
        branch_length: 0.1,
        children: [
          {
            name: 'node-l1-l2-l3',
            branch_length: 1,
            children: [
              {name: 'l1', branch_length: 1},
              {
                name: 'node-l2-l3',
                branch_length: 0.2,
                children: [
                  {name: 'l2', branch_length: 1.2},
                  {name: 'l3', branch_length: 1.4},
                ]
              }
            ]
          },
          {
            name: 'node-l4-l5-l6',
            branch_length: 1.6,
            children: [
              {
                name: 'node-l4-l5',
                branch_length: 0.3,
                children: [
                  {name: 'l4', branch_length: 0.3},
                  {name: 'l5', branch_length: 1.0},
                ],
              },
              {name: 'l6', branch_length: 0.3},
            ]
          }
        ],
      },
      tgt: [
        {
          cutLevel: 0.05,
          leafNameList: [['l1', 'l2', 'l3', 'l4', 'l5', 'l6'],],
        },
        {cutLevel: 1.2, leafNameList: [['l1'], ['l2', 'l3'], ['l4', 'l5', 'l6']]},
        {cutLevel: 1.5, leafNameList: [['l1'], ['l2'], ['l3'], ['l4', 'l5', 'l6']]},
        {cutLevel: 1.8, leafNameList: [['l1'], ['l2'], ['l3'], ['l4', 'l5'], ['l6']]},
        {cutLevel: 2.05, leafNameList: [['l1'], ['l2'], ['l3'], ['l4'], ['l5']]},
        {
          cutLevel: 2.6,
          leafNameList: [['l3'], ['l5']],
          tree: {
            name: 'root',
            branch_length: 0.1,
            children: [
              {
                name: 'node-l1-l2-l3',
                branch_length: 1,
                children: [{
                  name: 'node-l2-l3',
                  branch_length: 0.2,
                  children: [{name: 'l3', branch_length: 0.07, cuttedLeafNameList: ['l3']},],
                }]
              },
              {
                name: 'node-l4-l5-l6',
                branch_length: 1.6,
                children: [{
                  name: 'node-l4-l5',
                  branch_length: 0.3,
                  children: [{name: 'l5', branch_length: 0.07, cuttedLeafNameList: ['l5']},],
                }]
              }]
          }
        },
        {
          cutLevel: 3.0,
          leafNameList: [],
          tree: null
        }
      ]
    },
  };

  function treeToNewick(node: NodeType): string {
    const isLeaf = !node.children || node.children.length == 0;

    if (isLeaf) {
      return Array<string>().concat(
        node.name,
        node.branch_length ? `:${node.branch_length}` : [],
      ).join('');
    } else {
      const childrenText = node.children!.map((childNode) => treeToNewick(childNode)).join(',');
      return Array<string>().concat(
        `(${childrenText})`,
        node.name,
        node.branch_length ? `:${node.branch_length}` : [],
      ).join('');
    }
  }

  // -- treeCutAsLeafs --

  test('treeCutAsLeafs1_0.05', async () => {
    const v = data[Tests.cut1];
    _testTreeCutAsLeafs(v.tree, v.tgt[0].cutLevel, v.tgt[0].leafNameList);
  });

  test('treeCutAsLeafs1_1.20', async () => {
    const v = data[Tests.cut1];
    _testTreeCutAsLeafs(v.tree, v.tgt[1].cutLevel, v.tgt[1].leafNameList);
  });

  test('treeCutAsLeafs1_1.50', async () => {
    const v = data[Tests.cut1];
    _testTreeCutAsLeafs(v.tree, v.tgt[2].cutLevel, v.tgt[2].leafNameList);
  });

  test('treeCutAsLeafs1_1.80', async () => {
    const v = data[Tests.cut1];
    _testTreeCutAsLeafs(v.tree, v.tgt[3].cutLevel, v.tgt[3].leafNameList);
  });

  test('treeCutAsLeafs1_2.05', async () => {
    const v = data[Tests.cut1];
    _testTreeCutAsLeafs(v.tree, v.tgt[4].cutLevel, v.tgt[4].leafNameList);
  });

  test('treeCutAsLeafs1_2.60', async () => {
    const v = data[Tests.cut1];
    _testTreeCutAsLeafs(v.tree, v.tgt[5].cutLevel, v.tgt[5].leafNameList);
  });

  test('treeCutAsLeafs1_3.00', async () => {
    const v = data[Tests.cut1];
    _testTreeCutAsLeafs(v.tree, v.tgt[6].cutLevel, v.tgt[6].leafNameList);
  });

  function _testTreeCutAsLeafs(tree: NodeType, cutHeight: number, tgtLeafNameList: string[][]) {
    const resClusterList: NodeType[] = treeCutAsLeafs(tree, cutHeight);
    const resLeafNameList = resClusterList
      .map((clusterRootNode) => getLeafList(clusterRootNode).map((l) => l.name));
    expectArray(resLeafNameList, tgtLeafNameList);
  }

  // -- treeCutAsTree --

  test('treeCutAsTree1_2.60', async () => {
    const v = data[Tests.cut1];
    _testTreeCutAsTree(v.tree, v.tgt[5].cutLevel, v.tgt[5].tree!);
  });

  function _testTreeCutAsTree(tree: NodeType, cutHeight: number, tgtTree: NodeType) {
    const resTree: NodeType | null = treeCutAsTree(tree, cutHeight);
    let k = 11;

    expectObject(resTree!, tgtTree);
  }
});