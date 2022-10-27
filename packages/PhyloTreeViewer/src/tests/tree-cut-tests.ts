import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as bio from '@datagrok-libraries/bio';

import {after, before, category, test, expect, expectArray, expectObject} from '@datagrok-libraries/utils/src/test';
import {getLeafList, treeCutAsLeafs, treeCutAsTree} from '../utils/tree-helper';
import {Node} from '@datagrok-libraries/bio/src/types';
import {ExtNode} from '../types';

/*
R code to plot picture explaining tests

nwk <- '((l1:1,(l2:1.2,l3:1.4)node-l2-l3:0.2)node-l1-l2-l3:1,((l4:0.3,l5:0.7)node-l4-l5:0.3,l6:0.3)node-l4-l5-l6:1.6)root:0.1;';
tree <- ape::read.tree(text=nwk);
ape::plot.phylo(tree, show.node.label = TRUE, root.edge = TRUE, x.lim=c(0,4.0));
axis(1);
abline(v=0.05, col='red')
abline(v=1.2,  col='green');
abline(v=1.5,  col='blue');

abline(v=1.8,  col='orange');
abline(v=2.05, col='violet');
abline(v=2.6,  col='brown');
abline(v=3.0,  col='gray');

*/
category('treeCut', () => {
  const enum Tests {
    cut1 = 'cut1',
  }

  const data: {
    [test: string]: { tree: bio.Node, tgt: { cutLevel: number, leafNameList: string[][], tree?: Node | null }[] }
  } = {
    [Tests.cut1]: {
      tree: new bio.Node(
        'root',
        0.1,
        [
          new bio.Node(
            'node-l1-l2-l3',
            1,
            [
              new bio.Node('l1', 1),
              new bio.Node(
                'node-l2-l3',
                0.2,
                [
                  new bio.Node('l2', 1.2),
                  new bio.Node('l3', 1.4),
                ])
            ]),
          new bio.Node(
            'node-l4-l5-l6',
            1.6,
            [
              new bio.Node(
                'node-l4-l5',
                0.3,
                [
                  new bio.Node('l4', 0.3),
                  new bio.Node('l5', 1.0),
                ]),
              new bio.Node('l6', 0.3),
            ]
          )
        ]
      ),
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
          tree: new bio.Node(
            'root',
            0.1,
            [
              new bio.Node(
                'node-l1-l2-l3',
                1,
                [new bio.Node(
                  'node-l2-l3',
                  0.2,
                  [new ExtNode('l3', 0.07, ['l3']),],
                )]
              ),
              new bio.Node(
                'node-l4-l5-l6',
                1.6,
                [new bio.Node(
                  'node-l4-l5',
                  0.3,
                  [new ExtNode('l5', 0.07, ['l5']),],
                )]
              )
            ]
          )
        },
        {
          cutLevel: 3.0,
          leafNameList: [['l5']],
          tree: null,
        },
        {
          cutLevel: 3.1,
          leafNameList: [],
          tree: null,
        }
      ]
    },
  };


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

  test('treeCutAsLeafs1_3.10', async () => {
    const v = data[Tests.cut1];
    _testTreeCutAsLeafs(v.tree, v.tgt[7].cutLevel, v.tgt[7].leafNameList);
  });

  function _testTreeCutAsLeafs(tree: bio.Node, cutHeight: number, tgtLeafNameList: string[][]) {
    const resClusterList: bio.Node[] = treeCutAsLeafs(tree, cutHeight);
    const resLeafNameList = resClusterList
      .map((clusterRootNode) => getLeafList(clusterRootNode).map((l) => l.name));
    expectArray(resLeafNameList, tgtLeafNameList);
  }

  // -- treeCutAsTree --

  test('treeCutAsTree1_2.60', async () => {
    const v = data[Tests.cut1];
    _testTreeCutAsTree(v.tree, v.tgt[5].cutLevel, v.tgt[5].tree!);
  });

  function _testTreeCutAsTree(tree: bio.Node, cutHeight: number, tgtTree: bio.Node) {
    const resTree: bio.Node | null = treeCutAsTree(tree, cutHeight);
    let k = 11;

    expectObject(resTree!, tgtTree);
  }
});