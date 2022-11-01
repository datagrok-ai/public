import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as bio from '@datagrok-libraries/bio';
import * as u from '@datagrok-libraries/utils';
import {TreeHelper} from '../utils/tree-helper';

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
u.category('treeCut', () => {
  const th: bio.ITreeHelper = new TreeHelper();

  const enum Tests {
    cut1 = 'cut1',
  }

  const data: {
    [test: string]:
      { tree: bio.NodeType, tgt: { cutLevel: number, leafNameList: string[][], tree?: bio.NodeType | null }[] }
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
          tree:
            {
              name: 'root',
              branch_length: 0.05,
              children: [],
              cuttedChildren: [
                {
                  name: 'root.cutted',
                  branch_length: 0.05,
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
                }
              ],
              cuttedLeafNameList: ['l1', 'l2', 'l3', 'l4', 'l5', 'l6'],
            } as bio.NodeCuttedType,
        },
        {
          cutLevel: 1.2,
          leafNameList: [['l1'], ['l2', 'l3'], ['l4', 'l5', 'l6']], tree: {
            name: 'root',
            branch_length: 0.1,
            children: [
              {
                name: 'node-l1-l2-l3',
                branch_length: 1,
                children: [
                  {
                    name: 'l1',
                    branch_length: 0.1,
                    children: [],
                    cuttedChildren: [{name: 'l1.cutted', branch_length: 0.9, children: []}],
                    cuttedLeafNameList: ['l1'],
                  } as bio.NodeCuttedType,
                  {
                    name: 'node-l2-l3',
                    branch_length: 0.1,
                    children: [],
                    cuttedChildren: [{
                      name: 'node-l2-l3.cutted',
                      branch_length: 0.1,
                      children: [
                        {name: 'l2', branch_length: 1.2},
                        {name: 'l3', branch_length: 1.4},
                      ],
                    }],
                    cuttedLeafNameList: ['l2', 'l3']
                  } as bio.NodeCuttedType
                ],
              },
              {
                name: 'node-l4-l5-l6',
                branch_length: 1.1,
                children: [],
                cuttedChildren: [{
                  name: 'node-l4-l5-l6.cutted',
                  branch_length: 0.5,
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
                }],
                cuttedLeafNameList: ['l4', 'l5', 'l6'],
              } as bio.NodeCuttedType
            ],
          }
        },
        {cutLevel: 1.5, leafNameList: [['l1'], ['l2'], ['l3'], ['l4', 'l5', 'l6']]},
        {cutLevel: 1.8, leafNameList: [['l1'], ['l2'], ['l3'], ['l4', 'l5'], ['l6']]},
        {cutLevel: 2.05, leafNameList: [['l1'], ['l2'], ['l3'], ['l4'], ['l5']]},
        {
          cutLevel: 2.6,
          leafNameList: [['l3'], ['l5']],
          tree: {
            name: 'root',
            branch_length: 0.1,
            children: [{
              name: 'node-l1-l2-l3',
              branch_length: 1,
              children: [{
                name: 'node-l2-l3',
                branch_length: 0.2,
                children: [{
                  name: 'l3',
                  branch_length: 1.3,
                  children: [],
                  cuttedChildren: [{
                    'name': 'l3.cutted',
                    'branch_length': 0.1,
                    'children': []
                  }],
                  cuttedLeafNameList: ['l3',],
                } as bio.NodeCuttedType,],
              }],
            }, {
              name: 'node-l4-l5-l6',
              branch_length: 1.6,
              children: [{
                name: 'node-l4-l5',
                branch_length: 0.3,
                children: [{
                  name: 'l5',
                  branch_length: 0.6,
                  children: [],
                  cuttedChildren: [
                    {
                      'name': 'l5.cutted',
                      'branch_length': 0.4,
                      'children': []
                    }
                  ],
                  cuttedLeafNameList: ['l5',],
                } as bio.NodeCuttedType,],
              }],
            }],
          }
        },
        {
          cutLevel: 3.0,
          leafNameList: [['l5']],
          tree: null
        },
        {
          cutLevel: 3.1,
          leafNameList: [],
          tree: null,
        },
      ]
    },
  };


  // -- treeCutAsLeafs --

  u.test('treeCutAsLeafs1_0.05', async () => {
    const v = data[Tests.cut1];
    _testTreeCutAsLeafs(v.tree, v.tgt[0].cutLevel, v.tgt[0].leafNameList);
  });

  u.test('treeCutAsLeafs1_1.20', async () => {
    const v = data[Tests.cut1];
    _testTreeCutAsLeafs(v.tree, v.tgt[1].cutLevel, v.tgt[1].leafNameList);
  });

  u.test('treeCutAsLeafs1_1.50', async () => {
    const v = data[Tests.cut1];
    _testTreeCutAsLeafs(v.tree, v.tgt[2].cutLevel, v.tgt[2].leafNameList);
  });

  u.test('treeCutAsLeafs1_1.80', async () => {
    const v = data[Tests.cut1];
    _testTreeCutAsLeafs(v.tree, v.tgt[3].cutLevel, v.tgt[3].leafNameList);
  });

  u.test('treeCutAsLeafs1_2.05', async () => {
    const v = data[Tests.cut1];
    _testTreeCutAsLeafs(v.tree, v.tgt[4].cutLevel, v.tgt[4].leafNameList);
  });

  u.test('treeCutAsLeafs1_2.60', async () => {
    const v = data[Tests.cut1];
    _testTreeCutAsLeafs(v.tree, v.tgt[5].cutLevel, v.tgt[5].leafNameList);
  });

  u.test('treeCutAsLeafs1_3.00', async () => {
    const v = data[Tests.cut1];
    _testTreeCutAsLeafs(v.tree, v.tgt[6].cutLevel, v.tgt[6].leafNameList);
  });

  u.test('treeCutAsLeafs1_3.10', async () => {
    const v = data[Tests.cut1];
    _testTreeCutAsLeafs(v.tree, v.tgt[7].cutLevel, v.tgt[7].leafNameList);
  });

  function _testTreeCutAsLeafs(tree: bio.NodeType, cutHeight: number, tgtLeafNameList: string[][]) {
    const resClusterList: bio.NodeType[] = th.treeCutAsLeaves(tree, cutHeight);
    const resLeafNameList = resClusterList
      .map((clusterRootNode) => th.getLeafList(clusterRootNode).map((l) => l.name));
    u.expectArray(resLeafNameList, tgtLeafNameList);
  }

  // -- treeCutAsTree --

  u.test('treeCutAsTree1_0.05', async () => {
    const v = data[Tests.cut1];
    _testTreeCutAsTree(v.tree, v.tgt[0].cutLevel, v.tgt[0].tree!);
  });

  u.test('treeCutAsTree1_1.20', async () => {
    const v = data[Tests.cut1];
    _testTreeCutAsTree(v.tree, v.tgt[1].cutLevel, v.tgt[1].tree!);
  });

  u.test('treeCutAsTree1_2.60', async () => {
    const v = data[Tests.cut1];
    _testTreeCutAsTree(v.tree, v.tgt[5].cutLevel, v.tgt[5].tree!);
  });

  function _testTreeCutAsTree(tree: bio.NodeType, cutHeight: number, tgtTree: bio.NodeType) {
    const resTree: bio.NodeType | null = th.treeCutAsTree(tree, cutHeight);
    //const newickStr = th.toNewick(resTree);

    u.expectObject(resTree!, tgtTree);
  }
});
