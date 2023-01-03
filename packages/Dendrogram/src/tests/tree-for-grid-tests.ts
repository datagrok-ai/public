import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {category, test, expect, expectArray, expectObject} from '@datagrok-libraries/utils/src/test';

import {NodeType} from '@datagrok-libraries/bio';
import {markupNode, MarkupNodeType} from '../viewers/tree-renderers/markup';

category('treeForGrid', () => {
  const enum Tests {
    markup1 = 'markup1',
  }

  const data: { [test: string]: { tree: NodeType, tgt: MarkupNodeType } } =
    {
      [Tests.markup1]: {
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
        tgt: {
          name: 'root',
          branch_length: 0.1,
          children: [
            {
              name: 'node-l1-l2-l3',
              branch_length: 1,
              children: [
                {name: 'l1', branch_length: 1, index: 0, subtreeLength: 1} as MarkupNodeType,
                {
                  name: 'node-l2-l3',
                  branch_length: 0.2,
                  children: [
                    {name: 'l2', branch_length: 1.2, index: 1, subtreeLength: 1.2} as MarkupNodeType,
                    {name: 'l3', branch_length: 1.4, index: 2, subtreeLength: 1.4} as MarkupNodeType,
                  ],
                  minIndex: 1, maxIndex: 2, subtreeLength: 1.6
                } as MarkupNodeType,
              ],
              minIndex: 0, maxIndex: 2, subtreeLength: 2.6
            } as MarkupNodeType,
            {
              name: 'node-l4-l5-l6',
              branch_length: 1.6,
              children: [
                {
                  name: 'node-l4-l5',
                  branch_length: 0.3,
                  children: [
                    {name: 'l4', branch_length: 0.3, index: 3, subtreeLength: 0.3} as MarkupNodeType,
                    {name: 'l5', branch_length: 1.0, index: 4, subtreeLength: 1.0} as MarkupNodeType,
                  ],
                  minIndex: 3, maxIndex: 4, subtreeLength: 1.3
                },
                {name: 'l6', branch_length: 0.3, index: 5, subtreeLength: 0.3} as MarkupNodeType,
              ],
              minIndex: 3, maxIndex: 5, subtreeLength: 2.9
            } as MarkupNodeType,
          ],
          minIndex: 0, maxIndex: 5, subtreeLength: 3.0
        } as MarkupNodeType
      }
    };

  test('markup1', async () => {
    const tree: NodeType = JSON.parse(JSON.stringify(data[Tests.markup1].tree)); // Deep copy
    const tgt: MarkupNodeType = data[Tests.markup1].tgt;

    markupNode(tree);

    expectObject(tree, tgt);
  });
});
