import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {category, test, expectObject} from '@datagrok-libraries/test/src/test';
import {ClusterMatrix} from '@datagrok-libraries/bio/src/trees';
import {NodeType} from '@datagrok-libraries/bio/src/trees';
import {TreeHelper} from '../utils/tree-helper';

const clusterMat1: ClusterMatrix =
  {mergeRow1: new Int32Array([-2, -1, -4]),
    mergeRow2: new Int32Array([-3, 1, 2]),
    heightsResult: new Float32Array([1, 2.5, 5.3333])};

const clusterMat2: ClusterMatrix =
  {mergeRow1: new Int32Array([-1, -2]),
    mergeRow2: new Int32Array([-3, 1]),
    heightsResult: new Float32Array([3, 4.5])};

const enum Tests{
    clust1,
    clust2
};

const testData: {[test: string]: {
    root: NodeType
}} = {
  [Tests.clust1]: {
    root: {
      name: '',
      children: [
        {
          name: '3',
          branch_length: 5.3333001136779785,
        },
        {
          name: '',
          children: [
            {
              name: '0',
              branch_length: 2.5,
            },
            {
              name: '',
              children: [
                {
                  name: '1',
                  branch_length: 1,
                },
                {
                  name: '2',
                  branch_length: 1,
                },
              ],
              branch_length: 1.5,
            },
          ],
          branch_length: 2.8333001136779785,
        },
      ],
      branch_length: 0,
    },
  },
  [Tests.clust2]: {
    root: {
      name: '',
      children: [
        {
          name: '1',
          branch_length: 4.5,
        },
        {
          name: '',
          children: [
            {
              name: '0',
              branch_length: 3,
            },
            {
              name: '2',
              branch_length: 3,
            },
          ],
          branch_length: 1.5,
        },
      ],
      branch_length: 0,
    },
  },
};

// Tests for converting matrix to tree
category('clusterMatrix', () => {
  const th = new TreeHelper();
  test('clustMatrixToTree1', async () => {
    const actualRoot = th.parseClusterMatrix(clusterMat1);
    const expectedRoot = testData[Tests.clust1].root;
    expectObject(actualRoot, expectedRoot);
  });

  test('clustMatrixToTree2', async () => {
    const actualRoot = th.parseClusterMatrix(clusterMat2);
    const expectedRoot = testData[Tests.clust2].root;
    expectObject(actualRoot, expectedRoot);
  });
});
