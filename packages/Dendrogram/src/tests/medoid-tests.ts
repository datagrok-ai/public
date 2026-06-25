import * as DG from 'datagrok-api/dg';

import {category, test, expect, expectArray} from '@datagrok-libraries/test/src/test';
import {TreeHelper} from '../utils/tree-helper';
import {DistanceMetric} from '@datagrok-libraries/bio/src/trees';
import {DistanceMatrixService} from '@datagrok-libraries/ml/src/distance-matrix';
import {NumberMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {DistanceAggregationMethods} from '@datagrok-libraries/ml/src/distance-matrix/types';

category('medoids', () => {
  const th = new TreeHelper();

  test('numeric-1d-two-clusters', async () => {
    // two well-separated clusters; the central member of each is the medoid (rank 1)
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.TYPE.FLOAT, 'x', [0, 1, 2, 10, 11, 12]),
    ]);
    const {rankByRow, avgDistByRow} = await th.calcMedoids(
      df, ['x'], [[0, 1, 2], [3, 4, 5]], DistanceMetric.Euclidean);
    expectArray(rankByRow, [2, 1, 3, 2, 1, 3]);
    expectArray(avgDistByRow, [1.5, 1, 1.5, 1.5, 1, 1.5]);
  });

  test('numeric-1d-rank-order', async () => {
    // |x_i - x_j| mean distances: [29.5, 25.75, 25.5, 25.75, 95.5] -> medoid is idx2 (value 6)
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.TYPE.FLOAT, 'x', [0, 5, 6, 7, 100]),
    ]);
    const {rankByRow, avgDistByRow} = await th.calcMedoids(
      df, ['x'], [[0, 1, 2, 3, 4]], DistanceMetric.Euclidean);
    expectArray(rankByRow, [4, 2, 1, 3, 5]); // ties (idx1, idx3) broken by ascending row index
    expect(avgDistByRow[2], 25.5);
  });

  test('singleton-and-pair', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.TYPE.FLOAT, 'x', [5, 100, 101]),
    ]);
    const {rankByRow, avgDistByRow} = await th.calcMedoids(
      df, ['x'], [[0], [1, 2]], DistanceMetric.Euclidean);
    expect(rankByRow[0], 1); // singleton is its own representative
    expect(avgDistByRow[0], 0);
    expectArray(rankByRow, [1, 1, 2]); // pair: first member on a tie
  });

  test('empty-cluster-ignored', async () => {
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.TYPE.FLOAT, 'x', [0, 1, 2]),
    ]);
    const {rankByRow} = await th.calcMedoids(
      df, ['x'], [[], [0, 1, 2]], DistanceMetric.Euclidean);
    expectArray(rankByRow, [2, 1, 3]);
  });

  test('multi-column-two-clusters', async () => {
    // exercises the per-column normalization + aggregation path
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromList(DG.TYPE.FLOAT, 'x', [0, 1, 2, 10, 11, 12]),
      DG.Column.fromList(DG.TYPE.FLOAT, 'y', [0, 1, 2, 10, 11, 12]),
    ]);
    const {rankByRow} = await th.calcMedoids(
      df, ['x', 'y'], [[0, 1, 2], [3, 4, 5]], DistanceMetric.Euclidean);
    expect(rankByRow[1], 1); // medoid of cluster A
    expect(rankByRow[4], 1); // medoid of cluster B
  });

  test('service-calcMedoids-direct', async () => {
    // validate the library method independently of feature encoding
    const values = [Float32Array.from([0, 1, 2, 10, 11, 12])];
    const service = new DistanceMatrixService(true, true);
    try {
      const res = await service.calcMedoids(values, [NumberMetricsNames.Difference],
        [[0, 1, 2], [3, 4, 5]], [{}], [1], DistanceAggregationMethods.EUCLIDEAN);
      expectArray(res[0].members, [0, 1, 2]);
      expectArray(res[0].ranks, [2, 1, 3]);
      expectArray(res[0].meanDistances, [1.5, 1, 1.5]);
      expect(res[1].ranks[1], 1); // member 4 is the medoid of the second cluster
    } finally {
      service.terminate();
    }
  });
});
