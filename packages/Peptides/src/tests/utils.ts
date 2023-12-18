import * as DG from 'datagrok-api/dg';

import {expect} from '@datagrok-libraries/utils/src/test';
import {
  createDimensinalityReducingWorker,
} from '@datagrok-libraries/ml/src/workers/dimensionality-reducing-worker-creator';
import {StringMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {IReduceDimensionalityResult} from '@datagrok-libraries/ml/src/reduce-dimensionality';

export enum TEST_COLUMN_NAMES {
  SEQUENCE = 'sequence',
  ACTIVITY = 'activity',
  CLUSTER = 'cluster',
}

/**
 * Tests if a table has non zero rows and columns.
 * @param table Target table.
 */
export function _testTableIsNotEmpty(table: DG.DataFrame): void {
  expect(table.columns.length > 0 && table.rowCount > 0, true);
}

/**
 * Tests if dimensionality reducer works for both the method and the measure chosen.
 * @param columnData Strings to process.
 * @param method Embedding method.
 * @param measure Measure to apply to a pair of strings.
 */
export async function _testDimensionalityReducer(
  columnData: Array<string>, method: StringMetrics, measure: string): Promise<void> {
  const cyclesCount = 100;

  const reduceDimRes: IReduceDimensionalityResult = await createDimensinalityReducingWorker(
    {data: columnData, metric: measure as StringMetrics}, method, {cycles: cyclesCount});
  const embcols = reduceDimRes.embedding;

  const [X, Y] = embcols as Array<Float32Array>;

  expect(X.every((v) => v !== null && !Number.isNaN(v)), true);
  expect(Y.every((v) => v !== null && !Number.isNaN(v)), true);
}
