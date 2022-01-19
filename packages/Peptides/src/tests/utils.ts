import * as DG from 'datagrok-api/dg';

import {expect} from '@datagrok-libraries/utils/src/test';
import {PeptideSimilaritySpaceWidget} from '../utils/peptide-similarity-space';
import {
  createDimensinalityReducingWorker,
} from '@datagrok-libraries/ml/src/workers/dimensionality-reducing-worker-creator';

/**
 * Tests if peptide space viewer is drawing without exceptions.
 *
 * @param {DG.DataFrame} table Demo table.
 * @param {DG.TableView} view Demo view.
 */
export async function _testViewerIsDrawing(table: DG.DataFrame, view: DG.TableView) {
  let noException = true;

  try {
    const widget = new PeptideSimilaritySpaceWidget(table.getCol('AlignedSequence'), view);
    await widget.draw();
  } catch (error) {
    noException = false;
  }
  expect(noException, true);
}

/**
 * Tests if dimensionality reducer works for both the method and the measure chosen.
 *
 * @param {Array<string>} columnData Strings to process.
 * @param {string} method Embedding method.
 * @param {string} measure Measure to apply to a pair of strings.
 */
export async function _testDimensionalityReducer(columnData: Array<string>, method: string, measure: string) {
  const cyclesCount = 100;
  const embcols = await createDimensinalityReducingWorker(columnData, method, measure, cyclesCount);
  const [X, Y] = embcols as Array<Float32Array>;

  expect(X.every((v) => v !== null && v !== NaN), true);
  expect(Y.every((v) => v !== null && v !== NaN), true);
}
