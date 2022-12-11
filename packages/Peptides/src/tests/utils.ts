import * as DG from 'datagrok-api/dg';

import {expect} from '@datagrok-libraries/utils/src/test';
import {PeptideSimilaritySpaceWidget, createPeptideSimilaritySpaceViewer} from '../utils/peptide-similarity-space';
import {
  createDimensinalityReducingWorker,
  IReduceDimensionalityResult,
} from '@datagrok-libraries/ml/src/workers/dimensionality-reducing-worker-creator';
import {StringMetrics} from '@datagrok-libraries/ml/src/typed-metrics';

/**
 * Tests if a table has non zero rows and columns.
 *
 * @param {DG.DataFrame} table Target table.
 */
export function _testTableIsNotEmpty(table: DG.DataFrame): void {
  expect(table.columns.length > 0 && table.rowCount > 0, true);
}

/**
 * Tests if peptide space viewer is drawing without exceptions.
 *
 * @param {DG.DataFrame} table Demo table.
 * @param {DG.TableView} view Demo view.
 */
export async function _testViewerIsDrawing(table: DG.DataFrame, view: DG.TableView): Promise<void> {
  const widget = new PeptideSimilaritySpaceWidget(table.getCol('AlignedSequence'), view);
  await widget.draw();
}

/**
 * Tests if dimensionality reducer works for both the method and the measure chosen.
 *
 * @param {Array<string>} columnData Strings to process.
 * @param {string} method Embedding method.
 * @param {string} measure Measure to apply to a pair of strings.
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

/**
 * Tests if PeptideSimilaritySpaceViewer works for both the method and the measure chosen.
 *
 * @export
 * @param {DG.DataFrame} table Table.
 * @param {DG.Column} alignedSequencesColumn Aligned sequences column.
 * @param {string} method Embedding method.
 * @param {string} measure Strings similarity measure.
 * @param {number} cyclesCount Number of embedding iterations.
 * @param {(DG.TableView | null)} view Viewer to show graphics on.
 * @param {(string | null)} [activityColumnName] Name of column with activity.
 */
export async function _testPeptideSimilaritySpaceViewer(table: DG.DataFrame, alignedSequencesColumn: DG.Column,
  method: string, measure: string, cyclesCount: number): Promise<void> {
  const viewer = await createPeptideSimilaritySpaceViewer(
    table, method, measure, cyclesCount, alignedSequencesColumn, undefined);
  const df = viewer.dataFrame;

  const axesNames = ['~X', '~Y', '~MW'];
  const axes = axesNames.map((v) => df.getCol(v).getRawData() as Float32Array);

  for (const ax of axes)
    expect(ax.every((v) => v !== null && !Number.isNaN(v)), true);
}
