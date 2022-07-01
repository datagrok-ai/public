import {/*before, after, */after, category, test} from '@datagrok-libraries/utils/src/test';
import {
  _testViewerIsDrawing,
  _testDimensionalityReducer,
  _testPeptideSimilaritySpaceViewer,
  _testTableIsNotEmpty,
} from './utils';
import {DimensionalityReducer} from '@datagrok-libraries/ml/src/reduce-dimensionality';
import {cleanAlignedSequencesColumn} from '../utils/peptide-similarity-space';
import {aligned1} from './test-data';

import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {StringMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {computeWeights} from '../viewers/peptide-space-viewer';
import {_package} from '../package-test';

let table: DG.DataFrame;
let view: DG.TableView;

category('Peptide space', async () => {
  table = DG.DataFrame.fromCsv(aligned1);

  test('test_table.is_not_empty', async () => {
    _testTableIsNotEmpty(table);
  });

  test('PeptideSimilaritySpaceWidget.is_drawing', async () => {
    view = grok.shell.addTableView(table);
    await _testViewerIsDrawing(table, view);
  });

  const alignedSequencesColumn = table.getCol('AlignedSequence');

  test('test_deminsionality_reducer', async () => {
    const columnData = cleanAlignedSequencesColumn(alignedSequencesColumn);
    for (const method of DimensionalityReducer.availableMethods) {
      for (const measure of DimensionalityReducer.availableMetricsByType('String')) {
        test(`peptide_space.DimensinalityReducer.${method}.${measure}.is_numeric`, async () => {
          await _testDimensionalityReducer(columnData, method as StringMetrics, measure);
        });
      }
    }
  });

  test('test_peptide_similarity_space_viewer', async () => {
    for (const method of DimensionalityReducer.availableMethods) {
      for (const measure of DimensionalityReducer.availableMetricsByType('String')) {
        test(`peptide_space.PeptideSimilaritySpaceViewer.${method}.${measure}.is_proper`, async () => {
          await _testPeptideSimilaritySpaceViewer(table, alignedSequencesColumn, method, measure, 100);//, view);
        });
      }
    }
  });

  after(async () => {
    view?.close();
  });
});

category('Peptide Space Performance', () => {
  test('test_compute_weights_performance', async () => {
    const table = DG.DataFrame.fromCsv(await _package.files.readAsText('peptides_large.csv'));
    const results: {[key: string]: {[key: string]: {[key: string]: number}}} = {};
    const slice_volumes = [1, 2, 3, 4, 5, 7, 10];
    const methods = DimensionalityReducer.availableMethods;
    const metrics = DimensionalityReducer.availableMetricsByType('String');
    const total_runs = slice_volumes.length * methods.length * metrics.length;
    console.log('Started Peptide Space Performance benchmark...');

    let run = 0;
    for (const slice of slice_volumes) {
      const bitset = DG.BitSet.create(table.rowCount, (i) => i < slice * 1000);
      const table_slice = table.clone(bitset);
      const col = table_slice.getCol('sequence');
      const methodObj: {[key: string]: {[key: string]: number}} = {};

      for (const method of methods) {
        const measureObj: {[key: string]: number} = {};

        for (const metric of metrics) {
          console.log(`Run ${run++}/${total_runs}`);

          const start = new Date();
          await computeWeights(table_slice, method, metric, 100, col);
          const stop = new Date();

          measureObj[metric] = stop.getTime() - start.getTime();
        }
        methodObj[method] = measureObj;
      }
      results[`${slice}k`] = methodObj;
    }
    console.log('Peptide Space Performance benchmark finished...');
    console.log(results);
  });
});
