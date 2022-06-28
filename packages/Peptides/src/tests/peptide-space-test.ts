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
import { computeWeights } from '../viewers/peptide-space-viewer';

export const _package = new DG.Package();

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
  })

  test('test_peptide_similarity_space_viewer', async () => {
    for (const method of DimensionalityReducer.availableMethods) {
      for (const measure of DimensionalityReducer.availableMetricsByType('String')) {
        test(`peptide_space.PeptideSimilaritySpaceViewer.${method}.${measure}.is_proper`, async () => {
          await _testPeptideSimilaritySpaceViewer(table, alignedSequencesColumn, method, measure, 100);//, view);
        });
      }
    }
  });

  test('test_compute_weights_performance', async () => {
    const table = DG.DataFrame.fromCsv(await _package.files.readAsText('steven_temp.csv'));
    for (const slice of [1, 2, 3]) {
      const bitset = DG.BitSet.create(table.rowCount, (i) => i < slice * 1000);
      const table_slice = table.clone(bitset);
      const col = table_slice.getCol('sequence');
      for (const method of DimensionalityReducer.availableMethods) {
        for (const measure of DimensionalityReducer.availableMetricsByType('String')) {
          const start = new Date();
          await computeWeights(table_slice, method, measure, 100, col);
          const stop = new Date();
          console.warn(`Table slice: ${slice}k; method: ${method}; measure: ${measure}; ` + 
            `Resulted in ${stop.getMilliseconds() - start.getMilliseconds()} ms`);
        }
      }
    }
  });

  after(async () => {
    view?.close();
  });
});
