import {category, test} from '@datagrok-libraries/utils/src/test';
import * as utils from './utils';
import {DimensionalityReducer} from '@datagrok-libraries/ml/src/reduce-dimensionality';
import {cleanAlignedSequencesColumn} from '../utils/peptide-similarity-space';
import {aligned1} from './test-data';

import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {StringMetrics} from '@datagrok-libraries/ml/src/typed-metrics';

let table: DG.DataFrame;
let view: DG.TableView;

category('Peptide space', async () => {
  table = DG.DataFrame.fromCsv(aligned1);

  test('test_table.is_not_empty', async () => {
    utils._testTableIsNotEmpty(table);
  });

  test('PeptideSimilaritySpaceWidget.is_drawing', async () => {
    view = grok.shell.addTableView(table);
    await utils._testViewerIsDrawing(table, view);
  });

  const alignedSequencesColumn = table.getCol('AlignedSequence');

  test('test_deminsionality_reducer', async () => {
    const columnData = cleanAlignedSequencesColumn(alignedSequencesColumn);
    for (const method of DimensionalityReducer.availableMethods) {
      for (const measure of DimensionalityReducer.availableMetricsByType('String')) {
        test(`peptide_space.DimensinalityReducer.${method}.${measure}.is_numeric`, async () => {
          await utils._testDimensionalityReducer(columnData, method as StringMetrics, measure);
        });
      }
    }
  });

  test('test_peptide_similarity_space_viewer', async () => {
    for (const method of DimensionalityReducer.availableMethods) {
      for (const measure of DimensionalityReducer.availableMetricsByType('String')) {
        test(`peptide_space.PeptideSimilaritySpaceViewer.${method}.${measure}.is_proper`, async () => {
          await utils._testPeptideSimilaritySpaceViewer(table, alignedSequencesColumn, method, measure, 100);//, view);
        });
      }
    }
  });
});
