import {/*before, after, */category, test} from '@datagrok-libraries/utils/src/test';
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
import { StringMetrics } from '@datagrok-libraries/ml/src/typed-metrics';

export const _package = new DG.Package();

let table: DG.DataFrame;
let view: DG.TableView;

//const table = await grok.data.loadTable(`${_package.webRoot}files/aligned.csv`);
//'/home/www/data/dev/packages/data/peptides/aligned.csv');
//console.log(table);
//const table = await grok.data.files.openTable('Demo:Files/bio/peptides.csv');

category('peptides', async () => {
  /*before(async () => {
    console.log(['before']);
    // const text = await _package.files.readAsText('aligned.csv');
    // console.log([text]);
    // table = DG.DataFrame.fromCsv(text);

    // const path = `${_package.webRoot}files/aligned.csv`;
    // console.log([path]);
    // table = await grok.data.loadTable(path);
    // console.log([table]);

    table = await grok.data.files.openTable('Demo:Files/bio/peptides.csv');
    view = grok.shell.addTableView(table);
  });*/

  //table = await grok.data.files.openTable('Demo:Files/bio/peptides.csv');
  table = DG.DataFrame.fromCsv(aligned1);
  view = grok.shell.addTableView(table);

  test('peptide_space.test_table.is_not_empty', async () => {
    _testTableIsNotEmpty(table);
  });

  test('peptide_space.PeptideSimilaritySpaceWidget.is_drawing', async () => {
    await _testViewerIsDrawing(table, view);
  });

  const alignedSequencesColumn = table.getCol('AlignedSequence');
  const columnData = cleanAlignedSequencesColumn(alignedSequencesColumn);

  for (const method of DimensionalityReducer.availableMethods) {
    for (const measure of DimensionalityReducer.availableMetricsByType('String')) {
      test(`peptide_space.DimensinalityReducer.${method}.${measure}.is_numeric`, async () => {
        await _testDimensionalityReducer(columnData, method as StringMetrics, measure);
      });
    }
  }

  for (const method of DimensionalityReducer.availableMethods) {
    for (const measure of DimensionalityReducer.availableMetricsByType('String')) {
      test(`peptide_space.PeptideSimilaritySpaceViewer.${method}.${measure}.is_proper`, async () => {
        await _testPeptideSimilaritySpaceViewer(table, alignedSequencesColumn, method, measure, 100);//, view);
      });
    }
  }

  /*after(async () => {
    view.close();
    grok.shell.closeTable(table!);
  });*/
});
