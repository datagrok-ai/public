import {category, test} from '@datagrok-libraries/utils/src/test';
import {
  _testViewerIsDrawing,
  _testDimensionalityReducer,
  _testPeptideSimilaritySpaceViewer,
} from './utils';
import {DimensionalityReducer} from '@datagrok-libraries/ml/src/reduce-dimensionality';
import {cleanAlignedSequencesColumn} from '../utils/peptide-similarity-space';

import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

//import {_package} from '../package';

category('peptides', async () => {
  /*console.log([_package.files]);
  const path = await _package.files.readAsText('aligned.csv');
  console.log([path]);
  const table = DG.DataFrame.fromCsv(path);*/
  //const table = await grok.data.loadTable('/home/www/data/dev/packages/data/peptides/aligned.csv');
  //`${_package.webRoot}files/aligned.csv`);
  //console.log(table);
  const table = await grok.data.files.openTable('Demo:Files/bio/peptides.csv');
  const view = grok.shell.v as DG.TableView;

  test('PeptideSimilaritySpaceWidget.is_drawing', async () => {
    await _testViewerIsDrawing(table, view);
  });

  const alignedSequencesColumn = table.getCol('AlignedSequence');
  const columnData = cleanAlignedSequencesColumn(alignedSequencesColumn);

  for (const method of DimensionalityReducer.availableMethods) {
    for (const measure of DimensionalityReducer.availableMetrics) {
      test(`DimensinalityReducer.${method}.${measure}.is_numeric`, async () => {
        await _testDimensionalityReducer(columnData, method, measure);
      });
    }
  }

  for (const method of DimensionalityReducer.availableMethods) {
    for (const measure of DimensionalityReducer.availableMetrics) {
      test(`PeptideSimilaritySpaceViewer.${method}.${measure}.is_proper`, async () => {
        await _testPeptideSimilaritySpaceViewer(table, alignedSequencesColumn, method, measure, 100);//, view);
      });
    }
  }
});
