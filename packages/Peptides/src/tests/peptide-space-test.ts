import {category, test} from '@datagrok-libraries/utils/src/test';
import {
  _testViewerIsDrawing,
  _testDimensionalityReducer,
} from './utils';
import {DimensionalityReducer} from '@datagrok-libraries/ml/src/reduce-dimensionality';
import {cleanAlignedSequencesColumn} from '../utils/peptide-similarity-space';

import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';


category('peptides', async () => {
  const table = await grok.data.files.openTable('Demo:TestJobs:Files:DemoFiles/bio/peptides.csv');
  const view = grok.shell.v as DG.TableView;

  test('peptides.peptide_space.viewer_is_drawing', async () => {
    await _testViewerIsDrawing(table, view);
  });

  const alignedSequencesColumn = table.getCol('AlignedSequence');
  const columnData = cleanAlignedSequencesColumn(alignedSequencesColumn);

  for (const method of DimensionalityReducer.availableMethods) {
    for (const measure of DimensionalityReducer.availableMetrics) {
      test(`peptides.dimensionality_reducer.${method}.${measure}`, async () => {
        await _testDimensionalityReducer(columnData, method, measure);
      });
    }
  }
});
