import {category, test} from '@datagrok-libraries/utils/src/test';
import {
  _testMSAIsCorrect,
  _testTableIsNotEmpty,
} from './utils';
import {aligned1} from './test-data';

import * as DG from 'datagrok-api/dg';
//import * as grok from 'datagrok-api/grok';

export const _package = new DG.Package();

let table: DG.DataFrame;

category('peptides', async () => {
  //table = await grok.data.files.openTable('Demo:Files/bio/peptides.csv');
  table = DG.DataFrame.fromCsv(aligned1);
  const alignedSequencesColumn = table.getCol('AlignedSequence');

  test('MSA.test_table.is_not_empty', async () => {
    _testTableIsNotEmpty(table);
  });

  test('MSA.is_correct', async () => {
    _testMSAIsCorrect(alignedSequencesColumn);
  });
});
