import {category, test} from '@datagrok-libraries/utils/src/test';
import {
  _testMSAIsCorrect,
  _testTableIsNotEmpty,
} from './utils';
import {aligned1} from './test-data';

import * as DG from 'datagrok-api/dg';
//import * as grok from 'datagrok-api/grok';

export const _package = new DG.Package();



category('MSA', async () => {
  //table = await grok.data.files.openTable('Demo:Files/bio/peptides.csv');
  const fromCsv = `seq
  FWRWYVKHP
  YNRWYVKHP
  MWRSWYCKHP`;
  const toCsv = `seq
  -F-W-R--W-Y-V-K-H-P
  -Y-N-R--W-Y-V-K-H-P
  -M-W-R-S-W-Y-C-K-H-P`;
  const table: DG.DataFrame = DG.DataFrame.fromCsv(fromCsv);
  const toTable: DG.DataFrame = DG.DataFrame.fromCsv(toCsv);
  const alignedSequencesColumn = toTable.getCol('seq');

  test('test_table.is_not_empty', async () => {
    await _testTableIsNotEmpty(table);
  });

  test('is_correct', async () => {
    await _testMSAIsCorrect(alignedSequencesColumn);
  });
});
