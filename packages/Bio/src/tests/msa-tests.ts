import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {category, expect, expectArray, test} from '@datagrok-libraries/utils/src/test';

import {runKalign} from '../utils/multiple-sequence-alignment';
//import * as grok from 'datagrok-api/grok';

export const _package = new DG.Package();


category('MSA', async () => {
  //table = await grok.data.files.openTable('Demo:Files/bio/peptides.csv');
  const fromCsv = `seq
FWRWYVKHP
YNRWYVKHP
MWRSWYCKHP`;
  const toCsv = `seq
FWR-WYVKHP
YNR-WYVKHP
MWRSWYCKHP`;

  const longFromCsv = `seq
FWRWYVKHPFWRWYVKHPFWRWYVKHPFWRWYVKHPFWRWYVKHPFWRWYVKHPFWRWYVKHPFWRWYVKHP
YNRWYVKHPYNRWYVKHPYNRWYVKHPYNRWYVKHPYNRWYVKHPYNRWYVKHPYNRWYVKHPYNRWYVKHP
MWRSWYCKHPMWRSWYCKHPMWRSWYCKHPMWRSWYCKHPMWRSWYCKHPMWRSWYCKHPMWRSWYCKHPMWRSWYCKHP`;

  const longToCsv = `seq
FWR-WYVKHPFWR-WYVKHPFWR-WYVKHPFWR-WYVKHPFWR-WYVKHPFWR-WYVKHPFWR-WYVKHPFWR-WYVKHP
YNR-WYVKHPYNR-WYVKHPYNR-WYVKHPYNR-WYVKHPYNR-WYVKHPYNR-WYVKHPYNR-WYVKHPYNR-WYVKHP
MWRSWYCKHPMWRSWYCKHPMWRSWYCKHPMWRSWYCKHPMWRSWYCKHPMWRSWYCKHPMWRSWYCKHPMWRSWYCKHP`;

  // test('test_table.is_not_empty', async () => {
  //   await _testTableIsNotEmpty(table);
  // });

  test('isCorrect', async () => {
    await _testMsaIsCorrect(fromCsv, toCsv);
  });

  test('isCorrectLong', async () => {
    await _testMsaIsCorrect(longFromCsv, longToCsv);
  });
});

async function _testMsaIsCorrect(srcCsv: string, tgtCsv: string): Promise<void> {
  const srcDf: DG.DataFrame = DG.DataFrame.fromCsv(srcCsv);
  const tgtDf: DG.DataFrame = DG.DataFrame.fromCsv(tgtCsv);

  const srcCol: DG.Column = srcDf.getCol('seq')!;
  const semType: string = await grok.functions
    .call('Bio:detectMacromolecule', {col: srcCol}) as unknown as string;
  if (semType)
    srcCol.semType = semType;

  const tgtCol: DG.Column = tgtDf.getCol('seq')!;
  const msaCol: DG.Column = await runKalign(srcCol, true);
  expectArray(msaCol.toList(), tgtCol.toList());
}
