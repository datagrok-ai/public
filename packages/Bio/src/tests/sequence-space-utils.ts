import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {expect} from '@datagrok-libraries/utils/src/test';
import {sequenceSpaceTopMenu} from '../package';

export async function _testSequenceSpaceReturnsResult(df: DG.DataFrame, algorithm: string, colName: string) {
  // await grok.data.detectSemanticTypes(df);
  const col: DG.Column = df.getCol(colName);
  const semType: string = await grok.functions.call('Bio:detectMacromolecule', {col: col});
  if (semType)
    col.semType = semType;

  const sp = await sequenceSpaceTopMenu(df, df.col(colName)!, algorithm, 'Levenshtein', true);
  expect(sp != null, true);
}