import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import { expect } from '@datagrok-libraries/utils/src/test';
import { sequenceSpaceTopMenu } from '../package';

export async function _testSequenceSpaceReturnsResult(df: DG.DataFrame, algorithm: string, colName: string) {
  await grok.data.detectSemanticTypes(df);
  const sp = await sequenceSpaceTopMenu(df, df.col(colName)!, algorithm, 'Levenshtein', true);
  expect(sp != null, true);
}