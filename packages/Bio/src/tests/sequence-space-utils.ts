import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {expect} from '@datagrok-libraries/utils/src/test';
import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';
import {BYPASS_LARGE_DATA_WARNING} from '@datagrok-libraries/ml/src/functionEditors/consts';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/types';

export async function _testSequenceSpaceReturnsResult(
  df: DG.DataFrame, algorithm: DimReductionMethods, colName: string
) {
  // await grok.data.detectSemanticTypes(df);
  const col: DG.Column = df.getCol(colName);
  const semType: string = await grok.functions.call('Bio:detectMacromolecule', {col: col});
  if (semType)
    col.semType = semType;

  const preprocessingFunc = DG.Func.find({package: 'Bio', name: 'macromoleculePreprocessingFunction'})[0];
  if (!preprocessingFunc)
    throw new Error('Preprocessing function not found');
  const sp = await grok.functions.call('Bio:sequenceSpaceTopMenu', {
    table: df, molecules: df.col(colName)!,
    methodName: algorithm, similarityMetric: MmDistanceFunctionsNames.LEVENSHTEIN,
    plotEmbeddings: true, preprocessingFunction: preprocessingFunc, options: {[BYPASS_LARGE_DATA_WARNING]: true}
  });
  // const sp = await sequenceSpaceTopMenu(df, df.col(colName)!, algorithm, MmDistanceFunctionsNames.LEVENSHTEIN, true,
  //   preprocessingFunc, {[BYPASS_LARGE_DATA_WARNING]: true});
  expect(sp != null, true);
}
