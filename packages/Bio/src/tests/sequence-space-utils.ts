import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {expect} from '@datagrok-libraries/test/src/test';
import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';
import {BYPASS_LARGE_DATA_WARNING} from '@datagrok-libraries/ml/src/functionEditors/consts';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/types';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';

import {readDataframe} from './utils';

export async function _testSequenceSpaceReturnsResult(
  df: DG.DataFrame, algorithm: DimReductionMethods, colName: string
) {
  // await grok.data.detectSemanticTypes(df);
  const col: DG.Column = df.getCol(colName);
  df.name = 'seqSpaceDf';
  const semType: string = await grok.functions.call('Bio:detectMacromolecule', {col: col});
  if (semType)
    col.semType = semType;

  const preprocessingFunc = DG.Func.find({package: 'Bio', name: 'macromoleculePreprocessingFunction'})[0];
  if (!preprocessingFunc)
    throw new Error('Preprocessing function not found');
  await grok.functions.call('Bio:sequenceSpaceTopMenu', {
    table: df, molecules: df.col(colName)!,
    methodName: algorithm, similarityMetric: MmDistanceFunctionsNames.LEVENSHTEIN,
    plotEmbeddings: true, preprocessingFunction: preprocessingFunc, options: {[BYPASS_LARGE_DATA_WARNING]: true}
  });
  // const sp = await sequenceSpaceTopMenu(df, df.col(colName)!, algorithm, MmDistanceFunctionsNames.LEVENSHTEIN, true,
  //   preprocessingFunc, {[BYPASS_LARGE_DATA_WARNING]: true});
  const tv = grok.shell.tableView(df.name);
  const sp = Array.from(tv?.viewers ?? [])[1];
  expect(sp != null);
}

/** Runs Sequence Space with the dialog's defaults on the first 100 rows of a sample file. */
export async function _testSequenceSpaceOnNotation(fileName: string, colName: string, units: NOTATION): Promise<void> {
  const src = await readDataframe(fileName);
  const df = src.clone(DG.BitSet.create(src.rowCount, (i) => i < 100));
  await grok.data.detectSemanticTypes(df);
  const col = df.getCol(colName);
  expect(col.meta.units, units);
  const tv = grok.shell.addTableView(df);
  await grok.functions.call('Bio:sequenceSpaceTopMenu', {
    table: df, molecules: col, methodName: DimReductionMethods.UMAP,
    similarityMetric: MmDistanceFunctionsNames.HAMMING, plotEmbeddings: true, clusterEmbeddings: true,
    options: {[BYPASS_LARGE_DATA_WARNING]: true},
  });
  for (const name of ['Embed_X_1', 'Embed_Y_1']) {
    const embedCol = df.getCol(name);
    expect(embedCol.stats.missingValueCount, 0, `${name} has missing values`);
    expect(embedCol.stats.uniqueCount >= 10, true, `${name} has ${embedCol.stats.uniqueCount} distinct values`);
  }
  const sp = Array.from(tv.viewers).find((v) => v.type === DG.VIEWER.SCATTER_PLOT);
  expect(sp?.props.xColumnName, 'Embed_X_1');
}
