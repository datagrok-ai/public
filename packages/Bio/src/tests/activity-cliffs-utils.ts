import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {waitUntil, expect} from '@datagrok-libraries/test/src/test';
import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';
import {BitArrayMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {BYPASS_LARGE_DATA_WARNING} from '@datagrok-libraries/ml/src/functionEditors/consts';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/types';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';

import {readDataframe} from './utils';

export async function _testActivityCliffsOpen(df: DG.DataFrame, drMethod: DimReductionMethods,
  seqColName: string, activityColName: string, similarityThr: number, tgtNumberCliffs: number,
  similarityMetric: MmDistanceFunctionsNames | BitArrayMetrics, preprocessingFunction: DG.Func,
): Promise<void> {
  await grok.data.detectSemanticTypes(df);
  (await grok.functions.call('Bio:activityCliffs', {
    table: df,
    molecules: df.getCol(seqColName),
    activities: df.getCol(activityColName),
    similarity: similarityThr,
    methodName: drMethod,
    similarityMetric: similarityMetric,
    preprocessingFunction: preprocessingFunction,
    options: {[`${BYPASS_LARGE_DATA_WARNING}`]: true},
    demo: false,
  })) as DG.Viewer | undefined;
  const scatterPlot = Array.from(grok.shell.tv.viewers)[1];
  expect(scatterPlot?.type === DG.VIEWER.SCATTER_PLOT, true);

  await waitUntil(() => {
    const link = Array.from(scatterPlot!.root.getElementsByClassName('scatter_plot_link'));
    if (link.length)
      return (link[0] as HTMLElement).innerText.toLowerCase() === `${tgtNumberCliffs} cliffs`;
    return true;
  }, 'incorrect cliffs link', 3000);
}

/** Runs Activity Cliffs with the dialog's defaults on the first 100 rows of a sample file. */
export async function _testActivityCliffsOnNotation(fileName: string, seqColName: string, units: NOTATION,
  preprocessingFunction: DG.Func,
): Promise<void> {
  const src = await readDataframe(fileName);
  const df = src.clone(DG.BitSet.create(src.rowCount, (i) => i < 100));
  await grok.data.detectSemanticTypes(df);
  expect(df.getCol(seqColName).meta.units, units);
  const tv = grok.shell.addTableView(df);
  await grok.functions.call('Bio:activityCliffs', {
    table: df, molecules: df.getCol(seqColName), activities: df.getCol('Activity'), similarity: 80,
    methodName: DimReductionMethods.UMAP, similarityMetric: MmDistanceFunctionsNames.HAMMING,
    preprocessingFunction: preprocessingFunction, options: {[BYPASS_LARGE_DATA_WARNING]: true}, demo: false,
  });
  expect(df.getCol('Embed_X_1').stats.missingValueCount, 0, 'Embed_X_1 has missing values');
  await waitUntil(() => {
    const sp = Array.from(tv.viewers).find((v) => v.type === DG.VIEWER.SCATTER_PLOT);
    const link = sp?.root.querySelector('.scatter_plot_link') as HTMLElement | null;
    return parseInt(link?.innerText ?? '') >= 1;
  }, `no cliffs found in ${seqColName}`, 10000);
  expect(df.columns.names().some((name) => /sali/i.test(name)), true, 'no SALI column');
}
