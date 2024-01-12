import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {expect} from '@datagrok-libraries/utils/src/test';
import {activityCliffs} from '../package';
import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';
import {BitArrayMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {BYPASS_LARGE_DATA_WARNING} from '@datagrok-libraries/ml/src/functionEditors/consts';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/types';

export async function _testActivityCliffsOpen(df: DG.DataFrame, drMethod: DimReductionMethods,
  seqColName: string, activityColName: string, similarityThr: number, tgtNumberCliffs: number,
  similarityMetric: MmDistanceFunctionsNames | BitArrayMetrics, preprocessingFunction: DG.Func,
): Promise<void> {
  await grok.data.detectSemanticTypes(df);
  const scatterPlot = await activityCliffs(
    df, df.getCol(seqColName), df.getCol(activityColName),
    similarityThr, drMethod, similarityMetric, preprocessingFunction, {[`${BYPASS_LARGE_DATA_WARNING}`]: true});
  // const scatterPlot = (await grok.functions.call('Bio:activityCliffs', {
  //   table: df, molecules: df.getCol(colName), activities: df.getCol('Activity'),
  //   similarity: 50, methodName: method
  // })) as DG.Viewer | undefined;

  // test scatter plot without activityCliffs passed
  // const scatterPlot = (await df.plot.fromType(DG.VIEWER.SCATTER_PLOT, {})) as DG.Viewer;
  // const libHelper: IMonomerLibHelper = (await grok.functions.call('Bio:getMonomerLibHelper'));
  // const k = 11;

  expect(scatterPlot != null, true);

  const cliffsLink = Array.from(scatterPlot!.root.children).find((el) => {
    const classList: string[] = el.className.split(' ');
    return ['ui-btn', 'ui-btn-ok'].every((reqClassName) => classList.includes(reqClassName));
  });
  expect((cliffsLink as HTMLElement).innerText.toLowerCase(), `${tgtNumberCliffs} cliffs`);
}
