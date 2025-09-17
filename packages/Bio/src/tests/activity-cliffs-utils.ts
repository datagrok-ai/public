import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as api from '../package-api';

import {awaitCheck, expect} from '@datagrok-libraries/utils/src/test';
import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';
import {BitArrayMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {BYPASS_LARGE_DATA_WARNING} from '@datagrok-libraries/ml/src/functionEditors/consts';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/types';

export async function _testActivityCliffsOpen(df: DG.DataFrame, drMethod: DimReductionMethods,
  seqColName: string, activityColName: string, similarityThr: number, tgtNumberCliffs: number,
  similarityMetric: MmDistanceFunctionsNames | BitArrayMetrics, preprocessingFunction: DG.Func,
): Promise<void> {
  await grok.data.detectSemanticTypes(df);
  await api.funcs.activityCliffs(
    df, 
    df.getCol(seqColName), 
    df.getCol(activityColName), 
    similarityThr, 
    drMethod, 
    similarityMetric, 
    preprocessingFunction, 
    {[`${BYPASS_LARGE_DATA_WARNING}`]: true}, 
    false) as DG.Viewer | undefined);
  const scatterPlot = Array.from(grok.shell.tv.viewers)[1];
  expect(scatterPlot?.type === DG.VIEWER.SCATTER_PLOT, true);

  await awaitCheck(() => {
    const link = Array.from(scatterPlot!.root.getElementsByClassName('scatter_plot_link'));
    if (link.length)
      return (link[0] as HTMLElement).innerText.toLowerCase() === `${tgtNumberCliffs} cliffs`;
    return true;
  }, 'incorrect cliffs link', 3000);
}
