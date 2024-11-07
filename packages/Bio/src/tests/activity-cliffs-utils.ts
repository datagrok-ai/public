import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {expect} from '@datagrok-libraries/utils/src/test';
import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';
import {BitArrayMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {BYPASS_LARGE_DATA_WARNING} from '@datagrok-libraries/ml/src/functionEditors/consts';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/types';

export async function _testActivityCliffsOpen(df: DG.DataFrame, drMethod: DimReductionMethods,
  seqColName: string, activityColName: string, similarityThr: number, tgtNumberCliffs: number,
  similarityMetric: MmDistanceFunctionsNames | BitArrayMetrics, preprocessingFunction: DG.Func,
): Promise<void> {
  await grok.data.detectSemanticTypes(df);
  const scatterPlot = (await grok.functions.call('Bio:activityCliffs', {
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
  expect(scatterPlot != null, true);

  const cliffsLink = Array.from(scatterPlot!.root.children).find((el) => {
    const classList: string[] = el.className.split(' ');
    return ['ui-btn', 'ui-btn-ok'].every((reqClassName) => classList.includes(reqClassName));
  });
  expect((cliffsLink as HTMLElement).innerText.toLowerCase(), `${tgtNumberCliffs} cliffs`);
}
