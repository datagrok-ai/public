import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {expect} from '@datagrok-libraries/utils/src/test';
import {activityCliffs, BYPASS_LARGE_DATA_WARNING} from '../package';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/reduce-dimensionality';

export async function _testActivityCliffsOpen(df: DG.DataFrame, numberCliffs: number, method: DimReductionMethods,
  colName: string) {
  await grok.data.detectSemanticTypes(df);
  const scatterPlot = await activityCliffs(
    df, df.getCol(colName), df.getCol('activity'),
    90, method, {[`${BYPASS_LARGE_DATA_WARNING}`]: true});
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
  expect((cliffsLink as HTMLElement).innerText.toLowerCase(), `${numberCliffs} cliffs`);
}
