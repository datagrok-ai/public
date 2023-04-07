import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {delay, expect} from '@datagrok-libraries/utils/src/test';
import {IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';

export async function _testActivityCliffsOpen(df: DG.DataFrame, numberCliffs: number, method: string, colName: string) {
  await grok.data.detectSemanticTypes(df);
  // const scatterPlot = await activityCliffs(
  //   df, df.col(colName)!, df.col('Activity')!,
  //   50, method);
  // const scatterPlot = (await grok.functions.call('Bio:activityCliffs', {
  //   table: df, molecules: df.getCol(colName), activities: df.getCol('Activity'),
  //   similarity: 50, methodName: method
  // })) as DG.Viewer | undefined;
  const libHelper: IMonomerLibHelper = (await grok.functions.call('Bio:getMonomerLibHelper'));
  const k = 11;

  // expect(scatterPlot != null, true);
  //
  // const cliffsLink = Array.from(scatterPlot!.root.children).find((el) => {
  //   const classList: string[] = el.className.split(' ');
  //   return ['ui-btn', 'ui-btn-ok'].every((reqClassName) => classList.includes(reqClassName));
  // });
  // expect((cliffsLink as HTMLElement).innerText.toLowerCase(), `${numberCliffs} cliffs`);
}
