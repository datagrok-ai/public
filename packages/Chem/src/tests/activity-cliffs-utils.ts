import * as DG from 'datagrok-api/dg';
import {delay, expect} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import { activityCliffs } from '../package';
import { createTableView, readDataframe } from './utils';


export async function _testActivityCliffsOpen(dfName: string, numberCliffs: number) {

  const actCliffsTableView = await createTableView(dfName);
   const scatterPlot = await activityCliffs(
    actCliffsTableView.dataFrame, 
    actCliffsTableView.dataFrame.col('smiles')!, 
    actCliffsTableView.dataFrame.col('Activity')!, 
     80, 
     't-SNE');

    expect(scatterPlot != null, true);

    const cliffsLink = Array.from(scatterPlot.root.children).filter(it => it.className === 'ui-btn ui-btn-ok');
    expect((cliffsLink[0] as HTMLElement).innerText, `${numberCliffs} cliffs`);
    actCliffsTableView.close();
}