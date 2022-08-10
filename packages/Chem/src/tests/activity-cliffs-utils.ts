import * as DG from 'datagrok-api/dg';
import {delay, expect} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import { activityCliffs } from '../package';


export async function _testActivityCliffsOpen(df: DG.DataFrame, numberCliffs: number) {
   const scatterPlot = await activityCliffs(
     df, 
     df.col('smiles')!, 
     df.col('Activity')!, 
     80, 
     't-SNE');

    expect(scatterPlot != null, true);

    const cliffsLink = Array.from(scatterPlot.root.children).filter(it => it.className === 'ui-btn ui-btn-ok');
    expect((cliffsLink[0] as HTMLElement).innerText, `${numberCliffs} cliffs`);
}