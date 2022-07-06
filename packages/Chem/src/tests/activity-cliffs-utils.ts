import * as DG from 'datagrok-api/dg';
import {delay, expect} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import { chemSpace, getEmbeddingColsNames } from '../analysis/chem-space';
import { getActivityCliffs } from '@datagrok-libraries/ml/src/viewers/activity-cliffs';
import { chemGetSimilarities } from '../chem-searches';
import { drawTooltip } from '../analysis/activity-cliffs';


export async function _testActivityCliffsOpen(df: DG.DataFrame) {
    const axesNames = getEmbeddingColsNames(df);
    const options = {
     'SPE': {cycles: 2000, lambda: 1.0, dlambda: 0.0005},
   };
   const scatterPlot = await getActivityCliffs(
     df, 
     df.col('smiles')!, 
     axesNames, 
     df.col('Activity')!, 
     80, 
     'Tanimoto',
     't-SNE',
     DG.SEMTYPE.MOLECULE,
     'smiles',
     chemSpace,
     chemGetSimilarities,
     drawTooltip);

    expect(scatterPlot != null, true);

    const cliffsLink = Array.from(scatterPlot.root.children).filter(it => it.className === 'ui-btn ui-btn-ok');
    expect((cliffsLink[0] as HTMLElement).innerText, '92 cliffs');
}