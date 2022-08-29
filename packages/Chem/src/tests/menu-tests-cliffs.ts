import { before, after, expect, category, test } from '@datagrok-libraries/utils/src/test';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import { _package } from '../package-test';
import { createTableView } from './utils';

import { activityCliffs } from '../package';
import { chemSpace } from '../analysis/chem-space';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import { getSimilaritiesMarix, getSimilaritiesMarixFromDistances } from '../utils/similarity-utils';
import { chemSpaceTopMenu } from '../package';
var { jStat } = require('jstat')

category('top menu activity cliffs', async () => {

  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  test('activityCliffsOpen', async () => {
    await _testActivityCliffsOpen('activity_cliffs.csv', 92);
  });

  test('activityCliffsWithEmptyRows', async () => {
    await _testActivityCliffsOpen('tests/activity_cliffs_empty_rows.csv', 91);
  });

  after(async () => {
  });

});

async function _testActivityCliffsOpen(dfName: string, numberCliffs: number) {

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
