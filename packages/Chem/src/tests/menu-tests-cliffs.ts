import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {_package} from '../package-test';
import {createTableView} from './utils';
import {activityCliffs} from '../package';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {before, after, expect, category, test} from '@datagrok-libraries/utils/src/test';
// const {jStat} = require('jstat');


category('top menu activity cliffs', async () => {
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  test('activityCliffsOpen.smiles', async () => {
    await _testActivityCliffsOpen('tests/activity_cliffs_test.csv', 2);
  });

  test('activityCliffsOpen.molV2000', async () => {
    await _testActivityCliffsOpen('tests/spgi-100.csv', 1, 'V2000');
  });

  test('activityCliffsOpen.molV3000', async () => {
    await _testActivityCliffsOpen('v3000_sample.csv', 185, 'V3000');
  });

  test('activityCliffsWithEmptyRows', async () => {
    await _testActivityCliffsOpen('tests/activity_cliffs_empty_rows.csv', 1);
  });

  after(async () => {
    grok.shell.closeAll();
  });
});

async function _testActivityCliffsOpen(dfName: string, numberCliffs: number, ver?: string) {
  const actCliffsTableView = await createTableView(dfName);
  if (ver === 'V3000') actCliffsTableView.dataFrame.rows.removeAt(51, 489);
  await activityCliffs(
    actCliffsTableView.dataFrame,
    actCliffsTableView.dataFrame.getCol(ver ? (ver === 'V2000' ? 'Structure' : 'molecule') : 'smiles'),
    actCliffsTableView.dataFrame.getCol(ver === 'V2000' ? 'Chemical Space X' : 'Activity'),
    80,
    't-SNE');
  let scatterPlot: DG.Viewer | null = null;
  for (const i of actCliffsTableView.viewers) {
    if (i.type == DG.VIEWER.SCATTER_PLOT)
      scatterPlot = i;
  }
  expect(scatterPlot != null, true);
  const cliffsLink = Array.from(scatterPlot!.root.children)
    .filter((it) => it.className === 'ui-btn ui-btn-ok scatter_plot_link cliffs_grid');
  expect((cliffsLink[0] as HTMLElement).innerText.toLowerCase(), `${numberCliffs} cliffs`);
  actCliffsTableView.close();
}
