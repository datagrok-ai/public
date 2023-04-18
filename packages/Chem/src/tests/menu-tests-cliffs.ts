import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {_package} from '../package-test';
import {createTableView} from './utils';
import {activityCliffs} from '../package';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {before, after, expect, category, test, awaitCheck} from '@datagrok-libraries/utils/src/test';
// const {jStat} = require('jstat');


category('top menu activity cliffs', async () => {
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  test('activityCliffsOpen.smiles', async () => {
    await _testActivityCliffsOpen('tests/activity_cliffs_test.csv', 'smiles', 'Activity', 2);
  });

  test('activityCliffsOpen.molV2000', async () => {
    await _testActivityCliffsOpen('tests/spgi-100.csv', 'Structure', 'Chemical Space X', 1);
  });

  test('activityCliffsOpen.molV3000', async () => {
    await _testActivityCliffsOpen('v3000_sample.csv', 'molecule', 'Activity', 185);
  });

  test('activityCliffs.emptyValues', async () => {
    await _testActivityCliffsOpen('tests/activity_cliffs_empty_rows.csv', 'smiles', 'Activity', 1);
  });

  test('activityCliffs.malformedData', async () => {
    DG.Balloon.closeAll();
    await _testActivityCliffsOpen('tests/Test_smiles_malformed.csv', 'canonical_smiles', 'FractionCSP3', 24);
    try {
      await awaitCheck(() => document.querySelector('.d4-balloon-content')?.innerHTML ===
        '2 molecules with indexes 2,9 are possibly malformed and are not included in analysis', 'cannot find warning balloon', 1000);
    } finally {
      grok.shell.closeAll();
      DG.Balloon.closeAll();
    }
  });

  after(async () => {
    grok.shell.closeAll();
    DG.Balloon.closeAll();
  });
});

async function _testActivityCliffsOpen(dfName: string, molCol: string, activityCol: string, numberCliffs: number) {
  const actCliffsTableView = await createTableView(dfName);
  if (molCol === 'molecule') actCliffsTableView.dataFrame.rows.removeAt(51, 489);
  await activityCliffs(
    actCliffsTableView.dataFrame,
    actCliffsTableView.dataFrame.getCol(molCol),
    actCliffsTableView.dataFrame.getCol(activityCol),
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
