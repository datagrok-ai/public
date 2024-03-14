import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {_package} from '../package-test';
import {readDataframe} from './utils';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {before, after, expect, category, test, awaitCheck} from '@datagrok-libraries/utils/src/test';
import {BitArrayMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {getActivityCliffs} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';
import {createPropPanelElement, createTooltipElement} from '../analysis/activity-cliffs';
import {MALFORMED_DATA_WARNING_CLASS} from '../constants';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/types';
// const {jStat} = require('jstat');


category('top menu activity cliffs', async () => {
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  test('activityCliffsOpen.smiles', async () => {
    const df = DG.Test.isInBenchmark ?
      await grok.data.files.openTable('System:AppData/Chem/tests/smiles_10K_with_activities.csv') :
      await readDataframe('tests/activity_cliffs_test.csv');
    await _testActivityCliffsOpen(df, 'smiles', 'Activity', DG.Test.isInBenchmark ? 78 : 2);
  });

  test('activityCliffsOpen.molV2000', async () => {
    await _testActivityCliffsOpen(await readDataframe('tests/spgi-100.csv'), 'Structure', 'Chemical Space X', 1);
  });

  test('activityCliffsOpen.molV3000', async () => {
    await _testActivityCliffsOpen(await readDataframe('v3000_sample.csv'), 'molecule', 'Activity', 185);
  });

  test('activityCliffs.emptyValues', async () => {
    await _testActivityCliffsOpen(await readDataframe('tests/activity_cliffs_empty_rows.csv'),
      'smiles', 'Activity', 1);
  });

  test('activityCliffs.malformedData', async () => {
    DG.Balloon.closeAll();
    await _testActivityCliffsOpen(await readDataframe('tests/Test_smiles_malformed.csv'),
      'canonical_smiles', 'FractionCSP3', 24);
    try {
      await awaitCheck(() => document.querySelector(`.${MALFORMED_DATA_WARNING_CLASS}`)?.innerHTML ===
        '2 molecules with indexes 31,41 are possibly malformed and are not included in analysis',
      'cannot find warning balloon', 5000);
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

async function _testActivityCliffsOpen(df: DG.DataFrame, molCol: string, activityCol: string, numberCliffs: number) {
  await grok.data.detectSemanticTypes(df);
  const actCliffsTableView = grok.shell.addTableView(df);
  await awaitCheck(() => actCliffsTableView.name?.toLowerCase() === df.name?.toLowerCase(),
    'Activity cliffs table view hasn\'t been created', 1000);
  if (molCol === 'molecule') actCliffsTableView.dataFrame.rows.removeAt(51, 489);
  const encodingFunc = DG.Func.find({name: 'getFingerprints', package: 'Chem'})[0];
  await getActivityCliffs(df, df.col(molCol)!,
    ['Embed_X_1', 'Embed_Y_1'], 'Activity cliffs', actCliffsTableView.dataFrame.getCol(activityCol),
    80, BitArrayMetricsNames.Tanimoto, DimReductionMethods.UMAP, {}, DG.SEMTYPE.MOLECULE,
    {'units': df.col(molCol)!.tags['units']}, encodingFunc, createTooltipElement, createPropPanelElement);
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
