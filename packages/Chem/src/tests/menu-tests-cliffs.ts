import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {_package} from '../package-test';
import {readDataframe} from './utils';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {before, after, expect, category, test, awaitCheck, delay} from '@datagrok-libraries/utils/src/test';
import {BitArrayMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {getActivityCliffs} from '@datagrok-libraries/ml/src/viewers/activity-cliffs';
import {createPropPanelElement, createTooltipElement} from '../analysis/activity-cliffs';
import {MALFORMED_DATA_WARNING_CLASS} from '../constants';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/types';
import { activityCliffs } from '../package';


category('top menu activity cliffs', async () => {
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  test('activityCliffsOpenAndLayoutApply.smiles', async () => {
    const df = DG.Test.isInBenchmark ?
      await grok.data.files.openTable('System:AppData/Chem/tests/smiles_10K_with_activities.csv') :
      await readDataframe('tests/activity_cliffs_test.csv');
    await _testActivityCliffsOpen(df, 'smiles', 'Activity', DG.Test.isInBenchmark ? 78 : 2, !DG.Test.isInBenchmark);
  }, {timeout: 20000, benchmark: true});

  test('activityCliffsOpen.molV2000', async () => {
    await _testActivityCliffsOpen(await readDataframe('tests/spgi-100.csv'), 'Structure', 'Chemical Space X', 1);
  }, {timeout: 20000});

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

  test('activityCliffs_layout', async () => {
    await _testActivityCliffsOpen(await readDataframe('tests/spgi-100.csv'), 'Structure', 'Chemical Space X', 1);
  });

  after(async () => {
    grok.shell.closeAll();
    DG.Balloon.closeAll();
  });
});

async function _testActivityCliffsOpen(df: DG.DataFrame, molCol: string, activityCol: string,
  numberCliffs: number, layout?: boolean) {
  await grok.data.detectSemanticTypes(df);
  const actCliffsTableView = grok.shell.addTableView(df);
  await awaitCheck(() => actCliffsTableView.name?.toLowerCase() === df.name?.toLowerCase(),
    'Activity cliffs table view hasn\'t been created', 1000);
  if (molCol === 'molecule') actCliffsTableView.dataFrame.rows.removeAt(51, 489);
  const encodingFunc = DG.Func.find({name: 'getFingerprints', package: 'Chem'})[0];
  await activityCliffs(df, df.col(molCol)!, actCliffsTableView.dataFrame.getCol(activityCol), 80,
    DimReductionMethods.UMAP, BitArrayMetricsNames.Tanimoto, encodingFunc);
  let scatterPlot: DG.Viewer | null = null;
  for (const i of actCliffsTableView.viewers) {
    if (i.type == DG.VIEWER.SCATTER_PLOT)
      scatterPlot = i;
  }
  expect(scatterPlot != null, true);
  const checkScatterPlotInitialized = () => {
    const cliffsLink = Array.from(scatterPlot!.root.children)
    .filter((it) => it.className === 'ui-btn ui-btn-ok scatter_plot_link cliffs_grid');
    return !!cliffsLink.length && ((cliffsLink[0] as HTMLElement).innerText.toLowerCase() === `${numberCliffs} cliffs`)
  }
  await awaitCheck(() => {
    return checkScatterPlotInitialized();
  }, 'activity cliffs failed', 1000);
  if (layout) {
    const layout = actCliffsTableView.saveLayout();
    scatterPlot?.close();
    await delay(100);
    actCliffsTableView.loadLayout(layout);
    await delay(100); //to add GridViewer
    //waiting for layout to be applied
    await awaitCheck(() => {
      return checkScatterPlotInitialized();
    }, 'activity cliffs failed', 5000);
  }

  actCliffsTableView.close();
}
