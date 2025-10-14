import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {after, awaitCheck, before, category, test} from '@datagrok-libraries/utils/src/test';

import {BitArrayMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {
  getUserLibSettings, setUserLibSettings
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/types';
import {getHelmHelper, IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {expect} from '@datagrok-libraries/utils/src/test';
import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';
import {BitArrayMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {BYPASS_LARGE_DATA_WARNING} from '@datagrok-libraries/ml/src/functionEditors/consts';

import {_package} from '../package-test';


category('activityCliffs', async () => {
  let helmHelper: IHelmHelper;
  let monomerLibHelper: IMonomerLibHelper;
  /** Backup actual user's monomer libraries settings */
  let userLibSettings: UserLibSettings;
  const helmEncodingFunc = DG.Func.find({name: 'helmPreprocessingFunction', package: 'Bio'})[0];
  before(async () => {
    helmHelper = await getHelmHelper(); // init Helm package
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();

    // Test 'helm' requires default monomer library loaded
    await monomerLibHelper.loadMonomerLibForTests();
  });

  after(async () => {
    // UserDataStorage.put() replaces existing data
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadMonomerLib(true); // load user settings libraries
  });


  test('Helm', async () => {
    const df = await _package.files.readCsv('samples/HELM_50.csv');
    const _view = grok.shell.addTableView(df);

    await _testActivityCliffsOpen(df, DimReductionMethods.UMAP,
      'HELM', 'Activity', 65, 20, BitArrayMetricsNames.Tanimoto, helmEncodingFunc);
  });
});


async function _testActivityCliffsOpen(df: DG.DataFrame, drMethod: DimReductionMethods,
  seqColName: string, activityColName: string, similarityThr: number, tgtNumberCliffs: number,
  similarityMetric: MmDistanceFunctionsNames | BitArrayMetrics, preprocessingFunction: DG.Func,
): Promise<void> {
  await grok.data.detectSemanticTypes(df);
  const scatterPlot = (await grok.functions.call('Bio:activityCliffs', {
    table: df,
    molecules: df.getCol(seqColName),
    activities: df.getCol(activityColName),
    similarity: similarityThr,
    methodName: drMethod,
    similarityMetric: similarityMetric,
    preprocessingFunction: preprocessingFunction,
    options: {[`${BYPASS_LARGE_DATA_WARNING}`]: true},
    demo: false,
  })) as DG.Viewer | undefined;
  expect(scatterPlot != null, true);

  const checkScatterPlotInitialized = () => {
    const cliffsLink = scatterPlot!.root.getElementsByClassName('scatter_plot_link');
    return cliffsLink.length > 0 &&
      ((cliffsLink[0] as HTMLElement).innerText.toLowerCase() === `${tgtNumberCliffs} cliffs`);
  };
  await awaitCheck(() => {
    return checkScatterPlotInitialized();
  }, 'activity cliffs failed');
}
