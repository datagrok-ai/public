import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test} from '@datagrok-libraries/utils/src/test';

import {readDataframe} from './utils';
import {_testActivityCliffsOpen} from './activity-cliffs-utils';

import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';
import {BitArrayMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {
  getUserLibSettings, setUserLibSettings, setUserLibSettingsForTests
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/types';
import {getHelmHelper, IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';

import {_package} from '../package-test';


category('activityCliffs', async () => {
  let viewList: DG.ViewBase[] = [];
  let dfList: DG.DataFrame[] = [];

  let helmHelper: IHelmHelper;
  let monomerLibHelper: IMonomerLibHelper;
  /** Backup actual user's monomer libraries settings */
  let userLibSettings: UserLibSettings;
  const seqEncodingFunc = DG.Func.find({name: 'macromoleculePreprocessingFunction', package: 'Bio'})[0];
  const helmEncodingFunc = DG.Func.find({name: 'helmPreprocessingFunction', package: 'Bio'})[0];
  before(async () => {
    helmHelper = await getHelmHelper(); // init Helm package
    monomerLibHelper = await getMonomerLibHelper();
    userLibSettings = await getUserLibSettings();

    // Test 'helm' requires default monomer library loaded
    await setUserLibSettingsForTests();
    await monomerLibHelper.loadLibraries(true); // load default libraries

    viewList = [];
    dfList = [];
  });

  after(async () => {
    // for (const df of dfList) grok.shell.closeTable(df);
    // for (const view of viewList) view.close();

    // UserDataStorage.put() replaces existing data
    await setUserLibSettings(userLibSettings);
    await monomerLibHelper.loadLibraries(true); // load user settings libraries
  });

  test('activityCliffsOpens', async () => {
    const actCliffsDf = await readDataframe(
      DG.Test.isInBenchmark ? 'test/peptides_motif-with-random_10000.csv' : 'tests/100_3_clustests.csv',
    );
    dfList.push(actCliffsDf);
    const actCliffsTableView = grok.shell.addTableView(actCliffsDf);
    viewList.push(actCliffsTableView);
    const cliffsNum = DG.Test.isInBenchmark ? 6 : 3;

    await _testActivityCliffsOpen(actCliffsDf, DimReductionMethods.UMAP,
      'sequence', 'Activity', 90, cliffsNum, MmDistanceFunctionsNames.LEVENSHTEIN, seqEncodingFunc);
  });

  test('activityCliffsWithEmptyRows', async () => {
    const actCliffsDfWithEmptyRows = await readDataframe('tests/100_3_clustests_empty_vals.csv');
    dfList.push(actCliffsDfWithEmptyRows);
    const actCliffsTableViewWithEmptyRows = grok.shell.addTableView(actCliffsDfWithEmptyRows);
    viewList.push(actCliffsTableViewWithEmptyRows);

    await _testActivityCliffsOpen(actCliffsDfWithEmptyRows, DimReductionMethods.UMAP,
      'sequence', 'Activity', 90, 3, MmDistanceFunctionsNames.LEVENSHTEIN, seqEncodingFunc);
  });

  test('Helm', async () => {
    const df = await _package.files.readCsv('samples/HELM_50.csv');
    const _view = grok.shell.addTableView(df);

    await _testActivityCliffsOpen(df, DimReductionMethods.UMAP,
      'HELM', 'Activity', 65, 20, BitArrayMetricsNames.Tanimoto, helmEncodingFunc);
  });
});
