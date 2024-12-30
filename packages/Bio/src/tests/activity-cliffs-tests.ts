import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {after, before, category, test} from '@datagrok-libraries/utils/src/test';

import {readDataframe} from './utils';
import {_testActivityCliffsOpen} from './activity-cliffs-utils';

import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';
import {BitArrayMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {
  getUserLibSettings, setUserLibSettings
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/types';
import {getHelmHelper, IHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';

import {_package} from '../package-test';


category('activityCliffs', async () => {
  let helmHelper: IHelmHelper;
  let monomerLibHelper: IMonomerLibHelper;
  /** Backup actual user's monomer libraries settings */
  let userLibSettings: UserLibSettings;
  const seqEncodingFunc = DG.Func.find({name: 'macromoleculePreprocessingFunction', package: 'Bio'})[0];
  const helmEncodingFunc = DG.Func.find({name: 'helmPreprocessingFunction', package: 'Bio'})[0];
  before(async () => {
    const helmPackInstalled = DG.Func.find({package: 'Helm', name: 'getHelmHelper'}).length;
    if (helmPackInstalled)
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

  test('activityCliffsOpens', async () => {
    const testData = !DG.Test.isInBenchmark ?
      {fileName: 'tests/100_3_clustests.csv', tgt: {cliffCount: 3}} :
      {fileName: 'tests/peptides_with_random_motif_1600.csv', tgt: {cliffCount: 64}};
    const actCliffsDf = await readDataframe(testData.fileName);
    const actCliffsTableView = grok.shell.addTableView(actCliffsDf);

    await _testActivityCliffsOpen(actCliffsDf, DimReductionMethods.UMAP,
      'sequence', 'Activity', 90, testData.tgt.cliffCount, MmDistanceFunctionsNames.LEVENSHTEIN, seqEncodingFunc);
  }, {benchmark: true});

  test('activityCliffsWithEmptyRows', async () => {
    const actCliffsDfWithEmptyRows = await readDataframe('tests/100_3_clustests_empty_vals.csv');
    const actCliffsTableViewWithEmptyRows = grok.shell.addTableView(actCliffsDfWithEmptyRows);

    await _testActivityCliffsOpen(actCliffsDfWithEmptyRows, DimReductionMethods.UMAP,
      'sequence', 'Activity', 90, 3, MmDistanceFunctionsNames.LEVENSHTEIN, seqEncodingFunc);
  });

  test('Helm', async () => {
    const helmPackInstalled = DG.Func.find({package: 'Helm', name: 'getHelmHelper'}).length;
    if (helmPackInstalled) {
      const df = await _package.files.readCsv('samples/HELM_50.csv');
      const _view = grok.shell.addTableView(df);
  
      await _testActivityCliffsOpen(df, DimReductionMethods.UMAP,
        'HELM', 'Activity', 65, 20, BitArrayMetricsNames.Tanimoto, helmEncodingFunc);
    }
  });
});
