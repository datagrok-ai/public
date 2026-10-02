import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {after, before, category, expectArray, test} from '@datagrok-libraries/test/src/test';

import {readDataframe} from './utils';
import {_testActivityCliffsOpen} from './activity-cliffs-utils';

import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';
import {BitArrayMetricsNames, NumberMetricsNames} from '@datagrok-libraries/ml/src/typed-metrics';
import {BYPASS_LARGE_DATA_WARNING} from '@datagrok-libraries/ml/src/functionEditors/consts';
import {getEmbeddingColsNames, multiColReduceDimensionality, releaseEmbeddingColsNames}
  from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/reduce-dimensionality';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/types/monomer-library';
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
  }, {benchmark: true, skipReason: 'Fails'});

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

  test('embeddingColsNamesSkipTakenNames', async () => {
    const df = DG.DataFrame.create(3);
    df.columns.addNewFloat('Embed_X_2');
    df.columns.addNewFloat('Embed_Y_2');
    expectArray(getEmbeddingColsNames(df), ['Embed_X_3', 'Embed_Y_3']);
  });

  test('embeddingColsNamesOfOverlappingRuns', async () => {
    const df = DG.DataFrame.create(3);
    const first = getEmbeddingColsNames(df);
    expectArray(getEmbeddingColsNames(df), ['Embed_X_2', 'Embed_Y_2']);
    releaseEmbeddingColsNames(df, first);
    expectArray(getEmbeddingColsNames(df), ['Embed_X_1', 'Embed_Y_1']);
  });

  test('embeddingColsNamesOfFailedRun', async () => {
    const df = DG.DataFrame.create(3);
    const col = df.columns.addNewFloat('v');
    const names = getEmbeddingColsNames(df);
    df.columns.addNewFloat(names[1]);
    await multiColReduceDimensionality(df, [col], DimReductionMethods.UMAP, [NumberMetricsNames.Difference], [1],
      [null], 'MANHATTAN', false, false, {preprocessingFuncArgs: [{}]}, {embedColsNames: names});
    for (const name of names)
      df.columns.remove(name);
    expectArray(getEmbeddingColsNames(df), names);
  });

  test('activityCliffsAlongsideSequenceSpace', async () => {
    const df = await readDataframe('tests/100_3_clustests.csv');
    grok.shell.addTableView(df);
    await grok.data.detectSemanticTypes(df);
    const options = {[BYPASS_LARGE_DATA_WARNING]: true};
    await Promise.all([
      grok.functions.call('Bio:activityCliffs', {table: df, molecules: df.getCol('sequence'),
        activities: df.getCol('Activity'), similarity: 90, methodName: DimReductionMethods.UMAP,
        similarityMetric: MmDistanceFunctionsNames.LEVENSHTEIN, preprocessingFunction: seqEncodingFunc, options}),
      grok.functions.call('Bio:sequenceSpaceTopMenu', {table: df, molecules: df.getCol('sequence'),
        methodName: DimReductionMethods.UMAP, similarityMetric: MmDistanceFunctionsNames.LEVENSHTEIN,
        plotEmbeddings: true, preprocessingFunction: seqEncodingFunc, options}),
    ]);
    expectArray(df.columns.names().filter((it) => it.startsWith('Embed_')).sort(),
      ['Embed_X_1', 'Embed_X_2', 'Embed_Y_1', 'Embed_Y_2']);
  });
});
