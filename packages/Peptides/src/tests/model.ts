import * as DG from 'datagrok-api/dg';

import {category, test, before, expect, expectFloat} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import {PeptidesModel} from '../model';
import {startAnalysis} from '../widgets/peptides';
import {scaleActivity} from '../utils/misc';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {COLUMNS_NAMES, SCALING_METHODS} from '../utils/constants';

category('Model: Settings', () => {
  let df: DG.DataFrame;
  let model: PeptidesModel;
  let activityCol: DG.Column<number>;
  let sequenceCol: DG.Column<string>;
  let clusterCol: DG.Column<any>;
  let scaledActivityCol: DG.Column<number>;

  const mutationCliffsDefaultParams = {maxMutations: 1, minActivityDelta: 0};
  const mutationCliffsTestParams = {maxMutations: 2, minActivityDelta: 0.5};

  before(async () => {
    df = DG.DataFrame.fromCsv(await _package.files.readAsText('tests/HELM_small.csv'));
    activityCol = df.getCol('activity');
    sequenceCol = df.getCol('sequence');
    sequenceCol.semType = DG.SEMTYPE.MACROMOLECULE;
    sequenceCol.setTag(DG.TAGS.UNITS, NOTATION.HELM);
    scaledActivityCol = scaleActivity(activityCol, SCALING_METHODS.NONE);
    clusterCol = df.getCol('cluster');
    const tempModel = await startAnalysis(activityCol, sequenceCol, clusterCol, df, scaledActivityCol,
      SCALING_METHODS.NONE);
    if (tempModel === null)
      throw new Error('Model is null');
    model = tempModel;
  });

  test('Activity scaling', async () => {
    const getError = (row: number, method: SCALING_METHODS): string =>
      `Activity mismatch at row ${row} for scaling method '${method}'`;
    const tolerance = 0.0001;
    const origActivityData =
      model.df.getCol(model.settings.activityColumnName!).getRawData();
    const scaledActivity = model.df.getCol(COLUMNS_NAMES.ACTIVITY_SCALED);
    const dfLen = model.df.rowCount;

    // Check initial 'none' scaling
    let scaledActivityData = scaledActivity.getRawData();
    for (let i = 0; i < dfLen; i++)
      expectFloat(scaledActivityData[i], origActivityData[i], tolerance, getError(i, SCALING_METHODS.NONE));

    // Check 'lg' scaling
    model.settings = {scaling: SCALING_METHODS.LG};
    scaledActivityData = scaledActivity.getRawData();
    for (let i = 0; i < dfLen; i++)
      expectFloat(scaledActivityData[i], Math.log10(origActivityData[i]), tolerance, getError(i, SCALING_METHODS.LG));

    // Check '-lg' scaling
    model.settings = {scaling: SCALING_METHODS.MINUS_LG};
    scaledActivityData = scaledActivity.getRawData();
    for (let i = 0; i < dfLen; i++) {
      expectFloat(scaledActivityData[i], -Math.log10(origActivityData[i]), tolerance,
        getError(i, SCALING_METHODS.MINUS_LG));
    }

    // Check 'none' scaling
    model.settings = {scaling: SCALING_METHODS.NONE};
    scaledActivityData = scaledActivity.getRawData();
    for (let i = 0; i < dfLen; i++)
      expectFloat(scaledActivityData[i], origActivityData[i], tolerance, getError(i, SCALING_METHODS.NONE));
  });

  test('Bidirectional analysis', async () => {
    // Check that bidirectional analysis is disabled by default
    expect(model.settings.isBidirectional, false, 'Bidirectional analysis is enabled by default');

    // Check that bidirectional analysis can be enabled
    model.settings = {isBidirectional: true};
    expect(model.settings.isBidirectional, true, 'Bidirectional analysis is disabled after enabling');

    // Check that bidirectional analysis can be disabled
    model.settings = {isBidirectional: false};
    expect(model.settings.isBidirectional, false, 'Bidirectional analysis is enabled after disabling');
  });

  test('Mutation Cliffs', async () => {
    // Check default mutation cliffs parameters
    expect(model.settings.maxMutations, mutationCliffsDefaultParams.maxMutations, `Max mutations mismatch: expected ` +
      `${mutationCliffsDefaultParams.maxMutations}, actual ${model.settings.maxMutations}`);
    expect(model.settings.minActivityDelta, mutationCliffsDefaultParams.minActivityDelta, `Min activity delta ` +
      `mismatch: expected ${mutationCliffsDefaultParams.minActivityDelta}, actual ${model.settings.minActivityDelta}`);

    // Check test mutation cliffs parameters
    model.settings = {maxMutations: mutationCliffsTestParams.maxMutations,
      minActivityDelta: mutationCliffsTestParams.minActivityDelta};
    expect(model.settings.maxMutations, mutationCliffsTestParams.maxMutations, `Max mutations mismatch: expected ` +
      `${mutationCliffsTestParams.maxMutations}, actual ${model.settings.maxMutations}`);
    expect(model.settings.minActivityDelta, mutationCliffsTestParams.minActivityDelta, `Min activity delta ` +
      `mismatch: expected ${mutationCliffsTestParams.minActivityDelta}, actual ${model.settings.minActivityDelta}`);
  });

  test('Include columns', async () => {

  }, {skipReason: 'Not implemented yet'});

  test('Dendrogram', async () => {

  }, {skipReason: 'Not implemented yet'});
});
