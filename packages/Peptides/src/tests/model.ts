import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, test, before, expect, expectFloat, awaitCheck} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import {PeptidesModel, VIEWER_TYPE, getAggregatedColName} from '../model';
import {startAnalysis} from '../widgets/peptides';
import {scaleActivity} from '../utils/misc';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {COLUMNS_NAMES, SCALING_METHODS} from '../utils/constants';
import {LogoSummaryTable} from '../viewers/logo-summary';
import {TEST_COLUMN_NAMES} from './utils';

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
    activityCol = df.getCol(TEST_COLUMN_NAMES.ACTIVITY);
    sequenceCol = df.getCol(TEST_COLUMN_NAMES.SEQUENCE);
    sequenceCol.semType = DG.SEMTYPE.MACROMOLECULE;
    sequenceCol.setTag(DG.TAGS.UNITS, NOTATION.HELM);
    scaledActivityCol = scaleActivity(activityCol, SCALING_METHODS.NONE);
    clusterCol = df.getCol(TEST_COLUMN_NAMES.CLUSTER);
    const tempModel = await startAnalysis(activityCol, sequenceCol, clusterCol, df, scaledActivityCol,
      SCALING_METHODS.NONE);
    if (tempModel === null)
      throw new Error('Model is null');
    model = tempModel;
    let overlayInit = false;
    model._analysisView!.grid.onAfterDrawOverlay.subscribe(() => overlayInit = true);

    // Ensure grid finished initializing to prevent Unhandled exceptions
    let accrodionInit = false;
    grok.events.onAccordionConstructed.subscribe((_) => accrodionInit = true);
    await awaitCheck(() => model!.df.currentRowIdx === 0, 'Grid cell never finished initializing', 2000);
    await awaitCheck(() => grok.shell.o instanceof DG.Column, 'Shell object never changed', 2000);
    await awaitCheck(() => accrodionInit, 'Accordion never finished initializing', 2000);
    await awaitCheck(() => overlayInit, 'Overlay never finished initializing', 2000);
  });

  test('Activity scaling', async () => {
    const getError = (row: number, method: SCALING_METHODS): string =>
      `Activity mismatch at row ${row} for scaling method '${method}'`;
    const tolerance = 0.0001;
    const origActivityData = model.df.getCol(model.settings.activityColumnName!).getRawData();
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
    const columnName = 'rank';
    const testColumns = {[columnName]: DG.AGG.AVG};

    // Include column
    model.settings = {columns: testColumns};

    expect(Object.keys(model.settings.columns!)[0], Object.keys(testColumns)[0], 'Expected to include column ' +
      `'${Object.keys(testColumns)[0]}' but '${Object.keys(model.settings.columns!)[0]}' is included instead`);
    expect(model.settings.columns![columnName], testColumns[columnName], `Expected to aggregate column ` +
      `'${Object.keys(testColumns)[0]}' with '${testColumns[columnName]}' but aggregated with ` +
      `'${model.settings.columns![columnName]}' instead`);

    const lstViewer = model.findViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE) as LogoSummaryTable;
    const aggColName = getAggregatedColName(testColumns[columnName], columnName);
    expect(lstViewer.viewerGrid.col(aggColName) !== null, true, `Expected to include column '${columnName}' in ` +
      `${VIEWER_TYPE.LOGO_SUMMARY_TABLE} but it is absent`);

    // Remove column
    model.settings = {columns: {}};
    expect(Object.keys(model.settings.columns!).length, 0,
      `Expected to remove all column aggregations but columns {${Object.keys(model.settings.columns!).join(' & ')}} ` +
      `are still included`);

    expect(lstViewer.viewerGrid.col(aggColName) === null, true, `Expected to remove column '${columnName}' from ` +
      `${VIEWER_TYPE.LOGO_SUMMARY_TABLE} but it is still present`);
  });

  test('Dendrogram', async () => {
    // Enable dendrogram
    model.settings = {showDendrogram: true};
    expect(model.settings.showDendrogram, true, 'Dendrogram is disabled after enabling');

    await awaitCheck(() => model.findViewer(VIEWER_TYPE.DENDROGRAM) !== null,
      'Dendrogram is not present in the view after 5s delay', 5000);

    // Disable dendrogram
    model.settings = {showDendrogram: false};
    expect(model.settings.showDendrogram, false, 'Dendrogram is enabled after disabling');
    expect(model.findViewer(VIEWER_TYPE.DENDROGRAM) === null, true,
      'Dendrogram is present in the view after disabling');
  }, {skipReason: 'Need to find a way to replace _package variable to call for Bio function with tests'});
}, {clear: false});
