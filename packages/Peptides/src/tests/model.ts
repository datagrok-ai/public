import * as DG
  from 'datagrok-api/dg';

import {
  after,
  awaitCheck,
  before,
  category,
  delay,
  expect,
  expectFloat,
  test,
} from '@datagrok-libraries/test/src/test';
import {
  _package,
} from '../package-test';
import {
  PeptidesModel,
  VIEWER_TYPE,
} from '../model';
import {
  startAnalysis,
} from '../widgets/peptides';
import {
  scaleActivity,
} from '../utils/misc';
import {
  NOTATION,
} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {
  COLUMNS_NAMES,
  LST_COLUMN_NAMES,
  SCALING_METHODS,
} from '../utils/constants';
import {
  LogoSummaryTable,
} from '../viewers/logo-summary';
import {
  analyzeDemoPeptides,
  loadDemoPeptides,
  TEST_COLUMN_NAMES,
} from './utils';
import {
  getAggregatedColName,
} from '../utils/statistics';
import {
  MonomerPosition,
} from '../viewers/sar-viewer';
import {PeptideUtils} from '../peptideUtils';

category('Model: Settings', () => {
  let df: DG.DataFrame;
  let model: PeptidesModel;
  let activityCol: DG.Column<number>;
  let sequenceCol: DG.Column<string>;
  let clusterCol: DG.Column<any>;
  let scaledActivityCol: DG.Column<number>;

  const mutationCliffsDefaultParams = {
    maxMutations: 1,
    minActivityDelta: 0,
  };
  const mutationCliffsTestParams = {
    maxMutations: 2,
    minActivityDelta: 0.5,
  };

  before(async () => {
    await PeptideUtils.loadComponents();
    df = DG.DataFrame.fromCsv(await _package.files.readAsText('tests/HELM_small.csv'));
    df.name = 'HELM_small';
    //df.id ??= `HELM_small-${Date.now()}`;
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

    await delay(500);
  });

  after(async () => await delay(3000));

  test('Activity scaling', async () => {
    const getError = (row: number, method: SCALING_METHODS): string =>
      `Activity mismatch at row ${row} for scaling method '${method}'`;
    const tolerance = 0.0001;
    const origActivityData = model.df.getCol(model.settings?.activityColumnName!).getRawData();
    const scaledActivity = model.df.getCol(COLUMNS_NAMES.ACTIVITY);
    const dfLen = model.df.rowCount;

    // Check initial 'none' scaling
    let scaledActivityData = scaledActivity.getRawData();
    for (let i = 0; i < dfLen; i++)
      expectFloat(scaledActivityData[i], origActivityData[i], tolerance, getError(i, SCALING_METHODS.NONE));


    // Check 'lg' scaling
    scaledActivityData = scaleActivity(activityCol, SCALING_METHODS.LG).getRawData();
    for (let i = 0; i < dfLen; i++)
      expectFloat(scaledActivityData[i], Math.log10(origActivityData[i]), tolerance, getError(i, SCALING_METHODS.LG));


    // Check '-lg' scaling
    scaledActivityData = scaleActivity(activityCol, SCALING_METHODS.MINUS_LG).getRawData();
    for (let i = 0; i < dfLen; i++) {
      expectFloat(scaledActivityData[i], -Math.log10(origActivityData[i]), tolerance,
        getError(i, SCALING_METHODS.MINUS_LG));
    }

    // Check 'none' scaling
    scaledActivityData = scaledActivityCol.getRawData();
    for (let i = 0; i < dfLen; i++)
      expectFloat(scaledActivityData[i], origActivityData[i], tolerance, getError(i, SCALING_METHODS.NONE));
  });

  test('Mutation Cliffs', async () => {
    // Check default mutation cliffs parameters
    const sarViewer = model.findViewer(VIEWER_TYPE.SEQUENCE_VARIABILITY_MAP) as MonomerPosition;
    expect(sarViewer.maxMutations, mutationCliffsDefaultParams.maxMutations, `Max mutations mismatch: expected ` +
      `${mutationCliffsDefaultParams.maxMutations}, actual ${sarViewer.maxMutations}`);
    expect(sarViewer.minActivityDelta, mutationCliffsDefaultParams.minActivityDelta, `Min activity delta ` +
      `mismatch: expected ${mutationCliffsDefaultParams.minActivityDelta}, actual ${sarViewer.minActivityDelta}`);

    // Check test mutation cliffs parameters
    sarViewer.maxMutations = mutationCliffsTestParams.maxMutations;
    sarViewer.minActivityDelta = mutationCliffsTestParams.minActivityDelta;
    expect(sarViewer.maxMutations, mutationCliffsTestParams.maxMutations, `Max mutations mismatch: expected ` +
      `${mutationCliffsTestParams.maxMutations}, actual ${sarViewer.maxMutations}`);
    expect(sarViewer.minActivityDelta, mutationCliffsTestParams.minActivityDelta, `Min activity delta ` +
      `mismatch: expected ${mutationCliffsTestParams.minActivityDelta}, actual ${sarViewer.minActivityDelta}`);
  });

  test('Include columns', async () => {
    const columnName = 'rank';
    const testColumns = {[columnName]: DG.AGG.AVG};
    model.settings = {columns: testColumns};

    // Include column
    const lstViewer = model.findViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE) as LogoSummaryTable;
    // lstViewer.setOptions({columns: testColumns});
    const aggColName = getAggregatedColName(DG.AGG.AVG, columnName);
    expect(lstViewer.viewerGrid.col(aggColName) !== null, true, `Expected to include column '${columnName}' in ` +
      `${VIEWER_TYPE.LOGO_SUMMARY_TABLE} but it is absent`);

    // Remove column
    model.settings = {columns: {}};
    expect(Object.keys(model.settings?.columns!).length, 0,
      `Expected to remove all column aggregations but columns {${Object.keys(model.settings?.columns!).join(' & ')}} ` +
      `are still included`);

    expect(lstViewer.viewerGrid.col(aggColName) === null, true, `Expected to remove column '${columnName}' from ` +
      `${VIEWER_TYPE.LOGO_SUMMARY_TABLE} but it is still present`);
  });

  test('Dendrogram', async () => {
    // Enable dendrogram
    model.settings = {showDendrogram: true};
    expect(model.settings?.showDendrogram, true, 'Dendrogram is disabled after enabling');

    await awaitCheck(() => model.dendrogramCloseButton !== null,
      'Dendrogram is not present in the view after 5s delay', 5000);

    // Disable dendrogram
    model.settings = {showDendrogram: false};
    expect(model.settings?.showDendrogram, false, 'Dendrogram is enabled after disabling');
    expect(model.dendrogramCloseButton === null, true,
      'Dendrogram is present in the view after disabling');
  }, {skipReason: 'Need to find a way to replace _package variable to call for Bio function with tests'});
}, {clear: false});

category('Model: MCL clusters', () => {
  before(async () => {
    await PeptideUtils.loadComponents();
  });

  /** Checks every Logo Summary Table row against the source rows of its cluster and the rest. */
  function expectClusterSummary(lst: LogoSummaryTable): void {
    const df = lst.dataFrame;
    const clusters = df.getCol(lst.clustersColumnName);
    const activity = df.getCol(COLUMNS_NAMES.ACTIVITY);
    const summary = lst.logoSummaryTable;
    const pValues = summary.getCol(LST_COLUMN_NAMES.P_VALUE);
    const mean = (values: number[]): number => values.reduce((sum, value) => sum + value, 0) / values.length;
    let members = 0;
    for (let row = 0; row < summary.rowCount; row++) {
      const cluster = String(summary.get(LST_COLUMN_NAMES.CLUSTER, row));
      const selected: number[] = [];
      const rest: number[] = [];
      for (let i = 0; i < df.rowCount; i++)
        (String(clusters.get(i)) === cluster ? selected : rest).push(activity.get(i));
      members += selected.length;
      expect(summary.get(LST_COLUMN_NAMES.MEMBERS, row), selected.length, `members of cluster ${cluster}`);
      const meanDifference = rest.length === 0 ? 0 : mean(selected) - mean(rest);
      expectFloat(summary.get(LST_COLUMN_NAMES.MEAN_DIFFERENCE, row), meanDifference,
        Math.max(Math.abs(meanDifference) * 1e-4, 1e-9), `mean difference of cluster ${cluster}`);
      const pValue = pValues.get(row);
      expect(selected.length > 1 && rest.length > 1 ? pValue >= 0 && pValue <= 1 : pValues.isNone(row), true,
        `p-value ${pValue} of cluster ${cluster}, which has no t-test with fewer than two rows on a side`);
    }
    expect(members, df.rowCount, 'members total');
  }

  // 96 splits the first 200 peptides into five clusters on some stands and six on others
  for (const [threshold, fewest, most] of [[10, 1, 1], [75, 2, 2], [96, 5, 6]]) {
    test(`Similarity threshold ${threshold}`, async () => {
      const model = await analyzeDemoPeptides(await loadDemoPeptides(200), threshold);
      expect(model.findViewer(VIEWER_TYPE.MCL)?.getWidgetStatus().values?.['completed threshold'], threshold,
        'completed MCL threshold');
      const lst = model.findViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE) as LogoSummaryTable;
      const clusterCount = lst.logoSummaryTable.rowCount;
      expect(clusterCount >= fewest && clusterCount <= most, true,
        `${clusterCount} clusters at threshold ${threshold}, expected ${fewest}-${most}`);
      expectClusterSummary(lst);
    }, {timeout: 120000});
  }
});
