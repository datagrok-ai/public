import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import {PeptidesModel, VIEWER_TYPE} from '../model';
import {scaleActivity} from '../utils/misc';
import {startAnalysis} from '../widgets/peptides';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import * as C from '../utils/constants';
import {getSettingsDialog, PANES_INPUTS, SETTINGS_PANES} from '../widgets/settings';
import {getDistributionWidget} from '../widgets/distribution';
import {mutationCliffsWidget} from '../widgets/mutation-cliffs';
import {TEST_COLUMN_NAMES} from './utils';
import wu from 'wu';
import {CLUSTER_TYPE, LogoSummaryTable} from '../viewers/logo-summary';
import {MonomerPosition} from '../viewers/sar-viewer';
import {PeptideUtils} from '../peptideUtils';

category('Widgets: Settings', () => {
  let df: DG.DataFrame;
  let model: PeptidesModel;
  let activityCol: DG.Column<number>;
  let sequenceCol: DG.Column<string>;
  let clusterCol: DG.Column<any>;
  let scaledActivityCol: DG.Column<number>;

  before(async () => {
    await PeptideUtils.loadSeqHelper();

    df = DG.DataFrame.fromCsv(await _package.files.readAsText('tests/HELM_small.csv'));
    activityCol = df.getCol(TEST_COLUMN_NAMES.ACTIVITY);
    sequenceCol = df.getCol(TEST_COLUMN_NAMES.SEQUENCE);
    sequenceCol.semType = DG.SEMTYPE.MACROMOLECULE;
    sequenceCol.setTag(DG.TAGS.UNITS, NOTATION.HELM);
    scaledActivityCol = scaleActivity(activityCol, C.SCALING_METHODS.NONE);
    clusterCol = df.getCol(TEST_COLUMN_NAMES.CLUSTER);
    const tempModel = await startAnalysis(activityCol, sequenceCol, clusterCol, df, scaledActivityCol,
      C.SCALING_METHODS.NONE);
    if (tempModel === null)
      throw new Error('Model is null');
    model = tempModel;

    await delay(500);
  });

  after(async () => await delay(3000));

  test('UI', async () => {
    const settingsElements = getSettingsDialog(model);

    // Check number of panes
    const panes = settingsElements.accordion.panes.map((pane) => pane.name);
    expect(panes.length, 5, `Expected 5 panes, got ${settingsElements.accordion.panes.length}`);
    for (const paneName of Object.values(SETTINGS_PANES))
      expect(panes.includes(paneName), true, `Pane ${paneName} is missing`);

    // Check inputs in each pane
    for (const paneName of Object.values(SETTINGS_PANES)) {
      const paneInputs = settingsElements.inputs[paneName].map((input) => input.caption);
      for (const inputName of Object.values(PANES_INPUTS[paneName]))
        expect(paneInputs.includes(inputName), true, `Input ${inputName} is missing from ${paneName}`);
    }
  });
});

category('Widgets: Distribution panel', () => {
  let df: DG.DataFrame;
  let model: PeptidesModel;
  let activityCol: DG.Column<number>;
  let sequenceCol: DG.Column<string>;
  let clusterCol: DG.Column<any>;
  let scaledActivityCol: DG.Column<number>;

  before(async () => {
    await PeptideUtils.loadSeqHelper();

    df = DG.DataFrame.fromCsv(await _package.files.readAsText('tests/HELM_small.csv'));
    activityCol = df.getCol(TEST_COLUMN_NAMES.ACTIVITY);
    sequenceCol = df.getCol(TEST_COLUMN_NAMES.SEQUENCE);
    sequenceCol.semType = DG.SEMTYPE.MACROMOLECULE;
    sequenceCol.setTag(DG.TAGS.UNITS, NOTATION.HELM);
    scaledActivityCol = scaleActivity(activityCol, C.SCALING_METHODS.NONE);
    clusterCol = df.getCol(TEST_COLUMN_NAMES.CLUSTER);
    const tempModel = await startAnalysis(activityCol, sequenceCol, clusterCol, df, scaledActivityCol,
      C.SCALING_METHODS.NONE);
    if (tempModel === null)
      throw new Error('Model is null');
    model = tempModel;

    await delay(500);
  });

  after(async () => await delay(3000));

  test('UI', async () => {
    model.df.selection.set(0, true);
    await delay(1000);
    const lstViewer = model.findViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE) as LogoSummaryTable | null;
    getDistributionWidget(model.df, {
      peptideSelection: DG.BitSet.create(model.df.rowCount), columns: model.settings!.columns!,
      activityCol: scaledActivityCol, clusterSelection: lstViewer!.clusterSelection, clusterColName: clusterCol.name,
      monomerPositionSelection: model.webLogoSelection,
    });
  });
});

category('Widgets: Mutation cliffs', () => {
  let df: DG.DataFrame;
  let model: PeptidesModel;
  let activityCol: DG.Column<number>;
  let sequenceCol: DG.Column<string>;
  let clusterCol: DG.Column<any>;
  let scaledActivityCol: DG.Column<number>;

  before(async () => {
    await PeptideUtils.loadSeqHelper();

    df = DG.DataFrame.fromCsv(await _package.files.readAsText('tests/HELM_small.csv'));
    activityCol = df.getCol(TEST_COLUMN_NAMES.ACTIVITY);
    sequenceCol = df.getCol(TEST_COLUMN_NAMES.SEQUENCE);
    sequenceCol.semType = DG.SEMTYPE.MACROMOLECULE;
    sequenceCol.setTag(DG.TAGS.UNITS, NOTATION.HELM);
    scaledActivityCol = scaleActivity(activityCol, C.SCALING_METHODS.NONE);
    clusterCol = df.getCol(TEST_COLUMN_NAMES.CLUSTER);
    const tempModel = await startAnalysis(activityCol, sequenceCol, clusterCol, df, scaledActivityCol,
      C.SCALING_METHODS.NONE);
    if (tempModel === null)
      throw new Error('Model is null');
    model = tempModel;

    await delay(500);
  });

  after(async () => await delay(3000));

  test('UI', async () => {
    const sarViewer = model.findViewer(VIEWER_TYPE.SEQUENCE_VARIABILITY_MAP) as MonomerPosition;
    sarViewer.keyPressed = true; //required to emulate cell selection
    sarViewer._viewerGrid!.dataFrame.currentCell = sarViewer._viewerGrid?.dataFrame.cell(0, '1')!;
    await delay(1000);
    mutationCliffsWidget(model.df, {
      mutationCliffs: sarViewer.mutationCliffs!,
      mutationCliffsSelection: sarViewer.mutationCliffsSelection,
      gridColumns: model.analysisView.grid.columns,
      sequenceColumnName: sarViewer.sequenceColumnName,
      positionColumns: sarViewer.positionColumns,
      activityCol: scaledActivityCol,
    });
  });
});

category('Widgets: Actions', () => {
  let df: DG.DataFrame;
  let model: PeptidesModel;
  let activityCol: DG.Column<number>;
  let sequenceCol: DG.Column<string>;
  let clusterCol: DG.Column<any>;
  let scaledActivityCol: DG.Column<number>;

  before(async () => {
    await PeptideUtils.loadSeqHelper();

    df = DG.DataFrame.fromCsv(await _package.files.readAsText('tests/HELM_small.csv'));
    await df.meta.detectSemanticTypes();
    activityCol = df.getCol(TEST_COLUMN_NAMES.ACTIVITY);
    sequenceCol = df.getCol(TEST_COLUMN_NAMES.SEQUENCE);
    sequenceCol.semType = DG.SEMTYPE.MACROMOLECULE;
    sequenceCol.setTag(DG.TAGS.UNITS, NOTATION.HELM);
    scaledActivityCol = scaleActivity(activityCol, C.SCALING_METHODS.NONE);
    clusterCol = df.getCol(TEST_COLUMN_NAMES.CLUSTER);
    const tempModel = await startAnalysis(activityCol, sequenceCol, clusterCol, df, scaledActivityCol,
      C.SCALING_METHODS.NONE);
    if (tempModel === null)
      throw new Error('Model is null');
    model = tempModel;

    await delay(500);
  });

  after(async () => await delay(3000));

  test('New view', async () => {
    // Set compound bitset: filter out 2 rows and select 1 among them
    const filter = model.df.filter;
    filter.setAll(false, false);
    filter.set(0, true, false);
    filter.set(1, true, false);

    const selection = model.df.selection;
    selection.set(0, true, false);
    selection.set(1, true, false);

    const newViewId = model.createNewView();
    const currentTable = grok.shell.t;

    expect(currentTable.getTag(C.TAGS.MULTIPLE_VIEWS), '1', 'Current table is expected to have multiple views tag');
    expect(currentTable.getTag(DG.TAGS.ID), newViewId, 'Current table is expected to have the same UUID as new view');
    expect(currentTable.rowCount, 2, 'Current table is expected to have 2 rows');

    await delay(500);
    const currentTableModel = currentTable.temp[PeptidesModel.modelName] as PeptidesModel;
    const lstViewer = currentTableModel.findViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE);
    expect(lstViewer !== null, true, 'New view is expected to have Logo Summary Table viewer attached');
  });

  test('Custom clusters', async () => {
    // Set compound bitset: filter out 2 rows and select 1 among them
    const filter = model.df.filter;
    filter.setAll(false, false);
    filter.set(0, true, false);
    filter.set(1, true, false);

    const selection = model.df.selection;
    selection.set(0, true, false);

    const lstViewer = model.findViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE) as LogoSummaryTable | null;
    if (lstViewer === null)
      throw new Error('Logo summary table viewer is not found');

    // Check that custom clusters are not created yet
    expect(wu(lstViewer.customClusters).toArray().length, 0, 'Expected to have 0 custom clusters before creating one');

    // Create custom cluster
    lstViewer.clusterFromSelection();
    const customClusterList = wu(lstViewer.customClusters).toArray();
    expect(customClusterList.length, 1, 'Expected to have 1 custom cluster');
    const clustName = customClusterList[0].name;
    expect(model.df.col(clustName) !== null, true, 'Expected to have custom cluster column in the table');
    expect(lstViewer.viewerGrid.table.getCol(C.LST_COLUMN_NAMES.CLUSTER).categories.indexOf(clustName) !== -1, true,
      'Expected to have custom cluster in the Logo Summary Table');

    // Remove custom cluster
    lstViewer.modifyClusterSelection({monomerOrCluster: clustName, positionOrClusterType: CLUSTER_TYPE.CUSTOM});
    lstViewer.removeCluster();
    expect(wu(lstViewer.customClusters).toArray().length, 0, 'Expected to have 0 custom clusters after removing one');
    expect(model.df.col(clustName) === null, true,
      'Expected to have no custom cluster column in the table');
    expect(lstViewer.viewerGrid.table.getCol(C.LST_COLUMN_NAMES.CLUSTER).categories.indexOf(clustName) === -1, true,
      'Expected to have no custom cluster in the Logo Summary Table');
  });
}, {clear: false});
