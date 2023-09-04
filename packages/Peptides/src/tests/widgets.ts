import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, test, before, expect, awaitCheck} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import {PeptidesModel, VIEWER_TYPE} from '../model';
import {scaleActivity} from '../utils/misc';
import {startAnalysis} from '../widgets/peptides';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import * as C from '../utils/constants';
import {PANES_INPUTS, SETTINGS_PANES, getSettingsDialog} from '../widgets/settings';
import {getDistributionWidget} from '../widgets/distribution';
import {mutationCliffsWidget} from '../widgets/mutation-cliffs';
import {TEST_COLUMN_NAMES} from './utils';
import wu from 'wu';
import {LogoSummaryTable} from '../viewers/logo-summary';

category('Widgets: Settings', () => {
  let df: DG.DataFrame;
  let model: PeptidesModel;
  let activityCol: DG.Column<number>;
  let sequenceCol: DG.Column<string>;
  let clusterCol: DG.Column<any>;
  let scaledActivityCol: DG.Column<number>;

  before(async () => {
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

  test('UI', async () => {
    const settingsElements = getSettingsDialog(model);

    // Check number of panes
    const panes = settingsElements.accordion.panes.map((pane) => pane.name);
    expect(panes.length, 4, `Expected 4 panes, got ${settingsElements.accordion.panes.length}`);
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

  test('UI', async () => {
    getDistributionWidget(model.df, model);
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

  test('UI', async () => {
    mutationCliffsWidget(model.df, model);
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

  test('New view', async () => {
    // Set compound bitset: filter out 2 rows and select 1 among them
    const filter = model.df.filter;
    filter.setAll(false, false);
    filter.set(0, true, false);
    filter.set(1, true, false);

    const selection = model.df.selection;
    selection.set(0, true, false);

    const newViewId = model.createNewView();
    const currentTable = grok.shell.t;

    expect(currentTable.getTag(C.TAGS.MULTIPLE_VIEWS), '1', 'Current table is expected to have multiple views tag');
    expect(currentTable.getTag(C.TAGS.UUID), newViewId, 'Current table is expected to have the same UUID as new view');
    expect(currentTable.rowCount, 1, 'Current table is expected to have 1 row');

    await awaitCheck(() => currentTable.currentRowIdx === 0, 'Grid never finished initializing', 2000);

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
    expect(wu(model.customClusters).toArray().length, 0, 'Expected to have 0 custom clusters before creating one');

    // Create custom cluster
    lstViewer.clusterFromSelection();
    const customClusterList = wu(model.customClusters).toArray();
    expect(customClusterList.length, 1, 'Expected to have 1 custom cluster');
    const clustName = customClusterList[0].name;
    expect(model.df.col(clustName) !== null, true,
      'Expected to have custom cluster column in the table');
    expect(lstViewer.viewerGrid.table.getCol(C.LST_COLUMN_NAMES.CLUSTER).categories.indexOf(clustName) !== -1, true,
      'Expected to have custom cluster in the Logo Summary Table');

    // Remove custom cluster
    model.modifyClusterSelection(clustName);
    lstViewer.removeCluster();
    expect(wu(model.customClusters).toArray().length, 0, 'Expected to have 0 custom clusters after removing one');
    expect(model.df.col(clustName) === null, true,
      'Expected to have no custom cluster column in the table');
    expect(lstViewer.viewerGrid.table.getCol(C.LST_COLUMN_NAMES.CLUSTER).categories.indexOf(clustName) === -1, true,
      'Expected to have no custom cluster in the Logo Summary Table');
  });
}, {clear: false});
