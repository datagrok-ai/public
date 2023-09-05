import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {awaitCheck, before, category, expect, test, testViewer} from '@datagrok-libraries/utils/src/test';
import {aligned1} from './test-data';
import {CLUSTER_TYPE, PeptidesModel, VIEWER_TYPE} from '../model';
import {_package} from '../package-test';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {scaleActivity} from '../utils/misc';
import {startAnalysis} from '../widgets/peptides';
import {MONOMER_POSITION_MODE, MonomerPosition, MostPotentResidues, showTooltip} from '../viewers/sar-viewer';
import {SCALING_METHODS} from '../utils/constants';
import {LST_PROPERTIES, LogoSummaryTable} from '../viewers/logo-summary';
import {PositionHeight} from '@datagrok-libraries/bio/src/viewers/web-logo';
import {TEST_COLUMN_NAMES} from './utils';

category('Viewers: Basic', () => {
  const df = DG.DataFrame.fromCsv(aligned1);
  const viewers = DG.Func.find({package: 'Peptides', tags: ['viewer']}).map((f) => f.friendlyName);
  for (const v of viewers) {
    test(v, async () => {
      await testViewer(v, df.clone(), {detectSemanticTypes: true});
    }, {skipReason: 'GROK-11534'});
  }
});

category('Viewers: Monomer-Position', () => {
  let df: DG.DataFrame;
  let model: PeptidesModel;
  let activityCol: DG.Column<number>;
  let sequenceCol: DG.Column<string>;
  let clusterCol: DG.Column<any>;
  let scaledActivityCol: DG.Column<number>;
  let mpViewer: MonomerPosition;

  before(async () => {
    df = DG.DataFrame.fromCsv(await _package.files.readAsText('tests/HELM_small.csv'));
    activityCol = df.getCol(TEST_COLUMN_NAMES.ACTIVITY);
    sequenceCol = df.getCol(TEST_COLUMN_NAMES.SEQUENCE);
    sequenceCol.semType = DG.SEMTYPE.MACROMOLECULE;
    sequenceCol.setTag(DG.TAGS.UNITS, NOTATION.HELM);
    scaledActivityCol = scaleActivity(activityCol, SCALING_METHODS.NONE);
    clusterCol = df.getCol(TEST_COLUMN_NAMES.CLUSTER);
    const tempModel = await startAnalysis(
      activityCol, sequenceCol, clusterCol, df, scaledActivityCol, SCALING_METHODS.NONE);
    if (tempModel === null)
      throw new Error('Model is null');
    model = tempModel;
    let overlayInit = false;
    model._analysisView!.grid.onAfterDrawOverlay.subscribe(() => overlayInit = true);

    mpViewer = model.findViewer(VIEWER_TYPE.MONOMER_POSITION) as MonomerPosition;

    // Ensure grid finished initializing to prevent Unhandled exceptions
    let accrodionInit = false;
    grok.events.onAccordionConstructed.subscribe((_) => accrodionInit = true);
    await awaitCheck(() => model!.df.currentRowIdx === 0, 'Grid cell never finished initializing', 2000);
    await awaitCheck(() => grok.shell.o instanceof DG.Column, 'Shell object never changed', 2000);
    await awaitCheck(() => accrodionInit, 'Accordion never finished initializing', 2000);
    await awaitCheck(() => overlayInit, 'Overlay never finished initializing', 2000);
  });

  test('Tooltip', async () => {
    const cellCoordinates = {col: '9', row: 6};
    const gc = mpViewer.viewerGrid.cell(cellCoordinates.col, cellCoordinates.row);
    const mp = mpViewer.getMonomerPosition(gc);
    expect(showTooltip(mp, 0, 0, model), true,
      `Tooltip is not shown for grid cell at column '${cellCoordinates.col}', row ${cellCoordinates.row}`);
  });

  test('Modes', async () => {
    if (mpViewer === null)
      throw new Error('Monomer-Position viewer doesn\'t exist');

    expect(mpViewer.mode, MONOMER_POSITION_MODE.MUTATION_CLIFFS,
      `Default Monomer-Position mode is not ${MONOMER_POSITION_MODE.MUTATION_CLIFFS}`);

    mpViewer.mode = MONOMER_POSITION_MODE.INVARIANT_MAP;
    expect(mpViewer.mode, MONOMER_POSITION_MODE.INVARIANT_MAP,
      `Monomer-Position mode is not ${MONOMER_POSITION_MODE.INVARIANT_MAP} after switching`);

    mpViewer.mode = MONOMER_POSITION_MODE.MUTATION_CLIFFS;
    expect(mpViewer.mode, MONOMER_POSITION_MODE.MUTATION_CLIFFS,
      `Monomer-Position mode is not ${MONOMER_POSITION_MODE.MUTATION_CLIFFS} after switching`);
  });
}, {clear: false});

category('Viewers: Most Potent Residues', () => {
  let df: DG.DataFrame;
  let model: PeptidesModel;
  let activityCol: DG.Column<number>;
  let sequenceCol: DG.Column<string>;
  let clusterCol: DG.Column<any>;
  let scaledActivityCol: DG.Column<number>;
  let mprViewer: MostPotentResidues;

  before(async () => {
    df = DG.DataFrame.fromCsv(await _package.files.readAsText('tests/HELM_small.csv'));
    activityCol = df.getCol(TEST_COLUMN_NAMES.ACTIVITY);
    sequenceCol = df.getCol(TEST_COLUMN_NAMES.SEQUENCE);
    sequenceCol.semType = DG.SEMTYPE.MACROMOLECULE;
    sequenceCol.setTag(DG.TAGS.UNITS, NOTATION.HELM);
    scaledActivityCol = scaleActivity(activityCol, SCALING_METHODS.NONE);
    clusterCol = df.getCol(TEST_COLUMN_NAMES.CLUSTER);
    const tempModel = await startAnalysis(
      activityCol, sequenceCol, clusterCol, df, scaledActivityCol, SCALING_METHODS.NONE);
    if (tempModel === null)
      throw new Error('Model is null');
    model = tempModel;
    let overlayInit = false;
    model._analysisView!.grid.onAfterDrawOverlay.subscribe(() => overlayInit = true);

    mprViewer = model.findViewer(VIEWER_TYPE.MOST_POTENT_RESIDUES) as MostPotentResidues;

    // Ensure grid finished initializing to prevent Unhandled exceptions
    let accrodionInit = false;
    grok.events.onAccordionConstructed.subscribe((_) => accrodionInit = true);
    await awaitCheck(() => model!.df.currentRowIdx === 0, 'Grid cell never finished initializing', 2000);
    await awaitCheck(() => grok.shell.o instanceof DG.Column, 'Shell object never changed', 2000);
    await awaitCheck(() => accrodionInit, 'Accordion never finished initializing', 2000);
    await awaitCheck(() => overlayInit, 'Overlay never finished initializing', 2000);
  });

  test('Tooltip', async () => {
    const cellCoordinates = {col: 'Diff', row: 6};
    const gc = mprViewer.viewerGrid.cell(cellCoordinates.col, cellCoordinates.row);
    const mp = mprViewer.getMonomerPosition(gc);
    expect(showTooltip(mp, 0, 0, model), true,
      `Tooltip is not shown for grid cell at column '${cellCoordinates.col}', row ${cellCoordinates.row}`);
  });
});

category('Viewers: Logo Summary Table', () => {
  let df: DG.DataFrame;
  let model: PeptidesModel;
  let activityCol: DG.Column<number>;
  let sequenceCol: DG.Column<string>;
  let clusterCol: DG.Column<any>;
  let scaledActivityCol: DG.Column<number>;
  let lstViewer: LogoSummaryTable;

  before(async () => {
    df = DG.DataFrame.fromCsv(await _package.files.readAsText('tests/HELM_small.csv'));
    activityCol = df.getCol(TEST_COLUMN_NAMES.ACTIVITY);
    sequenceCol = df.getCol(TEST_COLUMN_NAMES.SEQUENCE);
    sequenceCol.semType = DG.SEMTYPE.MACROMOLECULE;
    sequenceCol.setTag(DG.TAGS.UNITS, NOTATION.HELM);
    scaledActivityCol = scaleActivity(activityCol, SCALING_METHODS.NONE);
    clusterCol = df.getCol(TEST_COLUMN_NAMES.CLUSTER);
    const tempModel = await startAnalysis(
      activityCol, sequenceCol, clusterCol, df, scaledActivityCol, SCALING_METHODS.NONE);
    if (tempModel === null)
      throw new Error('Model is null');
    model = tempModel;
    let overlayInit = false;
    model._analysisView!.grid.onAfterDrawOverlay.subscribe(() => overlayInit = true);

    lstViewer = model.findViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE) as LogoSummaryTable;

    // Ensure grid finished initializing to prevent Unhandled exceptions
    let accrodionInit = false;
    grok.events.onAccordionConstructed.subscribe((_) => accrodionInit = true);
    await awaitCheck(() => model!.df.currentRowIdx === 0, 'Grid cell never finished initializing', 2000);
    await awaitCheck(() => grok.shell.o instanceof DG.Column, 'Shell object never changed', 2000);
    await awaitCheck(() => accrodionInit, 'Accordion never finished initializing', 2000);
    await awaitCheck(() => overlayInit, 'Overlay never finished initializing', 2000);
  });

  test('Properties', async () => {
    // Change Logo Summary Table Web Logo Mode property to full
    const webLogoMode = lstViewer.getProperty(LST_PROPERTIES.WEB_LOGO_MODE);
    webLogoMode!.set(lstViewer, PositionHeight.full);
    expect(lstViewer.webLogoMode, PositionHeight.full,
      `Web Logo Mode property is not changed to ${PositionHeight.full}, got ${lstViewer.webLogoMode} instead`);

    // Change Logo Summary Table Members Ratio Threshold proprty to 0
    const threshold = 0;
    const membersRatioThreshold = lstViewer.getProperty(LST_PROPERTIES.MEMBERS_RATIO_THRESHOLD);
    membersRatioThreshold!.set(lstViewer, threshold);
    expect(lstViewer.membersRatioThreshold, threshold,
      `Members Ratio Threshold property is not changed to 0, got ${lstViewer.membersRatioThreshold} instead`);
    expect(lstViewer.viewerGrid.table.filter.anyTrue, true, `Expected to filter out all rows, but ` +
      `${lstViewer.viewerGrid.table.filter.trueCount} rows are left unfiltered`);
  });

  test('Tooltip', async () => {
    const cluster = '0';
    const tooltipElement = lstViewer.showTooltip({name: cluster, type: CLUSTER_TYPE.ORIGINAL}, 0, 0);
    expect(tooltipElement !== null, true, `Tooltip is not shown for cluster '${cluster}'`);
  });
}, {clear: false});
