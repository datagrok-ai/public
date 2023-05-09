import * as DG from 'datagrok-api/dg';

import {before, category, expect, test, testViewer} from '@datagrok-libraries/utils/src/test';
import {aligned1} from './test-data';
import {PeptidesModel, VIEWER_TYPE} from '../model';
import {_package} from '../package-test';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {scaleActivity} from '../utils/misc';
import {startAnalysis} from '../widgets/peptides';
import {MONOMER_POSITION_MODE, MonomerPosition, MostPotentResiduesViewer, showTooltip} from '../viewers/sar-viewer';
import {SCALING_METHODS} from '../utils/constants';
import {LST_PROPERTIES, LogoSummaryTable} from '../viewers/logo-summary';
import {PositionHeight} from '@datagrok-libraries/bio/src/viewers/web-logo';


category('Viewers: Basic', () => {
  const df = DG.DataFrame.fromCsv(aligned1);
  const viewers = DG.Func.find({package: 'Peptides', tags: ['viewer']}).map((f) => f.friendlyName);
  for (const v of viewers) {
    test(v, async () => {
      await testViewer(v, df.clone(), true);
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
    activityCol = df.getCol('activity');
    sequenceCol = df.getCol('sequence');
    sequenceCol.semType = DG.SEMTYPE.MACROMOLECULE;
    sequenceCol.setTag(DG.TAGS.UNITS, NOTATION.HELM);
    scaledActivityCol = scaleActivity(activityCol, SCALING_METHODS.NONE);
    clusterCol = df.getCol('cluster');
    const tempModel = await startAnalysis(
      activityCol, sequenceCol, clusterCol, df, scaledActivityCol, SCALING_METHODS.NONE);
    if (tempModel === null)
      throw new Error('Model is null');
    model = tempModel;
    mpViewer = model.findViewer(VIEWER_TYPE.MONOMER_POSITION) as MonomerPosition;
  });

  test('Tooltip', async () => {
    const cellCoordinates = {col: '9', row: 6};
    const gc = mpViewer.viewerGrid.cell(cellCoordinates.col, cellCoordinates.row);
    expect(showTooltip(gc, 0, 0, model), true,
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
});

category('Viewers: Most Potent Residues', () => {
  let df: DG.DataFrame;
  let model: PeptidesModel;
  let activityCol: DG.Column<number>;
  let sequenceCol: DG.Column<string>;
  let clusterCol: DG.Column<any>;
  let scaledActivityCol: DG.Column<number>;
  let mprViewer: MostPotentResiduesViewer;

  before(async () => {
    df = DG.DataFrame.fromCsv(await _package.files.readAsText('tests/HELM_small.csv'));
    activityCol = df.getCol('activity');
    sequenceCol = df.getCol('sequence');
    sequenceCol.semType = DG.SEMTYPE.MACROMOLECULE;
    sequenceCol.setTag(DG.TAGS.UNITS, NOTATION.HELM);
    scaledActivityCol = scaleActivity(activityCol, SCALING_METHODS.NONE);
    clusterCol = df.getCol('cluster');
    const tempModel = await startAnalysis(
      activityCol, sequenceCol, clusterCol, df, scaledActivityCol, SCALING_METHODS.NONE);
    if (tempModel === null)
      throw new Error('Model is null');
    model = tempModel;
    mprViewer = model.findViewer(VIEWER_TYPE.MOST_POTENT_RESIDUES) as MostPotentResiduesViewer;
  });

  test('Tooltip', async () => {
    const cellCoordinates = {col: 'Diff', row: 6};
    const gc = mprViewer.viewerGrid.cell(cellCoordinates.col, cellCoordinates.row);
    expect(showTooltip(gc, 0, 0, model), true,
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
    activityCol = df.getCol('activity');
    sequenceCol = df.getCol('sequence');
    sequenceCol.semType = DG.SEMTYPE.MACROMOLECULE;
    sequenceCol.setTag(DG.TAGS.UNITS, NOTATION.HELM);
    scaledActivityCol = scaleActivity(activityCol, SCALING_METHODS.NONE);
    clusterCol = df.getCol('cluster');
    const tempModel = await startAnalysis(
      activityCol, sequenceCol, clusterCol, df, scaledActivityCol, SCALING_METHODS.NONE);
    if (tempModel === null)
      throw new Error('Model is null');
    model = tempModel;
    lstViewer = model.findViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE) as LogoSummaryTable;
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
    const tooltipElement = lstViewer.showTooltip(cluster, 0, 0);
    expect(tooltipElement !== null, true, `Tooltip is not shown for cluster '${cluster}'`);
  });
});
