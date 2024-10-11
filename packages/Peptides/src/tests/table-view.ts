import * as DG from 'datagrok-api/dg';

import {after, awaitCheck, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import {PeptidesModel, VIEWER_TYPE} from '../model';
import {startAnalysis} from '../widgets/peptides';
import {initSelection, scaleActivity} from '../utils/misc';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {COLUMNS_NAMES, SCALING_METHODS} from '../utils/constants';
import {TEST_COLUMN_NAMES} from './utils';
import {CLUSTER_TYPE, LogoSummaryTable} from '../viewers/logo-summary';
import {MonomerPosition} from '../viewers/sar-viewer';

category('Table view', () => {
  let df: DG.DataFrame;
  let model: PeptidesModel;
  let activityCol: DG.Column<number>;
  let sequenceCol: DG.Column<string>;
  let clusterCol: DG.Column<any>;
  let scaledActivityCol: DG.Column<number>;
  const scaling = 'none' as SCALING_METHODS;

  const firstPair = {monomerOrCluster: 'Aze', positionOrClusterType: '10', mcCount: 3, imCount: 1};
  const secondPair = {monomerOrCluster: 'meI', positionOrClusterType: '1', mcCount: 2, imCount: 10};
  const firstCluster = {monomerOrCluster: '0', positionOrClusterType: CLUSTER_TYPE.ORIGINAL, count: 3};
  const secondCluster = {monomerOrCluster: '1', positionOrClusterType: CLUSTER_TYPE.ORIGINAL, count: 3};

  before(async () => {
    df = DG.DataFrame.fromCsv(await _package.files.readAsText('tests/HELM_small.csv'));
    activityCol = df.getCol(TEST_COLUMN_NAMES.ACTIVITY);
    sequenceCol = df.getCol(TEST_COLUMN_NAMES.SEQUENCE);
    sequenceCol.semType = DG.SEMTYPE.MACROMOLECULE;
    sequenceCol.setTag(DG.TAGS.UNITS, NOTATION.HELM);
    scaledActivityCol = scaleActivity(activityCol, scaling);
    clusterCol = df.getCol(TEST_COLUMN_NAMES.CLUSTER);
    const tempModel = await startAnalysis(activityCol, sequenceCol, clusterCol, df, scaledActivityCol, scaling);
    if (tempModel === null)
      throw new Error('Model is null');
    model = tempModel;

    await delay(500);
  });

  after(async () => await delay(3000));

  test('Visible columns', async () => {
    const gridCols = model.analysisView.grid.columns;
    const posCols = model.positionColumns!.map((col) => col.name);
    for (let colIdx = 1; colIdx < gridCols.length; colIdx++) {
      const col = gridCols.byIndex(colIdx)!;
      const tableColName = col.column!.name;
      const expectedVisibility = posCols.includes(tableColName) || (tableColName === COLUMNS_NAMES.ACTIVITY);
      expect(col.visible, expectedVisibility, `Column ${tableColName} is visible === ${col.visible} but should be ` +
        `${expectedVisibility}`);
    }
  });

  //TODO: split into separate tests for Mutation Cliffs and Logo Summary Table
  test('Mutation Cliffs selection', async () => {
    const selection = model.df.selection;

    const mpViewer = model.findViewer(VIEWER_TYPE.SEQUENCE_VARIABILITY_MAP) as MonomerPosition;
    await awaitCheck(() => mpViewer.mutationCliffs !== null, 'mutation cliffs haven\'t been generated', 2000);

    for (const [position, selectedMonomers] of Object.entries(mpViewer.mutationCliffsSelection))
      expect(selectedMonomers.length, 0, `Selection is not empty for position ${position} after initialization`);

    // Select first monomer-position pair
    mpViewer.modifyMutationCliffsSelection(firstPair);
    expect(mpViewer.mutationCliffsSelection[firstPair.positionOrClusterType].includes(firstPair.monomerOrCluster), true,
      `Monomer ${firstPair.monomerOrCluster} is not selected at position ${firstPair.positionOrClusterType}`);
    expect(selection.trueCount, firstPair.mcCount, `Selection count is not equal to ${firstPair.mcCount} ` +
      `for monomer ${firstPair.monomerOrCluster} at position ${firstPair.positionOrClusterType}`);

    // Select second monomer-position pair
    mpViewer.modifyMutationCliffsSelection(secondPair, {shiftPressed: true, ctrlPressed: false});
    expect(mpViewer.mutationCliffsSelection[secondPair.positionOrClusterType].includes(secondPair.monomerOrCluster),
      true, `Monomer ${secondPair.monomerOrCluster} is not selected at position ${secondPair.positionOrClusterType}`);
    expect(selection.trueCount, secondPair.mcCount + firstPair.mcCount, `Selection count is not equal ` +
      `to ${secondPair.mcCount + firstPair.mcCount} for monomer ${secondPair.monomerOrCluster} at ` +
      `position ${secondPair.positionOrClusterType}`);

    // Deselect second monomer-position pair
    mpViewer.modifyMutationCliffsSelection(secondPair, {shiftPressed: true, ctrlPressed: true});
    expect(mpViewer.mutationCliffsSelection[secondPair.positionOrClusterType].includes(secondPair.monomerOrCluster),
      false, `Monomer ${secondPair.monomerOrCluster} is still selected at position ` +
      `${secondPair.positionOrClusterType} after deselection`);
    expect(selection.trueCount, firstPair.mcCount, `Selection count is not equal to ${firstPair.mcCount} ` +
      `for monomer ${firstPair.monomerOrCluster} at position ${firstPair.positionOrClusterType}`);

    // Clear monomer-position selection
    mpViewer.mutationCliffsSelection = initSelection(mpViewer.positionColumns);
    for (const [position, selectedMonomers] of Object.entries(mpViewer.mutationCliffsSelection)) {
      expect(selectedMonomers.length, 0, `Selection is not empty for position ${position} after clearing ` +
        `monomer-position selection`);
    }
    expect(selection.trueCount, 0, `Selection count is not equal to 0 after clearing monomer-position selection`);

    const lstViewer = model.findViewer(VIEWER_TYPE.LOGO_SUMMARY_TABLE) as LogoSummaryTable | null;
    expect(lstViewer !== null, true, `Couldn't find Logo Summary Table viewer`);
    // Select first cluster
    lstViewer!.modifyClusterSelection(firstCluster);
    expect(lstViewer!.clusterSelection[firstCluster.positionOrClusterType].includes(firstCluster.monomerOrCluster),
      true, `Cluster ${firstCluster.monomerOrCluster} is not selected`);
    expect(selection.trueCount, firstCluster.count, `Selection count is not equal to ${firstCluster.count} for ` +
      `cluster ${firstCluster.monomerOrCluster}`);

    // Select second cluster
    lstViewer!.modifyClusterSelection(secondCluster, {shiftPressed: true, ctrlPressed: false});
    expect(lstViewer!.clusterSelection[secondCluster.positionOrClusterType].includes(secondCluster.monomerOrCluster),
      true, `Cluster ${secondCluster.monomerOrCluster} is not selected`);
    expect(selection.trueCount, firstCluster.count + secondCluster.count, `Selection count is not equal to ` +
      `${firstCluster.count + secondCluster.count} for cluster ${firstCluster.monomerOrCluster} and cluster 
      ${secondCluster.monomerOrCluster}`);

    // Deselect first cluster
    lstViewer!.modifyClusterSelection(firstCluster, {shiftPressed: true, ctrlPressed: true});
    expect(lstViewer!.clusterSelection[firstCluster.positionOrClusterType].includes(firstCluster.monomerOrCluster),
      false, `Cluster ${firstCluster.monomerOrCluster} is still selected after deselection`);
    expect(selection.trueCount, secondCluster.count, `Selection count is not equal to ${secondCluster.count} for ` +
      `cluster ${secondCluster.monomerOrCluster} after deselection of cluster ${firstCluster.monomerOrCluster}`);

    // Clear selection
    lstViewer!.initClusterSelection();
    expect(lstViewer!.isClusterSelectionEmpty, true, `Selection is not empty after clearing cluster selection`);
    expect(selection.trueCount, 0, `Selection count is not equal to 0 after clearing cluster selection`);
  });

  test('Invariant Map selection', async () => {
    const selection = model.df.selection;

    for (const [position, filteredMonomers] of Object.entries(model.webLogoSelection))
      expect(filteredMonomers.length, 0, `Filter is not empty for position ${position} after initialization`);

    // Select by second monomer-position pair
    model.modifyWebLogoSelection(secondPair);
    expect(model.webLogoSelection[secondPair.positionOrClusterType].includes(secondPair.monomerOrCluster), true,
      `Monomer ${secondPair.monomerOrCluster} is not filtered at position ${secondPair.positionOrClusterType}`);
    expect(selection.trueCount, secondPair.imCount, `Filter count is not equal to ${secondPair.imCount} ` +
      `for monomer ${secondPair.monomerOrCluster} at position ${secondPair.positionOrClusterType}`);

    // Select by first monomer-position pair
    model.modifyWebLogoSelection(firstPair, {shiftPressed: true, ctrlPressed: false});
    expect(model.webLogoSelection[firstPair.positionOrClusterType].includes(firstPair.monomerOrCluster), true,
      `Monomer ${firstPair.monomerOrCluster} is not filtered at position ${firstPair.positionOrClusterType}`);
    expect(selection.trueCount, secondPair.imCount, `Filter count is not equal to ${secondPair.imCount} ` +
      `for monomer ${firstPair.monomerOrCluster} at position ${firstPair.positionOrClusterType}`);

    // Deselect filter for second monomer-position pair
    model.modifyWebLogoSelection(secondPair, {shiftPressed: true, ctrlPressed: true});
    expect(model.webLogoSelection[secondPair.positionOrClusterType].includes(secondPair.monomerOrCluster),
      false, `Monomer ${secondPair.monomerOrCluster} is still filtered at position 
      ${secondPair.positionOrClusterType} after deselection`);
    expect(selection.trueCount, firstPair.imCount, `Filter count is not equal to ${firstPair.imCount} ` +
      `for monomer ${firstPair.monomerOrCluster} at position ${firstPair.positionOrClusterType} after deselection of ` +
      `monomer ${secondPair.monomerOrCluster} at position ${secondPair.positionOrClusterType}`);

    // Clear selection
    expect(model.positionColumns !== null, true, `Position columns are not initialized`);
    model.webLogoSelection = initSelection(model.positionColumns!);
    expect(selection.trueCount, 0, `Filter count is not equal to ${0} after clearing monomer-position filter`);

    for (const [position, filteredMonomers] of Object.entries(model.webLogoSelection)) {
      expect(filteredMonomers.length, 0, `Filter is not empty for position ${position} after clearing ` +
        `monomer-position filter`);
    }
  });
}, {clear: false});
