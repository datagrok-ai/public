import * as DG from 'datagrok-api/dg';

import {category, test, before, expect, delay, after} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import {CLUSTER_TYPE, PeptidesModel} from '../model';
import {startAnalysis} from '../widgets/peptides';
import {scaleActivity} from '../utils/misc';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {COLUMNS_NAMES, SCALING_METHODS} from '../utils/constants';
import {TEST_COLUMN_NAMES} from './utils';

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

  test('Tooltip', async () => {
    expect(model.showMonomerTooltip(firstPair.monomerOrCluster, 0, 0), true,
      `Couldn't structure for monomer ${firstPair.monomerOrCluster}`);
  }, {skipReason: 'Need to find a way to replace _package variable to call for Bio function with tests'});

  test('Visible columns', async () => {
    const gridCols = model.analysisView.grid.columns;
    const posCols = model.splitSeqDf.columns.names();
    for (let colIdx = 1; colIdx < gridCols.length; colIdx++) {
      const col = gridCols.byIndex(colIdx)!;
      const tableColName = col.column!.name;
      const expectedVisibility = posCols.includes(tableColName) || (tableColName === COLUMNS_NAMES.ACTIVITY_SCALED);
      expect(col.visible, expectedVisibility, `Column ${tableColName} is visible === ${col.visible} but should be ` +
        `${expectedVisibility}`);
    }
  });

  //TODO: split into separate tests for Mutation Cliffs and Logo Summary Table
  test('Mutation Cliffs selection', async () => {
    const selection = model.df.selection;

    for (const [position, selectedMonomers] of Object.entries(model.mutationCliffsSelection))
      expect(selectedMonomers.length, 0, `Selection is not empty for position ${position} after initialization`);

    // Select first monomer-position pair
    model.modifyMutationCliffsSelection(firstPair);
    await delay(500);
    expect(model.mutationCliffsSelection[firstPair.positionOrClusterType].includes(firstPair.monomerOrCluster), true,
      `Monomer ${firstPair.monomerOrCluster} is not selected at position ${firstPair.positionOrClusterType}`);
    expect(selection.trueCount, firstPair.mcCount, `Selection count is not equal to ${firstPair.mcCount} ` +
      `for monomer ${firstPair.monomerOrCluster} at position ${firstPair.positionOrClusterType}`);

    // Select second monomer-position pair
    model.modifyMutationCliffsSelection(secondPair, {shiftPressed: true, ctrlPressed: false});
    await delay(500);
    expect(model.mutationCliffsSelection[secondPair.positionOrClusterType].includes(secondPair.monomerOrCluster), true,
      `Monomer ${secondPair.monomerOrCluster} is not selected at position ${secondPair.positionOrClusterType}`);
    expect(selection.trueCount, secondPair.mcCount + firstPair.mcCount, `Selection count is not equal ` +
      `to ${secondPair.mcCount + firstPair.mcCount} for monomer ${secondPair.monomerOrCluster} at ` +
      `position ${secondPair.positionOrClusterType}`);

    // Deselect second monomer-position pair
    model.modifyMutationCliffsSelection(secondPair, {shiftPressed: true, ctrlPressed: true});
    await delay(500);
    expect(model.mutationCliffsSelection[secondPair.positionOrClusterType].includes(secondPair.monomerOrCluster), false,
      `Monomer ${secondPair.monomerOrCluster} is still selected at position ${secondPair.positionOrClusterType} after deselection`);
    expect(selection.trueCount, firstPair.mcCount, `Selection count is not equal to ${firstPair.mcCount} ` +
      `for monomer ${firstPair.monomerOrCluster} at position ${firstPair.positionOrClusterType}`);

    // Clear monomer-position selection
    model.initMutationCliffsSelection();
    await delay(500);
    for (const [position, selectedMonomers] of Object.entries(model.mutationCliffsSelection)) {
      expect(selectedMonomers.length, 0, `Selection is not empty for position ${position} after clearing ` +
        `monomer-position selection`);
    }
    expect(selection.trueCount, 0, `Selection count is not equal to 0 after clearing monomer-position selection`);

    // Select first cluster
    model.modifyClusterSelection(firstCluster);
    await delay(500);
    expect(model.clusterSelection[firstCluster.positionOrClusterType].includes(firstCluster.monomerOrCluster), true,
      `Cluster ${firstCluster.monomerOrCluster} is not selected`);
    expect(selection.trueCount, firstCluster.count, `Selection count is not equal to ${firstCluster.count} for ` +
      `cluster ${firstCluster.monomerOrCluster}`);

    // Select second cluster
    model.modifyClusterSelection(secondCluster, {shiftPressed: true, ctrlPressed: false});
    await delay(500);
    expect(model.clusterSelection[secondCluster.positionOrClusterType].includes(secondCluster.monomerOrCluster), true,
      `Cluster ${secondCluster.monomerOrCluster} is not selected`);
    expect(selection.trueCount, firstCluster.count + secondCluster.count, `Selection count is not equal to ` +
      `${firstCluster.count + secondCluster.count} for cluster ${firstCluster.monomerOrCluster} and cluster ${secondCluster.monomerOrCluster}`);

    // Deselect first cluster
    model.modifyClusterSelection(firstCluster, {shiftPressed: true, ctrlPressed: true});
    await delay(500);
    expect(model.clusterSelection[firstCluster.positionOrClusterType].includes(firstCluster.monomerOrCluster), false,
      `Cluster ${firstCluster.monomerOrCluster} is still selected after deselection`);
    expect(selection.trueCount, secondCluster.count, `Selection count is not equal to ${secondCluster.count} for ` +
      `cluster ${secondCluster.monomerOrCluster} after deselection of cluster ${firstCluster.monomerOrCluster}`);

    // Clear selection
    model.initClusterSelection();
    await delay(500);
    expect(model.isClusterSelectionEmpty, true, `Selection is not empty after clearing cluster selection`);
    expect(selection.trueCount, 0, `Selection count is not equal to 0 after clearing cluster selection`);
  });

  test('Invariant Map selection', async () => {
    const selection = model.df.selection;

    for (const [position, filteredMonomers] of Object.entries(model.invariantMapSelection))
      expect(filteredMonomers.length, 0, `Filter is not empty for position ${position} after initialization`);

    // Select by second monomer-position pair
    model.modifyInvariantMapSelection(secondPair);
    await delay(500);
    expect(model.invariantMapSelection[secondPair.positionOrClusterType].includes(secondPair.monomerOrCluster), true,
      `Monomer ${secondPair.monomerOrCluster} is not filtered at position ${secondPair.positionOrClusterType}`);
    expect(selection.trueCount, secondPair.imCount, `Filter count is not equal to ${secondPair.imCount} ` +
      `for monomer ${secondPair.monomerOrCluster} at position ${secondPair.positionOrClusterType}`);

    // Select by first monomer-position pair
    model.modifyInvariantMapSelection(firstPair, {shiftPressed: true, ctrlPressed: false});
    await delay(500);
    expect(model.invariantMapSelection[firstPair.positionOrClusterType].includes(firstPair.monomerOrCluster), true,
      `Monomer ${firstPair.monomerOrCluster} is not filtered at position ${firstPair.positionOrClusterType}`);
    expect(selection.trueCount, secondPair.imCount, `Filter count is not equal to ${secondPair.imCount} ` +
      `for monomer ${firstPair.monomerOrCluster} at position ${firstPair.positionOrClusterType}`);

    // Deselect filter for second monomer-position pair
    model.modifyInvariantMapSelection(secondPair, {shiftPressed: true, ctrlPressed: true});
    await delay(500);
    expect(model.invariantMapSelection[secondPair.positionOrClusterType].includes(secondPair.monomerOrCluster), false,
      `Monomer ${secondPair.monomerOrCluster} is still filtered at position ${secondPair.positionOrClusterType} after ` +
      `deselection`);
    expect(selection.trueCount, firstPair.imCount, `Filter count is not equal to ${firstPair.imCount} ` +
      `for monomer ${firstPair.monomerOrCluster} at position ${firstPair.positionOrClusterType} after deselection of ` +
      `monomer ${secondPair.monomerOrCluster} at position ${secondPair.positionOrClusterType}`);

    // Clear selection
    model.initInvariantMapSelection();
    await delay(500);
    expect(selection.trueCount, 0, `Filter count is not equal to ${0} after clearing monomer-position filter`);

    for (const [position, filteredMonomers] of Object.entries(model.invariantMapSelection)) {
      expect(filteredMonomers.length, 0, `Filter is not empty for position ${position} after clearing ` +
        `monomer-position filter`);
    }
  });
}, {clear: false});
