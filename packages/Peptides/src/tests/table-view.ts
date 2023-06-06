import * as DG from 'datagrok-api/dg';

import {category, test, before, expect} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import {PeptidesModel} from '../model';
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

  const firstMonomerPair = {monomer: 'N', position: '4', count: 7};
  const secondMonomerPair = {monomer: 'meI', position: '1', count: 10};
  const firstCluster = {name: '0', count: 3};
  const secondCluster = {name: '1', count: 3};

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
  });

  test('Tooltip', async () => {
    expect(model.showMonomerTooltip(firstMonomerPair.monomer, 0, 0), true,
      `Couldn't structure for monomer ${firstMonomerPair.monomer}`);
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

  test('Selection', async () => {
    const selection = model.df.selection;

    for (const [position, selectedMonomers] of Object.entries(model.monomerPositionSelection))
      expect(selectedMonomers.length, 0, `Selection is not empty for position ${position} after initialization`);

    // Select first monomer-position pair
    model.modifyMonomerPositionSelection(firstMonomerPair.monomer, firstMonomerPair.position, false);
    expect(model.monomerPositionSelection[firstMonomerPair.position].includes(firstMonomerPair.monomer), true,
      `Monomer ${firstMonomerPair.monomer} is not selected at position ${firstMonomerPair.position}`);
    expect(selection.trueCount, firstMonomerPair.count, `Selection count is not equal to ${firstMonomerPair.count} ` +
      `for monomer ${firstMonomerPair.monomer} at position ${firstMonomerPair.position}`);

    // Select second monomer-position pair
    model.modifyMonomerPositionSelection(secondMonomerPair.monomer, secondMonomerPair.position, false);
    expect(model.monomerPositionSelection[secondMonomerPair.position].includes(secondMonomerPair.monomer), true,
      `Monomer ${secondMonomerPair.monomer} is not selected at position ${secondMonomerPair.position}`);
    expect(selection.trueCount, secondMonomerPair.count, `Selection count is not equal to ${secondMonomerPair.count} ` +
      `for monomer ${secondMonomerPair.monomer} at position ${secondMonomerPair.position}`);

    // Deselect second monomer-position pair
    model.modifyMonomerPositionSelection(secondMonomerPair.monomer, secondMonomerPair.position, false);
    expect(model.monomerPositionSelection[secondMonomerPair.position].includes(secondMonomerPair.monomer), false,
      `Monomer ${secondMonomerPair.monomer} is still selected at position ${secondMonomerPair.position} after ` +
      `deselection`);
    expect(selection.trueCount, firstMonomerPair.count, `Selection count is not equal to ${firstMonomerPair.count} ` +
      `for monomer ${firstMonomerPair.monomer} at position ${firstMonomerPair.position}`);

    // Clear monomer-position selection
    model.initMonomerPositionSelection({cleanInit: true});
    for (const [position, selectedMonomers] of Object.entries(model.monomerPositionSelection)) {
      expect(selectedMonomers.length, 0, `Selection is not empty for position ${position} after clearing ` +
        `monomer-position selection`);
    }
    expect(selection.trueCount, 0, `Selection count is not equal to 0 after clearing monomer-position selection`);

    // Select first cluster
    model.modifyClusterSelection(firstCluster.name);
    expect(model.clusterSelection.includes(firstCluster.name), true,
      `Cluster ${firstCluster.name} is not selected`);
    expect(selection.trueCount, firstCluster.count, `Selection count is not equal to ${firstCluster.count} for ` +
      `cluster ${firstCluster.name}`);

    // Select second cluster
    model.modifyClusterSelection(secondCluster.name);
    expect(model.clusterSelection.includes(secondCluster.name), true,
      `Cluster ${secondCluster.name} is not selected`);
    expect(selection.trueCount, firstCluster.count + secondCluster.count, `Selection count is not equal to ` +
      `${firstCluster.count + secondCluster.count} for cluster ${firstCluster.name} and cluster ${secondCluster.name}`);

    // Deselect first cluster
    model.modifyClusterSelection(firstCluster.name);
    expect(model.clusterSelection.includes(firstCluster.name), false,
      `Cluster ${firstCluster.name} is still selected after deselection`);
    expect(selection.trueCount, secondCluster.count, `Selection count is not equal to ${secondCluster.count} for ` +
      `cluster ${secondCluster.name} after deselection of cluster ${firstCluster.name}`);

    // Clear selection
    model.initClusterSelection();
    expect(model.clusterSelection.length, 0, `Selection is not empty after clearing cluster selection`);
    expect(selection.trueCount, 0, `Selection count is not equal to 0 after clearing cluster selection`);
  });

  test('Filtering', async () => {
    const filter = model.df.filter;

    for (const [position, filteredMonomers] of Object.entries(model.monomerPositionFilter))
      expect(filteredMonomers.length, 0, `Filter is not empty for position ${position} after initialization`);

    // Filter by second monomer-position pair
    model.modifyMonomerPositionSelection(secondMonomerPair.monomer, secondMonomerPair.position, true);
    expect(model.monomerPositionFilter[secondMonomerPair.position].includes(secondMonomerPair.monomer), true,
      `Monomer ${secondMonomerPair.monomer} is not filtered at position ${secondMonomerPair.position}`);
    expect(filter.trueCount, secondMonomerPair.count, `Filter count is not equal to ${secondMonomerPair.count} ` +
      `for monomer ${secondMonomerPair.monomer} at position ${secondMonomerPair.position}`);

    // Filter by first monomer-position pair
    model.modifyMonomerPositionSelection(firstMonomerPair.monomer, firstMonomerPair.position, true);
    expect(model.monomerPositionFilter[firstMonomerPair.position].includes(firstMonomerPair.monomer), true,
      `Monomer ${firstMonomerPair.monomer} is not filtered at position ${firstMonomerPair.position}`);
    expect(filter.trueCount, firstMonomerPair.count, `Filter count is not equal to ${firstMonomerPair.count} ` +
      `for monomer ${firstMonomerPair.monomer} at position ${firstMonomerPair.position}`);

    // Deselect filter for second monomer-position pair
    model.modifyMonomerPositionSelection(secondMonomerPair.monomer, secondMonomerPair.position, true);
    expect(model.monomerPositionFilter[secondMonomerPair.position].includes(secondMonomerPair.monomer), false,
      `Monomer ${secondMonomerPair.monomer} is still filtered at position ${secondMonomerPair.position} after ` +
      `deselection`);
    expect(filter.trueCount, firstMonomerPair.count, `Filter count is not equal to ${firstMonomerPair.count} ` +
      `for monomer ${firstMonomerPair.monomer} at position ${firstMonomerPair.position} after deselection of ` +
      `monomer ${secondMonomerPair.monomer} at position ${secondMonomerPair.position}`);

    // Clear selection
    model.initMonomerPositionFilter({cleanInit: true});
    expect(filter.trueCount, df.rowCount, `Filter count is not equal to ${df.rowCount} after clearing ` +
      `monomer-position filter`);

    for (const [position, filteredMonomers] of Object.entries(model.monomerPositionFilter)) {
      expect(filteredMonomers.length, 0, `Filter is not empty for position ${position} after clearing ` +
        `monomer-position filter`);
    }
  });
});
