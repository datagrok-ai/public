import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, test, before, expect, awaitCheck} from '@datagrok-libraries/utils/src/test';
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

  const firstPair = {monomer: 'Aze', position: '10', mcCount: 3, imCount: 1};
  const secondPair = {monomer: 'meI', position: '1', mcCount: 2, imCount: 10};
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

    // Ensure grid finished initializing to prevent Unhandled exceptions
    let accrodionInit = false;
    grok.events.onAccordionConstructed.subscribe((_) => accrodionInit = true);
    await awaitCheck(() => model!.df.currentRowIdx === 0, 'Grid cell never finished initializing', 2000);
    await awaitCheck(() => grok.shell.o instanceof DG.Column, 'Shell object never changed', 2000);
    await awaitCheck(() => accrodionInit, 'Accordion never finished initializing', 2000);
  });

  test('Tooltip', async () => {
    expect(model.showMonomerTooltip(firstPair.monomer, 0, 0), true,
      `Couldn't structure for monomer ${firstPair.monomer}`);
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
    model.modifyMonomerPositionSelection(firstPair.monomer, firstPair.position, false);
    expect(model.mutationCliffsSelection[firstPair.position].includes(firstPair.monomer), true,
      `Monomer ${firstPair.monomer} is not selected at position ${firstPair.position}`);
    expect(selection.trueCount, firstPair.mcCount, `Selection count is not equal to ${firstPair.mcCount} ` +
      `for monomer ${firstPair.monomer} at position ${firstPair.position}`);

    // Select second monomer-position pair
    model.modifyMonomerPositionSelection(secondPair.monomer, secondPair.position, false);
    expect(model.mutationCliffsSelection[secondPair.position].includes(secondPair.monomer), true,
      `Monomer ${secondPair.monomer} is not selected at position ${secondPair.position}`);
    expect(selection.trueCount, secondPair.mcCount + firstPair.mcCount, `Selection count is not equal ` +
      `to ${secondPair.mcCount + firstPair.mcCount} for monomer ${secondPair.monomer} at ` +
      `position ${secondPair.position}`);

    // Deselect second monomer-position pair
    model.modifyMonomerPositionSelection(secondPair.monomer, secondPair.position, false);
    expect(model.mutationCliffsSelection[secondPair.position].includes(secondPair.monomer), false,
      `Monomer ${secondPair.monomer} is still selected at position ${secondPair.position} after deselection`);
    expect(selection.trueCount, firstPair.mcCount, `Selection count is not equal to ${firstPair.mcCount} ` +
      `for monomer ${firstPair.monomer} at position ${firstPair.position}`);

    // Clear monomer-position selection
    model.initMutationCliffsSelection({cleanInit: true});
    for (const [position, selectedMonomers] of Object.entries(model.mutationCliffsSelection)) {
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

  test('Invariant Map selection', async () => {
    const selection = model.df.selection;

    for (const [position, filteredMonomers] of Object.entries(model.invariantMapSelection))
      expect(filteredMonomers.length, 0, `Filter is not empty for position ${position} after initialization`);

    // Select by second monomer-position pair
    model.modifyMonomerPositionSelection(secondPair.monomer, secondPair.position, true);
    expect(model.invariantMapSelection[secondPair.position].includes(secondPair.monomer), true,
      `Monomer ${secondPair.monomer} is not filtered at position ${secondPair.position}`);
    expect(selection.trueCount, secondPair.imCount, `Filter count is not equal to ${secondPair.imCount} ` +
      `for monomer ${secondPair.monomer} at position ${secondPair.position}`);

    // Select by first monomer-position pair
    model.modifyMonomerPositionSelection(firstPair.monomer, firstPair.position, true);
    expect(model.invariantMapSelection[firstPair.position].includes(firstPair.monomer), true,
      `Monomer ${firstPair.monomer} is not filtered at position ${firstPair.position}`);
    expect(selection.trueCount, secondPair.imCount, `Filter count is not equal to ${secondPair.imCount} ` +
      `for monomer ${firstPair.monomer} at position ${firstPair.position}`);

    // Deselect filter for second monomer-position pair
    model.modifyMonomerPositionSelection(secondPair.monomer, secondPair.position, true);
    expect(model.invariantMapSelection[secondPair.position].includes(secondPair.monomer), false,
      `Monomer ${secondPair.monomer} is still filtered at position ${secondPair.position} after ` +
      `deselection`);
    expect(selection.trueCount, firstPair.imCount, `Filter count is not equal to ${firstPair.imCount} ` +
      `for monomer ${firstPair.monomer} at position ${firstPair.position} after deselection of ` +
      `monomer ${secondPair.monomer} at position ${secondPair.position}`);

    // Clear selection
    model.initInvariantMapSelection({cleanInit: true});
    expect(selection.trueCount, 0, `Filter count is not equal to ${0} after clearing monomer-position filter`);

    for (const [position, filteredMonomers] of Object.entries(model.invariantMapSelection)) {
      expect(filteredMonomers.length, 0, `Filter is not empty for position ${position} after clearing ` +
        `monomer-position filter`);
    }
  });
}, {clear: false});
