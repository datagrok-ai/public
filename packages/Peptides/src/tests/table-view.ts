import * as DG from 'datagrok-api/dg';

import {category, test, before, expect} from '@datagrok-libraries/utils/src/test';
import {_package} from '../package-test';
import {PeptidesModel} from '../model';
import {startAnalysis} from '../widgets/peptides';
import {ScalingMethods} from '../utils/types';
import {scaleActivity} from '../utils/misc';
import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {COLUMNS_NAMES} from '../utils/constants';

category('Table view', () => {
  let df: DG.DataFrame;
  let model: PeptidesModel;
  let activityCol: DG.Column<number>;
  let sequenceCol: DG.Column<string>;
  let clusterCol: DG.Column<any>;
  let scaledActivityCol: DG.Column<number>;
  const scaling: ScalingMethods = 'none';

  const firstPair = {monomer: 'N', position: '4', count: 7};
  const secondPair = {monomer: 'meI', position: '1', count: 10};

  before(async () => {
    df = DG.DataFrame.fromCsv(await _package.files.readAsText('tests/HELM_small.csv'));
    activityCol = df.getCol('activity');
    sequenceCol = df.getCol('sequence');
    sequenceCol.semType = DG.SEMTYPE.MACROMOLECULE;
    sequenceCol.setTag(DG.TAGS.UNITS, NOTATION.HELM);
    scaledActivityCol = scaleActivity(activityCol, scaling);
    clusterCol = df.getCol('cluster');
    const tempModel = await startAnalysis(activityCol, sequenceCol, clusterCol, df, scaledActivityCol, scaling);
    if (tempModel === null)
      throw new Error('Model is null');
    model = tempModel;
  });

  test('Tooltip', async () => {

  }, {skipReason: 'Not implemented yet'});

  test('Visible columns', async () => {
    const gridCols = model.analysisView.grid.columns;
    const posCols = model.splitSeqDf.columns.names();
    const visibleColumns = Object.keys(model.settings.columns ?? {});
    for (let colIdx = 1; colIdx < gridCols.length; colIdx++) {
      const col = gridCols.byIndex(colIdx)!;
      const tableColName = col.column!.name;
      expect(col.visible, posCols.includes(tableColName) || (tableColName === COLUMNS_NAMES.ACTIVITY_SCALED) ||
        visibleColumns.includes(tableColName));
    }
  });

  test('Selection', async () => {
    const selection = model.df.selection;

    // Select first monomer-position pair
    model.modifyMonomerPositionSelection(firstPair.monomer, firstPair.position, false);
    expect(model.mutationCliffsSelection[firstPair.position].includes(firstPair.monomer), true);
    expect(selection.trueCount, firstPair.count);

    // Select second monomer-position pair
    model.modifyMonomerPositionSelection(secondPair.monomer, secondPair.position, false);
    expect(model.mutationCliffsSelection[secondPair.position].includes(secondPair.monomer), true);
    expect(selection.trueCount, secondPair.count);

    // Deselect second monomer-position pair
    model.modifyMonomerPositionSelection(secondPair.monomer, secondPair.position, false);
    expect(model.mutationCliffsSelection[secondPair.position].includes(secondPair.monomer), false);
    expect(selection.trueCount, firstPair.count);

    // Clear selection
    model.initMutationCliffsSelection(true);
    expect(selection.trueCount, 0);
  });

  test('Filtering', async () => {
    const filter = model.df.filter;

    // Select second monomer-position pair
    model.modifyMonomerPositionSelection(secondPair.monomer, secondPair.position, true);
    expect(model.invariantMapSelection[secondPair.position].includes(secondPair.monomer), true);
    expect(filter.trueCount, secondPair.count);

    // Select first monomer-position pair
    model.modifyMonomerPositionSelection(firstPair.monomer, firstPair.position, true);
    expect(model.invariantMapSelection[firstPair.position].includes(firstPair.monomer), true);
    expect(filter.trueCount, firstPair.count);

    // Deselect second monomer-position pair
    model.modifyMonomerPositionSelection(secondPair.monomer, secondPair.position, true);
    expect(model.invariantMapSelection[secondPair.position].includes(secondPair.monomer), false);
    expect(filter.trueCount, firstPair.count);

    // Clear selection
    model.initInvariantMapSelection(true);
    expect(filter.trueCount, df.rowCount);
  });
});
