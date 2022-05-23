import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {PeptidesController} from '../peptides';
import * as C from '../utils/constants';

export async function substitutionsWidget(table: DG.DataFrame): Promise<DG.Widget> {
  const controller = await PeptidesController.getInstance(table);
  controller.init(table);
  const substTable = controller.getSubstitutions();

  if (!substTable)
    return new DG.Widget(ui.label('No substitution table generated'));

  const dfRowCount = substTable.rowCount;
  const aminoInputFrom = ui.stringInput('from', '');
  const aminoInputTo = ui.stringInput('to', '');
  const fromToMap: {[key: string]: DG.BitSet} = {};
  let aminoFrom = '';
  let aminoTo = '';
  const initialCol: DG.Column = substTable.getCol('Initial');
  const substitutedCol: DG.Column = substTable.getCol('Substituted');

  initialCol.semType = 'alignedSequenceDifference';
  initialCol.name = 'Substitution';
  const separator = initialCol.tags[C.TAGS.SEPARATOR];
  // substTable.columns.remove('Substituted');
  // const substCol = table.getCol('Substitution');

  substTable.filter.setAll(false);
  for (let i = 0; i < dfRowCount; ++i) {
    // const [from, to] = substCol!.get(i).split('#');
    const from = initialCol.get(i);
    const to = substitutedCol.get(i);
    const aminosFrom: [] = from.split(separator);
    const aminosTo: [] = to.split(separator);

    for (let j = 0; j < aminosFrom.length; ++j) {
      const aar: string = table.tags[C.TAGS.AAR];
      const pos = parseInt(table.tags[C.TAGS.POSITION]);
      const aarFrom = aminosFrom[j];
      const aarTo = aminosTo[j];
      if (j == pos && aarFrom != aminosTo[j] && (aarFrom == aar || aarTo == aar)) {
        const idx = `${getAmino(aarFrom)}#${getAmino(aarTo)}`;

        if (!(idx in fromToMap))
          fromToMap[idx] = DG.BitSet.create(dfRowCount);
        fromToMap[idx].set(i, true);
        substTable.filter.set(i, true);
      }
    }
  }

  for (let i = 0; i < initialCol.length; ++i) {
    const sequenceDifference = `${initialCol.get(i)}#${substitutedCol.get(i)}`;
    initialCol.set(i, sequenceDifference);
  }

  aminoInputFrom.onInput(() => {
    aminoFrom = getAmino(aminoInputFrom.value);
    const fromKey = `${aminoFrom}#${aminoTo}`;
    if (fromKey in fromToMap)
      substTable.selection.copyFrom(fromToMap[fromKey]);
  });

  aminoInputTo.onInput(() => {
    aminoTo = getAmino(aminoInputTo.value);
    const toKey = `${aminoFrom}#${aminoTo}`;
    if (toKey in fromToMap)
      substTable.selection.copyFrom(fromToMap[toKey]);
  });

  (substTable.columns as DG.ColumnList).remove('Substituted');
  const grid = substTable.plot.grid();
  grid.props.allowEdit = false;
  grid.root.style.width = 'auto';
  grid.root.style.height = '150px';
  return new DG.Widget(ui.divV([aminoInputFrom.root, aminoInputTo.root, grid.root]));
}

function getAmino(amino: string): string {
  return amino === '' ? '-' : amino;
}
