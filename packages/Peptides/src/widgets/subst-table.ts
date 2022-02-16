import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export function substTableWidget(table: DG.DataFrame): DG.Widget {
  if (!table?.temp['isReal'])
    return new DG.Widget(ui.label('No substitution'));

  const dfRowCount = table.rowCount;
  const aminoInputFrom = ui.stringInput('from', '');
  const aminoInputTo = ui.stringInput('to', '');
  const fromToMap: {[key: string]: DG.BitSet} = {};
  let aminoFrom = '';
  let aminoTo = '';
  const initialCol: DG.Column = table.columns.byName('Initial');
  const substitutedCol: DG.Column = table.columns.byName('Substituted');

  // for (let i = 0; i < initialCol.length; ++i) {
  //   const initialPeptide: string = initialCol.get(i);
  //   const substPeptide: string = substitutedCol.get(i);

  //   initialCol.set(i, initialPeptide + '#' + substPeptide);
  // }

  initialCol.semType = 'alignedSequenceDifference';
  initialCol.name = 'Substitution';
  table.columns.remove('Substituted');
  // const substCol = table.getCol('Substitution');

  for (let i = 0; i < dfRowCount; ++i) {
    // const [from, to] = substCol!.get(i).split('#');
    const from = initialCol.get(i);
    const to = substitutedCol.get(i);
    const aminosFrom: [] = from.split('-');
    const aminosTo: [] = to.split('-');

    for (let j = 0; j < aminosFrom.length; ++j) {
      if (aminosFrom[j] != aminosTo[j]) {
        const idx = (aminosFrom[j] === '' ? '-' : aminosFrom[j]) + '#' + (aminosTo[j] === '' ? '-' : aminosTo[j]);

        if (!(idx in fromToMap))
          fromToMap[idx] = DG.BitSet.create(dfRowCount);
        fromToMap[idx].set(i, true);
      }
    }
  }

  aminoInputFrom.onInput(() => {
    aminoFrom = aminoInputFrom.value;
    const fromKey = aminoFrom + '#' + aminoTo;
    if (fromKey in fromToMap)
      table.selection.copyFrom(fromToMap[fromKey]);
  });

  aminoInputTo.onInput(() => {
    aminoTo = aminoInputTo.value;
    const toKey = aminoFrom + '#' + aminoTo;
    if (toKey in fromToMap)
      table.selection.copyFrom(fromToMap[toKey]);
  });

  const grid = table.plot.grid();
  grid.props.allowEdit = false;
  grid.root.style.width = 'auto';
  return new DG.Widget(ui.divV([aminoInputFrom.root, aminoInputTo.root, grid.root]));
}
