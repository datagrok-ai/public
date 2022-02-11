import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export function substTableWidget(table: DG.DataFrame): DG.Widget {
  if (!table)
    return new DG.Widget(ui.label('No substitution'));

  const aminoInputFrom = ui.stringInput('from', '');
  const aminoInputTo = ui.stringInput('to', '');
  const fromToMap: {[key: string]: DG.BitSet} = {};
  let aminoFrom = '';
  let aminoTo = '';
  const initialCol: DG.Column = table.columns.byName('Initial');
  const substitutedCol: DG.Column = table.columns.byName('Substituted');

  for (let i = 0; i < initialCol.length; ++i) {
    const initialPeptide: string = initialCol.get(i);
    const substPeptide: string = substitutedCol.get(i);

    initialCol.set(i, initialPeptide + '#' + substPeptide);
  }

  initialCol.semType = 'alignedSequenceDifference';
  initialCol.name = 'Substitution';
  table.columns.remove('Substituted');
  const substCol = table.getCol('Substitution');

  for (let i = 0; i < substCol!.length; ++i) {
    const [from, to] = substCol!.get(i).split('#');
    const aminosFrom: [] = from.split('-');
    const aminosTo: [] = to.split('-');

    for (let j = 0; j < aminosFrom.length; ++j) {
      if (aminosFrom[j] != aminosTo[j]) {
        const idx = (aminosFrom[j] === '' ? '-' : aminosFrom[j]) + '#' + (aminosTo[j] === '' ? '-' : aminosTo[j]);

        if (!(idx in fromToMap))
          fromToMap[idx] = DG.BitSet.create(substCol.length);
        fromToMap[idx].set(i, true);
      }
    }
  }

  aminoInputFrom.onInput(() => {
    aminoFrom = aminoInputFrom.value;
    if (aminoFrom + '#' + aminoTo in fromToMap)
      table.selection.copyFrom(fromToMap[aminoFrom + '#' + aminoTo]);
  });

  aminoInputTo.onInput(() => {
    aminoTo = aminoInputTo.value;
    if (aminoFrom + '#' + aminoTo in fromToMap)
      table.selection.copyFrom(fromToMap[aminoFrom + '#' + aminoTo]);
  });

  const grid = table.plot.grid().root;
  grid.style.width = 'auto';
  return new DG.Widget(ui.divV([aminoInputFrom.root, aminoInputTo.root, grid]));
}
