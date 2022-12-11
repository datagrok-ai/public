import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as C from '../utils/constants';
import * as type from '../utils/types';
import {PeptidesModel} from '../model';
import {getSeparator} from '../utils/misc';

export function mutationCliffsWidget(table: DG.DataFrame, model: PeptidesModel): DG.Widget {
  const substInfo = model.substitutionsInfo;
  const currentCell = model.mutationCliffsSelection;
  const positions = Object.keys(currentCell);

  if (!positions.length)
    return new DG.Widget(ui.label('No mutations table generated'));

  const substitutionsArray: string[] = [];
  const deltaArray: number[] = [];
  const substitutedToArray: string[] = [];
  const alignedSeqCol = table.getCol(model.settings.sequenceColumnName!);
  const activityScaledCol = table.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED);
  const seenIndexes = new Map<number, number[]>();

  for (const pos of positions) {
    for (const aar of currentCell[pos]) {
      const substitutionsMap = substInfo.get(aar)?.get(pos) as Map<number, type.UTypedArray> | undefined;
      if (typeof substitutionsMap === 'undefined')
        continue;

      const posCol = table.getCol(pos);
      for (const [referenceIdx, indexArray] of substitutionsMap.entries()) {
        const forbiddentIndexes = seenIndexes.has(referenceIdx) ? seenIndexes.get(referenceIdx)! : [];
        const baseSequence = alignedSeqCol.get(referenceIdx);
        const baseActivity = activityScaledCol.get(referenceIdx);

        for (const subIdx of indexArray) {
          if (forbiddentIndexes.includes(subIdx))
            continue;

          if (!seenIndexes.has(subIdx))
            seenIndexes.set(subIdx, []);

          seenIndexes.get(subIdx)!.push(referenceIdx);
          substitutionsArray.push(`${baseSequence}#${alignedSeqCol.get(subIdx)}`);
          deltaArray.push(baseActivity - activityScaledCol.get(subIdx));
          substitutedToArray.push(posCol.get(subIdx));
        }
      }
    }
  }

  if (!substitutionsArray.length)
    return new DG.Widget(ui.label('No mutations table generated'));

  const substCol = DG.Column.fromStrings('Mutation', substitutionsArray);
  substCol.semType = C.SEM_TYPES.MACROMOLECULE_DIFFERENCE;
  substCol.tags[C.TAGS.SEPARATOR] = getSeparator(alignedSeqCol);
  substCol.tags[DG.TAGS.UNITS] = alignedSeqCol.tags[DG.TAGS.UNITS];
  substCol.tags[DG.TAGS.CELL_RENDERER] = 'MacromoleculeDifference';
  const toColName = '~to';
  const hiddenSubstToAarCol = DG.Column.fromStrings(toColName, substitutedToArray);
  const substTable =
    DG.DataFrame.fromColumns([substCol, DG.Column.fromList('double', 'Delta', deltaArray), hiddenSubstToAarCol]);

  const aminoToInput = ui.stringInput('Mutated to:', '', () => {
    const substitutedToAar = aminoToInput.stringValue;
    if (substitutedToAar != '')
      substTable.filter.init((idx) => hiddenSubstToAarCol.get(idx) === substitutedToAar);
    else
      substTable.filter.setAll(true);
  });

  const grid = substTable.plot.grid();
  grid.props.allowEdit = false;
  const gridRoot = grid.root;
  gridRoot.style.width = 'auto';
  gridRoot.style.height = '150px';
  return new DG.Widget(ui.divV([aminoToInput.root, gridRoot]));
}
