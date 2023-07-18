import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as C from '../utils/constants';
import * as type from '../utils/types';
import {PeptidesModel} from '../model';
import {getSeparator} from '../utils/misc';

export function mutationCliffsWidget(table: DG.DataFrame, model: PeptidesModel): DG.Widget {
  const filteredIndexes = table.filter.getSelectedIndexes();
  const substInfo = model.mutationCliffs;
  const currentCell = model.monomerPositionSelection;
  const positions = Object.keys(currentCell);

  if (!positions.length || substInfo === null)
    return new DG.Widget(ui.label('No mutations table generated'));

  const substitutionsArray: string[] = [];
  const deltaArray: number[] = [];
  const substitutedToArray: string[] = [];
  const alignedSeqCol = table.getCol(model.settings.sequenceColumnName!);
  const alignedSeqColCategories = alignedSeqCol.categories;
  const alignedSeqColData = alignedSeqCol.getRawData();
  const activityScaledCol = table.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED);
  const activityScaledColData = activityScaledCol.getRawData();
  const seenIndexes = new Map<number, number[]>();

  for (const pos of positions) {
    const posCol = table.getCol(pos);
    const posColCategories = posCol.categories;
    const posColData = posCol.getRawData();

    for (const aar of currentCell[pos]) {
      const substitutionsMap = substInfo.get(aar)?.get(pos) as Map<number, type.UTypedArray> | undefined;
      if (typeof substitutionsMap === 'undefined')
        continue;

      for (const [referenceIdx, indexArray] of substitutionsMap.entries()) {
        if (!filteredIndexes.includes(referenceIdx))
          continue;

        const forbiddentIndexes = seenIndexes.get(referenceIdx) ?? [];
        const baseSequence = alignedSeqColCategories[alignedSeqColData[referenceIdx]];
        const baseActivity = activityScaledColData[referenceIdx];

        for (const subIdx of indexArray) {
          if (forbiddentIndexes.includes(subIdx) || !filteredIndexes.includes(subIdx))
            continue;

          if (!seenIndexes.has(subIdx))
            seenIndexes.set(subIdx, []);
          const subSeq = alignedSeqColCategories[alignedSeqColData[subIdx]];

          seenIndexes.get(subIdx)!.push(referenceIdx);
          substitutionsArray.push(`${baseSequence}#${subSeq}`);
          deltaArray.push(baseActivity - activityScaledColData[subIdx]);
          substitutedToArray.push(posColCategories[posColData[subIdx]]);
        }
      }
    }
  }

  if (substitutionsArray.length === 0)
    return new DG.Widget(ui.label('No mutations table generated'));

  const substCol = DG.Column.fromStrings('Mutation', substitutionsArray);
  const toColName = '~to';
  const hiddenSubstToAarCol = DG.Column.fromStrings(toColName, substitutedToArray);
  const substTable =
    DG.DataFrame.fromColumns([substCol, DG.Column.fromList('double', 'Delta', deltaArray), hiddenSubstToAarCol]);

  const aminoToInput = ui.stringInput('Mutated to:', '', () => {
    const substitutedToAar = aminoToInput.stringValue;
    if (substitutedToAar !== '')
      substTable.filter.init((idx) => hiddenSubstToAarCol.get(idx) === substitutedToAar);
    else
      substTable.filter.setAll(true);
  });

  const grid = substTable.plot.grid();
  grid.props.allowEdit = false;
  const gridRoot = grid.root;
  gridRoot.style.width = 'auto';
  gridRoot.style.height = '150px';
  substCol.semType = C.SEM_TYPES.MACROMOLECULE_DIFFERENCE;
  substCol.tags[C.TAGS.SEPARATOR] = getSeparator(alignedSeqCol);
  substCol.tags[DG.TAGS.UNITS] = alignedSeqCol.tags[DG.TAGS.UNITS];
  substCol.tags[DG.TAGS.CELL_RENDERER] = 'MacromoleculeDifference';
  return new DG.Widget(ui.divV([aminoToInput.root, gridRoot]));
}
