import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as C from '../utils/constants';
import * as type from '../utils/types';
import {PeptidesModel} from '../model';

export async function substitutionsWidget(table: DG.DataFrame): Promise<DG.Widget> {
  const model = await PeptidesModel.getInstance(table);
  const substInfo = model.substitutionsInfo;
  // const currentCell = model.getCurrentAARandPos();
  const currentCell = model.currentSelection;
  const pos = Object.keys(currentCell)[0];
  const aar = currentCell[pos][0];

  if (currentCell.aar === currentCell.pos)
    return new DG.Widget(ui.label('No substitution table generated'));

  const substitutionsMap =
    substInfo.get(aar)?.get(pos) as Map<number, type.UTypedArray> | undefined;
  if (typeof substitutionsMap === 'undefined')
    return new DG.Widget(ui.label('No substitution table generated'));

  const substitutionsArray: string[] = [];
  const deltaArray: number[] = [];
  const substitutedToArray: string[] = [];
  // const alignedSeqCol = table.getCol(C.COLUMNS_NAMES.ALIGNED_SEQUENCE);
  // const activityScaledCol = table.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED);
  const alignedSeqCol = table.columns.bySemType(C.SEM_TYPES.ALIGNED_SEQUENCE)!;
  const activityScaledCol = table.columns.bySemType(C.SEM_TYPES.ACTIVITY_SCALED)!;
  const posCol = table.getCol(pos);
  for (const [referenceIdx, indexArray] of substitutionsMap.entries()) {
    const baseSequence = alignedSeqCol.get(referenceIdx);
    const baseActivity = activityScaledCol.get(referenceIdx);
    for (const subIdx of indexArray) {
      substitutionsArray.push(`${baseSequence}#${alignedSeqCol.get(subIdx)}`);
      deltaArray.push(baseActivity - activityScaledCol.get(subIdx));
      substitutedToArray.push(posCol.get(subIdx));
    }
  }
  const substCol = DG.Column.fromStrings('Substiutions', substitutionsArray);
  substCol.semType = C.SEM_TYPES.ALIGNED_SEQUENCE_DIFFERENCE;
  const toColName = '~to';
  const hiddenSubstToAarCol = DG.Column.fromStrings(toColName, substitutedToArray);
  const substTable =
    DG.DataFrame.fromColumns([substCol, DG.Column.fromList('double', 'Delta', deltaArray), hiddenSubstToAarCol]);

  const aminoToInput = ui.stringInput('Substituted to:', '', () => {
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
