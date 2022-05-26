import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {PeptidesController} from '../peptides';
import * as C from '../utils/constants';
import * as type from '../utils/types';

export async function substitutionsWidget(table: DG.DataFrame): Promise<DG.Widget> {
  const controller = await PeptidesController.getInstance(table);
  controller.init(table);
  const substInfo = controller.getSubstitutions();
  const currentCell = controller.getCurrentAARandPos();
  // const substTable: DG.DataFrame = DG.DataFrame.create();

  if (currentCell.aar === currentCell.pos)
    return new DG.Widget(ui.label('No substitution table generated'));

  const substitutionsMap = substInfo.get(currentCell.aar!)!.get(currentCell.pos!)! as Map<number, type.UTypedArray>;
  const substitutionsArray: string[] = [];
  const deltaArray: number[] = [];
  const substitutedToArray: string[] = [];
  const alignedSeqCol = table.getCol(C.COLUMNS_NAMES.ALIGNED_SEQUENCE);
  const activityScaledCol = table.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED);
  const posCol = table.getCol(currentCell.pos!);
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

  const grid = substTable.plot.grid();
  const aminoToInput = ui.stringInput('Substituted to:', '', () => {
    const substitutedToAar = aminoToInput.stringValue;
    if (substitutedToAar != '')
      substTable.filter.init((idx) => hiddenSubstToAarCol.get(idx) === substitutedToAar, false);
    else
      substTable.filter.setAll(true);
    grid.invalidate();
  });

  grid.props.allowEdit = false;
  const gridRoot = grid.root;
  gridRoot.style.width = 'auto';
  gridRoot.style.height = '150px';
  return new DG.Widget(ui.divV([aminoToInput.root, gridRoot]));
}
