import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import {model} from '../viewers/model';
import {splitAlignedPeptides} from '../utils/split-aligned';

/**
 * Manual sequence alignment widget.
 * 
 * @param {DG.Column} alignedSequenceCol Aligned sequence column.
 * @param {DG.DataFrame} currentDf Working table.
 * @returns {DG.Widget} Widget for manual sequence alignment.
 */
export function manualAlignmentWidget(alignedSequenceCol: DG.Column, currentDf: DG.DataFrame) {
  const sequenceInput = ui.textInput('', alignedSequenceCol.get(currentDf.currentRowIdx));
  (sequenceInput.input as HTMLElement).style.height = '50px';
  (sequenceInput.input as HTMLElement).style.overflow = 'hidden';

  const applyChangesBtn = ui.button('Apply', async () => {
    const newSequence = sequenceInput.value;
    const affectedRowIndex = currentDf.currentRowIdx;
    const [splitSequence] = splitAlignedPeptides(DG.Column.fromStrings('splitSequence', [newSequence]), false);

    alignedSequenceCol.set(affectedRowIndex, newSequence);
    for (const part of splitSequence.columns) {
      if (currentDf.col(part.name) !== null) {
        currentDf.set(part.name, affectedRowIndex, part.get(0));
      }
    }

    await model.updateDefault();
  });

  const resetBtn = ui.button(
    ui.iconFA('redo'),
    () => sequenceInput.value = alignedSequenceCol.get(currentDf.currentRowIdx),
    'Reset',
  );
  $(resetBtn).addClass('dt-snippet-editor-icon dt-reset-icon');

  return new DG.Widget(ui.divV([resetBtn, sequenceInput.root, applyChangesBtn], 'dt-textarea-box'));
}
