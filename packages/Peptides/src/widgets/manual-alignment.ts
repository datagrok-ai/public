import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

export function manualAlignmentWidget(alignedSequenceCol: DG.Column, currentDf: DG.DataFrame) {
  //TODO: update viewers right when the changes get applied
  const sequenceInput = ui.textInput('', alignedSequenceCol.get(currentDf.currentRowIdx));
  (sequenceInput.input as HTMLElement).style.height = '50px';
  (sequenceInput.input as HTMLElement).style.overflow = 'hidden';

  const applyChangesBtn = ui.button('Apply', () => {
    alignedSequenceCol.set(currentDf.currentRowIdx, sequenceInput.value);
  });

  const resetBtn = ui.button(
    ui.iconFA('redo'),
    () => sequenceInput.value = alignedSequenceCol.get(currentDf.currentRowIdx),
    'Reset',
  );
  $(resetBtn).addClass('dt-snippet-editor-icon dt-reset-icon');

  return new DG.Widget(ui.divV([resetBtn, sequenceInput.root, applyChangesBtn], 'dt-textarea-box'));
}
