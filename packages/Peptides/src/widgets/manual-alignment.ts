import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import '../styles.css';
import {PeptidesModel} from '../model';

/**
 * Allows to edit sequence and apply changes to the table and analysis.
 * @param alignedSequenceCol Aligned sequence column
 * @param currentDf Working table
 * @return Widget for manual sequence alignment
 */
export function manualAlignmentWidget(alignedSequenceCol: DG.Column<string>, currentDf: DG.DataFrame): DG.Widget {
  const rowIndex = currentDf.currentRowIdx;
  const columnName = currentDf.currentCol?.name;
  const sequenceInput = ui.input.textArea('Sequence', {value: alignedSequenceCol.get(rowIndex)!});
  $(sequenceInput.root).addClass('pep-textinput');

  const applyChangesBtn = ui.button('Apply', async () => {
    applyChangesBtn.disabled = true;
    try {
      await PeptidesModel.getInstance(currentDf).updateSequence(rowIndex, sequenceInput.value, columnName);
    } finally {
      applyChangesBtn.disabled = false;
    }
  }, 'Apply changes');

  const resetBtn = ui.button(ui.iconFA('redo'),
    () => sequenceInput.value = alignedSequenceCol.get(rowIndex)!, 'Reset');
  $(resetBtn).addClass('pep-snippet-editor-icon pep-reset-icon');

  return new DG.Widget(ui.divV([resetBtn, sequenceInput.root, applyChangesBtn], 'pep-textarea-box'));
}
