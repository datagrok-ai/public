import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import '../styles.css';
import {PeptidesModel} from '../model';
import {PeptideUtils} from '../peptideUtils';

/**
 * Allows to edit sequence and apply changes to the table and analysis.
 * @param alignedSequenceCol Aligned sequence column
 * @param currentDf Working table
 * @return Widget for manual sequence alignment
 */
export function manualAlignmentWidget(alignedSequenceCol: DG.Column<string>, currentDf: DG.DataFrame): DG.Widget {
  const sequenceInput = ui.input.textArea('', {value: alignedSequenceCol.get(currentDf.currentRowIdx)!});
  $(sequenceInput.root).addClass('pep-textinput');

  const applyChangesBtn = ui.button('Apply', async () => {
    const sh = PeptideUtils.getSeqHelper().getSeqHandler(alignedSequenceCol);
    const newSequence = sequenceInput.value;
    const splitSequence = sh.splitter(newSequence);
    const affectedRowIndex = currentDf.currentRowIdx;
    alignedSequenceCol.set(affectedRowIndex, newSequence);
    for (let i = 0; i < splitSequence.length; i++) {
      const part = splitSequence.getCanonical(i);
      if (currentDf.col(i.toString()) !== null)
        currentDf.set(i.toString(), affectedRowIndex, part);
    }
    const temp = grok.shell.o;
    grok.shell.o = null;

    const peptidesController = PeptidesModel.getInstance(currentDf);
    peptidesController.updateGrid();

    setTimeout(() => {
      grok.shell.o = temp;
    }, 100);
  }, 'Apply changes');

  const resetBtn = ui.button(ui.iconFA('redo'),
    () => sequenceInput.value = alignedSequenceCol.get(currentDf.currentRowIdx)!, 'Reset');
  $(resetBtn).addClass('pep-snippet-editor-icon pep-reset-icon');

  return new DG.Widget(ui.divV([resetBtn, sequenceInput.root, applyChangesBtn], 'pep-textarea-box'));
}
