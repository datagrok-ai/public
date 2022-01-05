import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as OCL from 'openchemlib/full.js';
import $ from 'cash-dom';
import {isMolBlock} from "../chem-utils";
import {oclMol} from "./drug-likeness";

export function molfileWidget(molStr: string): DG.Widget {
  const molfileStr = isMolBlock(molStr) ? molStr : oclMol(molStr).toMolfile();

  const molfileInput = ui.textInput('', molfileStr);
  (molfileInput.input as HTMLElement).style.height = '300px';
  (molfileInput.input as HTMLElement).style.overflow = 'hidden';

  const clipboardBtn = ui.button(ui.iconFA('copy'), () => {
    navigator.clipboard.writeText(molfileInput.stringValue);
    const copyIcon = clipboardBtn.removeChild(clipboardBtn.firstChild!);
    clipboardBtn.appendChild(ui.iconFA('clipboard-check'));
    setTimeout(() => {
      clipboardBtn.removeChild(clipboardBtn.firstChild!);
      clipboardBtn.appendChild(copyIcon);
    }, 1000);
  }, 'Copy');
  clipboardBtn.style.right = '30px';
  $(clipboardBtn).addClass('dt-snippet-editor-icon');

  const resetBtn = ui.button(ui.iconFA('redo'), () => molfileInput.value = molfileStr, 'Reset');
  $(resetBtn).addClass('dt-snippet-editor-icon dt-reset-icon');

  return new DG.Widget(ui.divV([clipboardBtn, resetBtn, molfileInput.root], 'dt-textarea-box'));
}
