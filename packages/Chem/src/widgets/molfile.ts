import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import {isMolBlock, _convertMolNotation} from '../utils/convert-notation-utils';
import {oclMol} from '../utils/chem-common-ocl';
import '../../css/chem.css';
import {getRdKitModule} from '../utils/chem-common-rdkit';

export function getPanelElements(molStr: string): [HTMLButtonElement, HTMLButtonElement, DG.InputBase] {
  const molfileStr = isMolBlock(molStr) ? molStr : oclMol(molStr).toMolfile();

  const molfileInput = ui.textInput('', molfileStr);
  molfileInput.input.style.height = '300px';
  molfileInput.input.style.overflow = 'hidden';

  const clipboardBtn = ui.button(ui.iconFA('copy'), async () => {
    await navigator.clipboard.writeText(molfileInput.stringValue);
    const copyIcon = clipboardBtn.removeChild(clipboardBtn.firstChild!);
    clipboardBtn.appendChild(ui.iconFA('clipboard-check'));
    setTimeout(() => {
      clipboardBtn.removeChild(clipboardBtn.firstChild!);
      clipboardBtn.appendChild(copyIcon);
    }, 1000);
  }, 'Copy');
  clipboardBtn.style.right = '30px';
  $(clipboardBtn).addClass('chem-snippet-editor-icon');

  const resetBtn = ui.button(ui.iconFA('redo'), () => molfileInput.value = molfileStr, 'Reset');
  $(resetBtn).addClass('chem-snippet-editor-icon dt-reset-icon');

  return [clipboardBtn, resetBtn, molfileInput];
}

export function molfileWidget(molStr: string): DG.Widget {
  const rdKitModule = getRdKitModule();
  try {
    molStr = _convertMolNotation(molStr, 'unknown', 'molblock', rdKitModule);
  } catch (e) {
    return new DG.Widget(ui.divText('Molecule is possible malformed'));
  }
  const panelElements: any[] = getPanelElements(molStr);
  panelElements[2] = panelElements[2].root;
  return new DG.Widget(ui.divV(panelElements, 'chem-textarea-box'));
}
