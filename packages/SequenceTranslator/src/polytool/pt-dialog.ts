/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import wu from 'wu';
import {Unsubscribable} from 'rxjs';

import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule/consts';
import {getHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {HelmAtom, HelmMol} from '@datagrok-libraries/bio/src/helm/types';
// import {FormsViewer} from '@datagrok-libraries/utils/src/viewers/forms-viewer';

import {RuleInputs, RULES_PATH, RULES_STORAGE_NAME} from './pt-rules';
import {addTransformedColumn} from './pt-conversion';

import {handleError} from './utils';
import {defaultErrorHandler} from '../utils/err-info';
import {getLibrariesList} from './utils';
import {getPtEnumeratorHelm, PT_HELM_EXAMPLE} from './pt-enumeration-helm';
import {getEnumerationChem, PT_CHEM_EXAMPLE} from './pt-enumeration-chem';
import {
  PolyToolEnumeratorParams, PolyToolEnumeratorType, PolyToolEnumeratorTypes, PolyToolPlaceholders
} from './types';

import {_package} from '../package';
import {PolyToolPlaceholdersInput} from './pt-placeholders-input';
import {InputBase} from 'datagrok-api/dg';
import {PT_ERROR_DATAFRAME, PT_UI_ADD_HELM, PT_UI_DIALOG_CONVERSION, PT_UI_DIALOG_ENUMERATION, PT_UI_GET_HELM, PT_UI_RULES_USED, PT_UI_USE_CHIRALITY, PT_WARNING_COLUMN} from './const';
import {PolyToolEnumerateDialog} from './pt-enumeration-helm-dialog';

export async function polyToolEnumerateHelmUI(cell?: DG.Cell): Promise<void> {
  const maxWidth = window.innerWidth;
  const maxHeight = window.innerHeight;

  try {
    const resizeInputs = () => {
      const contentHeight = $(dialog.root).find('div.d4-dialog-contents').get(0)!.clientHeight;

      const fitInputs: { [idx: number]: number } = {1: 1 /*, 3: 0.5*/};
      const fitInputsSumHeight = Object.values(fitInputs).reduce((sum, h) => sum + h, 0);

      const otherInputsHeight: number = dialog.inputs.filter((input, idx) => !(idx in fitInputs))
        .map((input) => input.root.offsetHeight).reduce((sum, h) => sum + h, 0);
      const remainFitHeight = contentHeight - otherInputsHeight - 38;
      dialog.inputs.forEach((input, idx) => {
        if (idx in fitInputs) {
          const inputFitHeight = remainFitHeight * fitInputs[idx] / fitInputsSumHeight;
          input.root.style.height = `${inputFitHeight}px`;
        }
      });
    };
    const [dialog, inputs] = await PolyToolEnumerateDialog.create2(cell, resizeInputs);

    let isFirstShow = true;
    ui.onSizeChanged(dialog.root).subscribe(() => {
      if (isFirstShow) {
        const dialogInputList = dialog.inputs;
        const dialogRootCash = $(dialog.root);
        const contentMaxHeight = maxHeight
          - dialogRootCash.find('div.d4-dialog-header').get(0)!.offsetHeight
          - dialogRootCash.find('div.d4-dialog-footer').get(0)!.offsetHeight;

        // dialog.inputs2.macromolecule.root.style.backgroundColor = '#CCFFCC';

        const dialogWidth = maxWidth * 0.7;
        const dialogHeight = maxHeight * 0.7;

        // Centered, but resizable dialog
        dialog.root.style.width = `${Math.min(maxWidth, dialogWidth)}px`;
        dialog.root.style.height = `${Math.min(maxHeight, dialogHeight)}px`;
        dialog.root.style.left = `${Math.floor((maxWidth - dialog.root.offsetWidth) / 2)}px`;
        dialog.root.style.top = `${Math.floor((maxHeight - dialog.root.offsetHeight) / 2)}px`;

        isFirstShow = false;
      }

      resizeInputs();
    });

    _package.logger.debug('PolyToolEnumerateHelmUI: dialog before show');
    const res = dialog.show({width: Math.max(350, maxWidth * 0.7), /* center: true,*/ resizable: true});
    _package.logger.debug('PolyToolEnumerateHelmUI: dialog after show');
    const k = 42;
  } catch (_err: any) {
    grok.shell.warning('To run PolyTool Enumeration, sketch the macromolecule and select monomers to vary');
  }
}

export function polyToolEnumerateChemUI(cell?: DG.Cell): void {
  getPolyToolEnumerationChemDialog(cell)
    .then((dialog) => {
      dialog.show({resizable: true});
    })
    .catch((_err: any) => {
      grok.shell.warning('To run PolyTool Enumeration, sketch the molecule and specify the R group to vary');
    });
}

export async function getPolyToolConversionDialog(targetCol?: DG.Column): Promise<DG.Dialog> {
  const targetColumns = grok.shell.t.columns.bySemTypeAll(DG.SEMTYPE.MACROMOLECULE);
  if (!targetColumns)
    throw new Error(PT_ERROR_DATAFRAME);

  const targetColumnInput = ui.input.column('Column', {
    table: grok.shell.t, value: targetColumns[0],
    filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MACROMOLECULE
  });

  targetColumnInput.value = targetCol ? targetCol : targetColumnInput.value;

  const generateHelmChoiceInput = ui.input.bool(PT_UI_GET_HELM, {value: true});
  ui.tooltip.bind(generateHelmChoiceInput.root, PT_UI_ADD_HELM);

  const chiralityEngineInput = ui.input.bool(PT_UI_USE_CHIRALITY, {value: false});
  const ruleInputs = new RuleInputs(RULES_PATH, RULES_STORAGE_NAME, '.json');
  const rulesHeader = ui.inlineText([PT_UI_RULES_USED]);
  ui.tooltip.bind(rulesHeader, 'Add or specify rules to use');
  const rulesForm = await ruleInputs.getForm();

  const div = ui.div([
    targetColumnInput,
    generateHelmChoiceInput,
    chiralityEngineInput,
    rulesHeader,
    rulesForm
  ]);

  const dialog = ui.dialog(PT_UI_DIALOG_CONVERSION)
    .add(div)
    .onOK(async () => {
      const pi = DG.TaskBarProgressIndicator.create('PolyTool converting');
      try {
        const sequencesCol = targetColumnInput.value;
        if (!sequencesCol) {
          grok.shell.warning(PT_WARNING_COLUMN);
          return;
        }

        const files = await ruleInputs.getActive();

        addTransformedColumn(sequencesCol!,
          generateHelmChoiceInput.value!,
          files,
          chiralityEngineInput.value!);
      } catch (err: any) {
        handleError(err);
      } finally {
        pi.close();
      }
    });

  return dialog;
}

async function getPolyToolEnumerationChemDialog(cell?: DG.Cell): Promise<DG.Dialog> {
  const [libList, helmHelper] = await Promise.all([
    getLibrariesList(), getHelmHelper()]);

  let molValue = PT_CHEM_EXAMPLE;//cell ? cell.value : PT_CHEM_EXAMPLE;
  const molInput = new DG.chem.Sketcher(DG.chem.SKETCHER_MODE.EXTERNAL);
  molInput.syncCurrentObject = false;
  // sketcher.setMolFile(col.tags[ALIGN_BY_SCAFFOLD_TAG]);
  molInput.onChanged.subscribe((_: any) => {
    molValue = molInput.getMolFile();
  });
  molInput.root.classList.add('ui-input-editor');
  molInput.root.style.marginTop = '3px';
  molInput.setMolFile(molValue);

  //const helmInput = helmHelper.createHelmInput('Macromolecule', {value: helmValue});
  const screenLibrary = ui.input.choice('Library to use', {value: null, items: libList});

  molInput.root.setAttribute('style', `min-width:250px!important;`);
  molInput.root.setAttribute('style', `max-width:250px!important;`);
  screenLibrary.input.setAttribute('style', `min-width:250px!important;`);

  const div = ui.div([
    molInput.root,
    screenLibrary.root
  ]);

  const cccSubs = grok.events.onCurrentCellChanged.subscribe(() => {
    const cell = grok.shell.tv.dataFrame.currentCell;

    if (cell.column.semType === DG.SEMTYPE.MOLECULE)
      molInput.setValue(cell.value);
  });

  // Displays the molecule from a current cell (monitors changes)
  const dialog = ui.dialog(PT_UI_DIALOG_ENUMERATION)
    .add(div)
    .onOK(async () => {
      try {
        const molString = molInput.getMolFile();

        if (molString === undefined || molString === '') {
          grok.shell.warning('PolyTool: no molecule was provided');
        } else if (!molString.includes('R#')) {
          grok.shell.warning('PolyTool: no R group was provided');
        } else {
          const molecules = await getEnumerationChem(molString, screenLibrary.value!);
          const molCol = DG.Column.fromStrings('Enumerated', molecules);
          const df = DG.DataFrame.fromColumns([molCol]);
          grok.shell.addTableView(df);
        }
      } catch (err: any) {
        defaultErrorHandler(err);
      } finally {
        cccSubs.unsubscribe();
      }
    }).onCancel(() => {
      cccSubs.unsubscribe();
    });

  return dialog;
}
