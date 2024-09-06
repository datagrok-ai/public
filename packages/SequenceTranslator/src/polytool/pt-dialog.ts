/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

import {getHelmHelper} from '@datagrok-libraries/bio/src/helm/helm-helper';

import {RuleInputs, RULES_PATH, RULES_STORAGE_NAME} from './pt-rules';
import {doPolyToolConvert} from './pt-conversion';

import {defaultErrorHandler} from '../utils/err-info';
import {getLibrariesList} from './utils';
import {getEnumerationChem, PT_CHEM_EXAMPLE} from './pt-enumeration-chem';
import {
  PT_ERROR_DATAFRAME, PT_UI_ADD_HELM, PT_UI_DIALOG_CONVERSION, PT_UI_DIALOG_ENUMERATION,
  PT_UI_GET_HELM, PT_UI_RULES_USED, PT_UI_USE_CHIRALITY, PT_WARNING_COLUMN
} from './const';

import {_package} from '../package';

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
      const ruleFileList = await ruleInputs.getActive();
      await polyToolConvertUI(targetColumnInput.value!, generateHelmChoiceInput.value!, chiralityEngineInput.value!, ruleFileList);
    });

  return dialog;
}

async function getPolyToolEnumerationChemDialog(cell?: DG.Cell): Promise<DG.Dialog> {
  const [libList, helmHelper] = await Promise.all([
    getLibrariesList(), getHelmHelper()]);

  const molStr = (cell && cell.rowIndex >= 0) ? cell.value : PT_CHEM_EXAMPLE;//cell ? cell.value : PT_CHEM_EXAMPLE;
  let molfileValue: string = await (async (): Promise<string> => {
    if (DG.chem.isMolBlock(molStr)) return molStr;
    return (await grok.functions.call('Chem:convertMolNotation', {
      molecule: molStr,
      sourceNotation: cell?.column.getTag(DG.TAGS.UNITS) ?? DG.chem.Notation.Unknown,
      targetNotation: DG.chem.Notation.MolBlock,
    }));
  })();

  const molInput = new DG.chem.Sketcher(DG.chem.SKETCHER_MODE.EXTERNAL);
  molInput.syncCurrentObject = false;
  // sketcher.setMolFile(col.tags[ALIGN_BY_SCAFFOLD_TAG]);
  molInput.onChanged.subscribe((_: any) => {
    molfileValue = molInput.getMolFile();
  });
  molInput.root.classList.add('ui-input-editor');
  molInput.root.style.marginTop = '3px';
  molInput.setMolFile(molfileValue);

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

export async function polyToolConvertUI(
  seqCol: DG.Column<string>, generateHelm: boolean, chiralityEngine: boolean, rules: string[]
): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('PolyTool converting');
  try {
    const table = seqCol.dataFrame;

    const [resHelmCol, resMolCol] = await doPolyToolConvert(seqCol,
      generateHelm,
      rules,
      chiralityEngine);
    resHelmCol.name = table.columns.getUnusedName(resHelmCol.name);
    resMolCol.name = table.columns.getUnusedName(resMolCol.name);
    table.columns.add(resHelmCol);
    table.columns.add(resMolCol);
  } finally {
    pi.close();
  }
}
