/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {RuleInputs, RULES_PATH, RULES_STORAGE_NAME} from './pt-rules';
import {addTransformedColumn} from './pt-conversion';

import {handleError} from './utils';
import {getLibrariesList, HelmInput} from './pt-enumeration';

const PT_ERROR_DATAFRAME = 'No dataframe with macromolecule columns open';
const PT_WARNING_COLUMN = 'No marcomolecule column chosen!';

const PT_UI_GET_HELM = 'Get HELM';
const PT_UI_ADD_HELM = 'Add HELM column';
const PT_UI_USE_CHIRALITY = 'Chirality engine';
const PT_UI_DIALOG_CONVERSION = 'Poly Tool Conversion';
const PT_UI_DIALOG_ENUMERATION = 'Poly Tool Enumeration';
const PT_UI_RULES_USED = 'Rules used';

export async function getPolyToolConversionDialog(): Promise<DG.Dialog> {
  const targetColumns = grok.shell.t.columns.bySemTypeAll(DG.SEMTYPE.MACROMOLECULE);
  if (!targetColumns)
    throw new Error(PT_ERROR_DATAFRAME);

  const targetColumnInput = ui.columnInput(
    'Column', grok.shell.t, targetColumns[0], null,
    {filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MACROMOLECULE}
  );

  const generateHelmChoiceInput = ui.boolInput(PT_UI_GET_HELM, true);
  ui.tooltip.bind(generateHelmChoiceInput.root, PT_UI_ADD_HELM);

  const chiralityEngineInput = ui.boolInput(PT_UI_USE_CHIRALITY, false);
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

export async function getPolyToolEnumerationDialog(): Promise<DG.Dialog> {
  const helmInput = await HelmInput.init();

  const libList = await getLibrariesList();
  const screenLibrary = ui.choiceInput('library to use', null, libList);

  const div = ui.div([
    helmInput.getDiv(),
    screenLibrary.root
  ]);

  const dialog = ui.dialog(PT_UI_DIALOG_ENUMERATION)
    .add(div)
    .onOK(() => {
      try {
        const helmString = helmInput.getHelmString();
        const helmSelections = helmInput.getHelmSelections();
        const lib = '';
        const 
      } catch (err: any) {

      } finally {

      }
    }).onCancel(() => {

    });

  return dialog;
}
