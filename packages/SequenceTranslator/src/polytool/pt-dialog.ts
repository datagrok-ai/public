/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {RuleInputs, RULES_PATH, RULES_STORAGE_NAME} from './pt-rules';
import {addTransformedColumn} from './pt-transformation';

const PT_ERROR_DATAFRAME = 'No dataframe with macromolecule columns open';
const PT_WARNING_COLUMN = 'No marcomolecule column chosen!';

const PT_UI_GET_HELM = 'Get HELM';
const PT_UI_ADD_HELM = 'Add HELM column';
const PT_UI_USE_CHIRALITY = 'Chirality engine';
const PT_UI_DIALOG_NAME = 'Poly Tool';
const PT_UI_RULES_USED = 'Rules used';

export async function getPolyToolDialog(): Promise<DG.Dialog> {
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
  const rulesForm = await ruleInputs.getForm();

  const div = ui.div([
    targetColumnInput,
    generateHelmChoiceInput,
    chiralityEngineInput,
    PT_UI_RULES_USED,
    rulesForm
  ]);

  const dialog = ui.dialog(PT_UI_DIALOG_NAME)
    .add(div)
    .onOK(async () => {
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
    });

  return dialog;
}
