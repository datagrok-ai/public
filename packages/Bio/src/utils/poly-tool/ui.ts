/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {addTransformedColumn} from './transformation';
import {RULES_PATH, RULES_STORAGE_NAME} from './transformation';
import {ActiveFiles} from '@datagrok-libraries/utils/src/settings/active-files-base';

class RuleInputs extends ActiveFiles {
  constructor(path: string, userStorageName: string, ext: string ) {
    super(path, userStorageName, ext);
  }
}

export async function getPolyToolDialog(): Promise<DG.Dialog> {
  const targetColumns = grok.shell.t.columns.bySemTypeAll(DG.SEMTYPE.MACROMOLECULE);
  if (!targetColumns)
    throw new Error('No dataframe with macromolecule columns open');

  const targetColumnInput = ui.columnInput(
    'Column', grok.shell.t, targetColumns[0], null,
    {filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MACROMOLECULE}
  );

  const generateHelmChoiceInput = ui.boolInput('Get HELM', true);
  ui.tooltip.bind(generateHelmChoiceInput.root, 'Add HELM column');

  const chiralityEngineInput = ui.boolInput('Chirality engine', false);
  const ruleInputs = new RuleInputs(RULES_PATH, RULES_STORAGE_NAME, '.csv');
  const rulesForm = await ruleInputs.getForm();

  const div = ui.div([
    targetColumnInput,
    generateHelmChoiceInput,
    chiralityEngineInput,
    'Rules used',
    rulesForm
  ]);

  const dialog = ui.dialog('Poly Tool')
    .add(div)
    .onOK(async () => {
      const molCol = targetColumnInput.value;
      if (!molCol) {
        grok.shell.warning('No marcomolecule column chosen!');
        return;
      }

      const files = await ruleInputs.getActive();

      addTransformedColumn(molCol!,
        generateHelmChoiceInput.value!,
        files,
        chiralityEngineInput.value!);
    });

  return dialog;
}
