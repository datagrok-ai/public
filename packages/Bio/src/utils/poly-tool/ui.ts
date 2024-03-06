/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {HELM_POLYMER_TYPE} from '@datagrok-libraries/bio/src/utils/const';
import {MonomerLibManager} from '../monomer-lib/lib-manager';
import {ALL_MONOMERS, CYCLIZATION_TYPE, TRANSFORMATION_TYPE} from './const';
import {addTransformedColumn} from './transformation';
import * as rxjs from 'rxjs';
//import {MetaData} from './types';

const RULE_PATH = 'System:AppData/Bio/polytool-rules/';

export function getPolyToolDialog(): DG.Dialog {
  //const monomerLib = MonomerLibManager.instance.getBioLib();
  const targetColumns = grok.shell.t.columns.bySemTypeAll(DG.SEMTYPE.MACROMOLECULE);
  if (!targetColumns)
    throw new Error('No dataframe with macromolecule columns open');

  const targetColumnInput = ui.columnInput(
    'Column', grok.shell.t, targetColumns[0], null,
    {filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MACROMOLECULE}
  );

  const generateHelmChoiceInput = ui.boolInput('Get HELM', true);
  ui.tooltip.bind(generateHelmChoiceInput.root, 'Add HELM column');

  // let rulesTable: DG.DataFrame = DG.DataFrame.create();

  const ruleFileInput = ui.button('ADD RULES', () => {
    DG.Utils.openFile({
      accept: '.csv',
      open: async (selectedFile) => {
        const content = await selectedFile.text();
        // rulesTable = DG.DataFrame.fromCsv(content);
        await grok.dapi.files.writeAsText(RULE_PATH + `${selectedFile.name}`, content);
        //console.log(df.toCsv());
      },
    });
  });

  const div = ui.div([
    targetColumnInput,
    generateHelmChoiceInput,
    ruleFileInput
  ]);

  const dialog = ui.dialog('Poly Tool')
    .add(div)
    .onOK(async () => {
      const molCol = targetColumnInput.value;
      if (!molCol) {
        grok.shell.warning('No marcomolecule column chosen!');
        return;
      }
      addTransformedColumn(molCol!, generateHelmChoiceInput.value!);
    }
    );

  return dialog;
}
