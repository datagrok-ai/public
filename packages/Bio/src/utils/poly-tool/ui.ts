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

  let rulesTable: DG.DataFrame = DG.DataFrame.create();

  const ruleFileInput = ui.button('ADD RULES', () => {
    DG.Utils.openFile({
      accept: '.csv',
      open: async (selectedFile) => {
        const content = await selectedFile.text();
        rulesTable = DG.DataFrame.fromCsv(content);
        //console.log(df.toCsv());
      },
    });
  });
  // dialog.addButton(
  //   'Add',
  //   () => eventManager.addLibraryFile(),
  //   undefined,
  //   'Upload new HELM monomer library'
  // );


  // const ruleFileInput = DG.Utils.openFile({
  //   accept: '.csv',
  //   open: async (selectedFile) => {
  //     const content = await selectedFile.text();
  //     const name = selectedFile.name;
  //     const df = DG.DataFrame.fromCsv(content);

  //     console.log(df.toCsv());
  //   },
  // });

  //grok.data.files.openTable('Samples:Files/chem/smiles_10K_with_activities.csv')
  // const file = await loadFileAsText(tableName);
  // const df = DG.DataFrame.fromCsv(file);
  // df.name = tableName.replace('.csv', '');

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
      addTransformedColumn(molCol!, rulesTable!, generateHelmChoiceInput.value!);
    }
    );

  return dialog;
}
