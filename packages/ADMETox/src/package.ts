/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { getModelsSingle, addForm, addCalculations } from './admet-analysis/admet-calculation';
import { ColumnInputOptions } from '@datagrok-libraries/utils/src/type-declarations';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: Biology | ADME/Tox
//tags: panel, chem, widgets
//input: semantic_value smiles { semType: Molecule }
//output: widget result
export function admetWidget(smiles: DG.SemanticValue): DG.Widget<any> {
  return getModelsSingle(smiles);
}

//top-menu: Chem | ADME/Tox | Calculate...
//name: ADME/Tox...
export async function admetoxCalculators() {
  const table = grok.shell.tv.dataFrame;
  const col = table.columns.bySemType(DG.SEMTYPE.MOLECULE);
  addCalculations(col!, table);
}

//name: admeFormEditor
//tags: editor
//input: funccall call
export async function admeFormEditor(call: DG.FuncCall) {
  const molColumns = grok.shell.tv.dataFrame.columns.bySemTypeAll(DG.SEMTYPE.MOLECULE);
  if (molColumns.length === 1) {
    call.func.prepare({molecules: molColumns[0]}).call(true);
  } else {
    const colInput = ui.columnInput('Molecules', grok.shell.tv.dataFrame, molColumns[0], null, {filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MOLECULE} as ColumnInputOptions);
    ui.dialog({title: 'Full profile'})
    .add(colInput)
    .onOK(async () => {
      call.func.prepare({molecules: colInput.value}).call(true);
    })
    .show();
  }
}

//top-menu: Chem | ADME/Tox | Full Profile
//description: Calculates all properties and visualizes on the form
//name: Full Profile...
//input: column molecules {semType: Molecule}
//editor: Admetox:admeFormEditor
export async function addFormViewer(molecules: DG.Column) {
  const df = grok.shell.tv.dataFrame;
  await addForm(molecules, df);
  const layout = await _package.files.readAsText('form-layout.json');
  const tableName = df.name;
  const modifiedLayout = layout.replaceAll("tableName", tableName).replaceAll("smilesColumn", molecules.name);
  const view = grok.shell.tv;
  view.loadLayout(DG.ViewLayout.fromJson(modifiedLayout));
}
