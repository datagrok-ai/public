/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { getModelsSingle, addForm, addTooltip, addColorCoding, addPredictions } from './admet-analysis/admet-calculation';

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

//top-menu: Chem | ADME/Tox | Calculations...
//name: ADME/Tox...
export async function admetoxCalculators() {
  const table = grok.shell.tv.dataFrame;
  addPredictions(table);
}

//top-menu: Chem | ADME/Tox | Add Form
export async function addFormViewer() {
  const df = grok.shell.tv.dataFrame;
  const col = df.columns.bySemType(DG.SEMTYPE.MOLECULE);
  if (col)
    await addForm(col, df);
  const layout = await _package.files.readAsText('form-layout.json');
  const tableName = df.name;
  const modifiedLayout = layout.replaceAll("tableName", tableName).replaceAll("smilesColumn", col!.name);
  let view = grok.shell.tv;
  view.loadLayout(DG.ViewLayout.fromJson(modifiedLayout));
  addTooltip();
}
