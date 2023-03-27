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
//input: string smiles { semType: Molecule }
//output: widget result
export function admetWidget(smiles: string): DG.Widget<any> {
  return getModelsSingle(smiles);
}

/*function addColorCoding(table: DG.DataFrame, column: DG.Column) {
  table.onCurrentColChanged.subscribe((_) => {
    let column = table.currentCol;
    table.onMetadataChanged.subscribe((_) => {
      if (column.colors.getType() === 'Conditional') {
        column.tags[DG.TAGS.COLOR_CODING_CONDITIONAL] = `{"0-0.5":"#f1b6b4","0.5-1":"#b4f1bc"}`;
      }
      if (column.colors.getType() === 'Linear') {
        column.tags[DG.TAGS.COLOR_CODING_LINEAR] = `["#f1b6b4", "#b4f1bc"]`;
      }
    })
  });
}*/

//top-menu: Chem | ADME/Tox | Calculations...
//name: ADME/Tox...
export async function admetoxCalculators() {
  let table = grok.shell.tv.dataFrame;
  let col = table.columns.bySemType(DG.SEMTYPE.MOLECULE);
  if (col) {
    addPredictions(col, table);
  }
}

//top-menu: Chem | ADME/Tox | Add Form...
//name: addForm
export async function testLayout() {
  const df = grok.shell.tv.dataFrame;
  const col = df.columns.bySemType(DG.SEMTYPE.MOLECULE);
  if (col)
    await addForm(col, df);
  const layout = await _package.files.readAsText('layout.json');
  const tableName = df.name;
  const modifiedLayout = layout.replaceAll("tableName", tableName).replaceAll("smilesColumn", col!.name);
  console.log(JSON.parse(modifiedLayout));
  //console.log(JSON.parse(modifiedLayout));
  let view = grok.shell.tv;
  view.loadLayout(DG.ViewLayout.fromJson(modifiedLayout));
  addColorCoding(df.columns.names());
  addTooltip();
}
