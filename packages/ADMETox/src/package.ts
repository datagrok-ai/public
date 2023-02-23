/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { addPredictions } from './admet-analysis/admet-calculation';
import { getModelsSingle, addForm } from './admet-analysis/admet-calculation';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: ADME/Tox
//tags: panel, chem, widgets
//input: string smiles { semType: Molecule }
//output: widget result
export function admetWidget(smiles: string): DG.Widget<any> {
  return getModelsSingle(smiles);
}

//top-menu: Chem | Analyze Structure | ADME/Tox Calculations...
//name: ADME/Tox...
export async function admetoxCalculators() {
  let table = grok.shell.tv.dataFrame;
  let col = table.columns.bySemType(DG.SEMTYPE.MOLECULE);
  if (col) {
    addPredictions(col, table).then(() => {
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
      })
    });
  }
}

//name: testDocker
export async function testDocker() {
  const dockerId = (await grok.dapi.dockerfiles.filter('admetox').first()).id;
  console.log(dockerId);
  const params: RequestInit = {
    method: 'POST',
    headers: {
      'Accept': 'text/csv',
      'Content-type': 'text/csv'
    },
    body: `smiles
    c1cc(O)ccc1`
  };
  
  const path = '/smiles/df_upload/?models=Ames';
  const response = await grok.dapi.dockerfiles.request(dockerId, path, params);
  console.log(response);
  if (response) {
    grok.shell.addTableView(DG.DataFrame.fromCsv(response));
  }
}

//top-menu: Chem | Analyze Structure | Add Form...
//name: testLayout
export async function testLayout() {
  const df = grok.shell.tv.dataFrame;
  const col = df.columns.bySemType(DG.SEMTYPE.MOLECULE);
  if (col)
    await addForm(col, df);
  const layout = await _package.files.readAsText('layout.json');
  const tableName = df.name;
  const modifiedLayout = layout.replaceAll("tableName", tableName);
  console.log(JSON.parse(modifiedLayout));
  //console.log(JSON.parse(modifiedLayout));
  let view = grok.shell.tv;
  view.loadLayout(DG.ViewLayout.fromJson(modifiedLayout));
}
