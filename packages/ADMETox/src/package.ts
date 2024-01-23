/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { getModelsSingle, addForm, addCalculations, addColorCoding, performPredictions } from './admet-analysis/admet-calculation';
import { ColumnInputOptions } from '@datagrok-libraries/utils/src/type-declarations';
import MALFORMED_MOL_V200 from '@datagrok-libraries/chem-meta/src/formats/molfile-v2k-const';
import { properties } from './admet-analysis/const';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: Biology | ADME/Tox
//tags: panel, chem, widgets
//input: semantic_value smiles { semType: Molecule }
//output: widget result
export async function admetWidget(smiles: DG.SemanticValue): Promise<DG.Widget<any>> {
  const molStr = await grok.functions.call('Chem:convertMolNotation',
    {molecule: smiles.value, sourceNotation: 'unknown', targetNotation: 'molblock'});
  if (molStr === MALFORMED_MOL_V200)
    return new DG.Widget(ui.divText('Molecule is possibly malformed'));
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
  const initialNames = df.columns.names();
  await addForm(molecules, df);
  const layout = await _package.files.readAsText('form-layout.json');
  const tableName = df.name;
  const modifiedLayout = layout.replaceAll("tableName", tableName).replaceAll("smilesColumn", molecules.name);
  const view = grok.shell.tv;
  view.loadLayout(DG.ViewLayout.fromJson(modifiedLayout));
  const finalNames = df.columns.names();
  const difference = finalNames.filter(element => !initialNames.includes(element));
  addColorCoding(difference);
  view.grid.invalidate();
}

//name: getAbsorption
//input: string property
//output: list<string> result
export function getModels(property: string): string[] {
  return properties[property].models
    .filter((model: any) => !model.skip)
    .map((model: any) => model.name);;
}

//top-menu: Chem | Analyze | ADMETox...
//name: ADMETox
//tags: HitTriageFunction
//input: dataframe table
//input: column molecules {semType: Molecule}
//input: list<string> absorption {choices: ADMETox:getModels('Absorption'); nullable: true}
//input: list<string> distribution {choices: ADMETox:getModels('Distribution'); nullable: true}
//input: list<string> metabolism {choices: ADMETox:getModels('Metabolism'); nullable: true}
//input: list<string> excretion {choices: ADMETox:getModels('Excretion'); nullable: true}
//input: bool addProbabilities = false
export async function admetox(
  table: DG.DataFrame, molecules: DG.Column, absorption: string[], distribution: string[], metabolism: string[], excretion: string[], addProbabilities: boolean
  ): Promise<void> {
    const resultString: string = [...absorption, ...distribution, ...metabolism, ...excretion].join(',');
    await performPredictions(molecules, table, resultString, addProbabilities);
}