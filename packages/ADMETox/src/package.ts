/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import { getModelsSingle, performChemicalPropertyPredictions, addSparklines } from './utils/admetox-utils';
import { properties } from './utils/admetox-utils';
import { AdmeticaBaseEditor } from './utils/admetox-editor';
import { _demoAdmetox } from './demo/demo-admetox';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: Biology | ADME/Tox
//tags: panel, chem, widgets
//input: semantic_value smiles { semType: Molecule }
//output: widget result
export async function admetWidget(semValue: DG.SemanticValue): Promise<DG.Widget<any>> {
  const smiles = await grok.functions.call('Chem:convertMolNotation',
    {molecule: semValue.value, sourceNotation: DG.chem.Notation.Unknown, targetNotation: DG.chem.Notation.Smiles});

  return await getModelsSingle(smiles, semValue);
}

//name: getModels
//input: string property
//output: list<string> result
export function getModels(property: string): string[] {
  return properties[property].models
    .filter((model: any) => !model.skip)
    .map((model: any) => model.name);;
}

//name: ADMETox
//tags: HitTriageFunction
//input: dataframe table
//input: column molecules {semType: Molecule}
//input: list<string> absorption {choices: ADMETox:getModels('Absorption'); nullable: true}
//input: list<string> distribution {choices: ADMETox:getModels('Distribution'); nullable: true}
//input: list<string> metabolism {choices: ADMETox:getModels('Metabolism'); nullable: true}
//input: list<string> excretion {choices: ADMETox:getModels('Excretion'); nullable: true}
export async function admetox(
  table: DG.DataFrame, molecules: DG.Column, absorption: string[], distribution: string[], metabolism: string[], excretion: string[], addProbabilities: boolean
  ): Promise<void> {
    const resultString: string = [...absorption, ...distribution, ...metabolism, ...excretion].join(',');
    await performChemicalPropertyPredictions(molecules, table, resultString);
}

//name: AdmeticaEditor
//tags: editor
//input: funccall call
export function AdmeticaEditor(call: DG.FuncCall): void {
  const funcEditor = new AdmeticaBaseEditor();
  ui.dialog({title: 'ADME/Tox'})
    .add(funcEditor.getEditor())
    .onOK(async () => {
      const params = funcEditor.getParams();
      call.func.prepare({
        table: params.table,
        molecules: params.col,
        templates: params.templatesName,
        models: params.models,
        addPiechart: params.addPiechart,
        addForm: params.addForm
      }).call(true);
    }).show();
}

//top-menu: Chem | ADME/Tox | Сalculate...
//name: Admetica
//input: dataframe table [Input data table]
//input: column molecules {type:categorical; semType: Molecule}
//input: string templates
//input: list<string> models
//input: bool addPiechart
//input: bool addForm
//editor: Admetox: AdmeticaEditor
export async function admetica(table: DG.DataFrame, molecules: DG.Column, templates: string, models: string[], addPiechart: boolean, addForm: boolean): Promise<void> {
  await performChemicalPropertyPredictions(molecules, table, models.join(','), templates, addPiechart, addForm);
}

//name: Demo Admetox
//meta.demoPath: Cheminformatics | ADMETox
export async function demoAdmetox(): Promise<void> {
  await _demoAdmetox();
}