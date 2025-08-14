import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: info
export function info() {
  return PackageFunctions.info();
}

//name: Biology | Admetica
//tags: panel, chem, widgets
//input: semantic_value smiles { semType: Molecule }
//output: widget result
export async function admeticaWidget(semValue: DG.SemanticValue) {
  return PackageFunctions.admeticaWidget(semValue);
}

//name: getModels
//input: string property 
//output: list<string> result
export async function getModels(property: string) {
  return PackageFunctions.getModels(property);
}

//name: AdmeticaHT
//tags: HitTriageFunction
//input: dataframe table 
//input: column molecules { semType: Molecule }
//input: list<string> absorption { choices: Admetica:getModels('Absorption'); nullable: true }
//input: list<string> distribution { choices: Admetica:getModels('Distribution'); nullable: true }
//input: list<string> metabolism { choices: Admetica:getModels('Metabolism'); nullable: true }
//input: list<string> excretion { choices: Admetica:getModels('Excretion'); nullable: true }
export async function admeticaHT(table: DG.DataFrame, molecules: DG.Column, absorption: string[], distribution: string[], metabolism: string[], excretion: string[]) {
  return PackageFunctions.admeticaHT(table, molecules, absorption, distribution, metabolism, excretion);
}

//name: AdmeticaEditor
//tags: editor
//input: funccall call 
export function admeticaEditor(call: DG.FuncCall) {
  return PackageFunctions.admeticaEditor(call);
}

//name: AdmeticaMenu
//input: dataframe table { description: Input data table }
//input: column molecules { type: categorical; semType: Molecule }
//input: string template 
//input: list<string> models 
//input: bool addPiechart 
//input: bool addForm 
//top-menu: Chem | Admetica | Сalculate...
//editor: Admetica: AdmeticaEditor
export async function admeticaMenu(table: DG.DataFrame, molecules: DG.Column, template: string, models: string[], addPiechart: boolean, addForm: boolean) {
  return PackageFunctions.admeticaMenu(table, molecules, template, models, addPiechart, addForm);
}

//name: admeProperty
//input: string molecule { semType: Molecule }
//input: string prop { choices: ['Caco2','Solubility','Lipophilicity','PPBR','VDss'] }
//output: double result
export async function admeProperty(molecule: string, prop: string) {
  return PackageFunctions.admeProperty(molecule, prop);
}

//name: Admetica
//tags: app
//output: view result
//meta.icon: images/vlaaivis.png
//meta.browsePath: Chem
export async function admeticaApp() {
  return PackageFunctions.admeticaApp();
}

//name: Admetica Demo
//description: Evaluating ADMET properties
//meta.demoPath: Cheminformatics | Admetica
export async function admeticaDemo() {
  return PackageFunctions.admeticaDemo();
}
