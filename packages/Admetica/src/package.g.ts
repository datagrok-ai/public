import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: info
export function info() : void {
  PackageFunctions.info();
}

//name: Biology | Admetica
//tags: panel, chem, widgets
//input: semantic_value smiles { semType: Molecule }
//output: widget result
export async function admeticaWidget(semValue: DG.SemanticValue) : Promise<any> {
  return await PackageFunctions.admeticaWidget(semValue);
}

//input: string property 
//output: list<string> result
export async function getModels(property: string) : Promise<string[]> {
  return await PackageFunctions.getModels(property);
}

//name: AdmeticaHT
//tags: HitTriageFunction
//input: dataframe table 
//input: column molecules { semType: Molecule }
//input: list<string> absorption { choices: Admetica:getModels('Absorption'); nullable: true }
//input: list<string> distribution { choices: Admetica:getModels('Distribution'); nullable: true }
//input: list<string> metabolism { choices: Admetica:getModels('Metabolism'); nullable: true }
//input: list<string> excretion { choices: Admetica:getModels('Excretion'); nullable: true }
export async function admeticaHT(table: DG.DataFrame, molecules: DG.Column, absorption: string[], distribution: string[], metabolism: string[], excretion: string[]) : Promise<void> {
  await PackageFunctions.admeticaHT(table, molecules, absorption, distribution, metabolism, excretion);
}

//name: AdmeticaEditor
//tags: editor
//input: funccall call 
export function admeticaEditor(call: DG.FuncCall) : void {
  PackageFunctions.admeticaEditor(call);
}

//name: AdmeticaMenu
//input: dataframe table { description: Input data table }
//input: categorical molecules { type: categorical; semType: Molecule }
//input: string template 
//input: list<string> models 
//input: bool addPiechart 
//input: bool addForm 
//top-menu: Chem | Admetica | Сalculate...
//editor: Admetica:AdmeticaEditor
export async function admeticaMenu(table: DG.DataFrame, molecules: DG.Column, template: string, models: string[], addPiechart: boolean, addForm: boolean) : Promise<void> {
  await PackageFunctions.admeticaMenu(table, molecules, template, models, addPiechart, addForm);
}

//input: string molecule { semType: Molecule }
//input: string prop { choices: ["Caco2","Solubility","Lipophilicity","PPBR","VDss"] }
//output: double result
export async function admeProperty(molecule: string, prop: string) : Promise<number> {
  return await PackageFunctions.admeProperty(molecule, prop);
}

//name: Admetica
//tags: app
//output: view result
//meta.icon: images/vlaaivis.png
//meta.browsePath: Chem
export async function admeticaApp() : Promise<any> {
  return await PackageFunctions.admeticaApp();
}

//name: Admetica Demo
//description: Evaluating ADMET properties
//meta.demoPath: Cheminformatics | Admetica
export async function admeticaDemo() : Promise<any> {
  return await PackageFunctions.admeticaDemo();
}
