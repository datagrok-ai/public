import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: info
export function info() : void {
  PackageFunctions.info();
}

//name: Biology | Admetica
//input: semantic_value smiles { semType: Molecule }
//output: widget result
//meta.role: widgets,panel
//meta.domain: chem
export async function admeticaWidget(semValue: DG.SemanticValue) : Promise<any> {
  return await PackageFunctions.admeticaWidget(semValue);
}

//input: string property 
//output: list<string> result
export async function getModels(property?: string) : Promise<string[]> {
  return await PackageFunctions.getModels(property);
}

//name: AdmeticaHT
//input: dataframe table 
//input: column molecules { semType: Molecule }
//input: list<string> absorption { choices: Admetica:getModels('Absorption'); nullable: true }
//input: list<string> distribution { choices: Admetica:getModels('Distribution'); nullable: true }
//input: list<string> metabolism { choices: Admetica:getModels('Metabolism'); nullable: true }
//input: list<string> excretion { choices: Admetica:getModels('Excretion'); nullable: true }
//meta.role: hitTriageFunction
export async function admeticaHT(table: DG.DataFrame, molecules: DG.Column, absorption: string[], distribution: string[], metabolism: string[], excretion: string[]) : Promise<void> {
  await PackageFunctions.admeticaHT(table, molecules, absorption, distribution, metabolism, excretion);
}

//name: AdmeticaEditor
//input: funccall call 
//meta.role: editor
export function admeticaEditor(call: DG.FuncCall) : void {
  PackageFunctions.admeticaEditor(call);
}

//name: AdmeticaMenu
//input: dataframe table { description: Input data table }
//input: column molecules { semType: Molecule }
//input: string template 
//input: list<string> models 
//input: bool addPiechart 
//input: bool addForm 
//top-menu: Chem | Admetica | Сalculate...
//editor: Admetica:AdmeticaEditor
export async function admeticaMenu(table: DG.DataFrame, molecules: DG.Column, template: string, models: string[], addPiechart: boolean, addForm: boolean) : Promise<void> {
  await PackageFunctions.admeticaMenu(table, molecules, template, models, addPiechart, addForm);
}

//input: column molecules { semType: Molecule }
//input: list<string> props { optional: true }
//output: dataframe result
//meta.vectorFunc: true
export async function getAdmeProperties(molecules: DG.Column, props?: string[]) : Promise<any> {
  return await PackageFunctions.getAdmeProperties(molecules, props);
}

//description: Predicts ADME properties for a given molecule.
//input: string molecule { semType: Molecule }
//output: dataframe result
export async function getAdmePropertiesSingle(molecule: string) : Promise<any> {
  return await PackageFunctions.getAdmePropertiesSingle(molecule);
}

//name: Admetica
//output: view result
//meta.icon: images/vlaaivis.png
//meta.browsePath: Chem
//meta.role: app
export async function runAdmeticaApplication() : Promise<any> {
  return await PackageFunctions.runAdmeticaApplication();
}

//name: Admetica Demo
//description: Evaluating ADMET properties
//meta.demoPath: Cheminformatics | Admetica
export async function admeticaDemo() : Promise<any> {
  return await PackageFunctions.admeticaDemo();
}
