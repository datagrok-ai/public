import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: info
export function info() : void {
  PackageFunctions.info();
}

//name: Biology | Admetica
//description: Panel with ADMET predictions for a molecule.
//input: semantic_value smiles { semType: Molecule; description: Molecule to predict. }
//output: widget result
//meta.role: widgets,panel
//meta.domain: chem
export async function admeticaWidget(semValue: DG.SemanticValue) : Promise<any> {
  return await PackageFunctions.admeticaWidget(semValue);
}

//description: Lists available ADMET properties.
//input: string property { description: Category: Absorption, Distribution, Metabolism, Excretion, Toxicity. }
//output: list<string> result
export async function getModels(property?: string) : Promise<string[]> {
  return await PackageFunctions.getModels(property);
}

//name: AdmeticaHT
//description: Runs ADMET predictions for Hit Triage.
//input: dataframe table { description: Table with molecules. }
//input: column molecules { semType: Molecule; description: Molecule column. }
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
//description: Predicts ADMET properties and appends result columns.
//input: dataframe table { description: Table with molecules. }
//input: column molecules { semType: Molecule; description: Molecule column. }
//input: string template { description: Optional JSON config. }
//input: list<string> models { description: Properties to compute. }
//input: bool addPiechart { description: Add a pie-chart column. }
//input: bool addForm { description: Add a form viewer. }
//top-menu: Chem | Admetica | Сalculate...
//editor: Admetica:AdmeticaEditor
export async function admeticaMenu(table: DG.DataFrame, molecules: DG.Column, template: string, models: string[], addPiechart: boolean, addForm: boolean) : Promise<void> {
  await PackageFunctions.admeticaMenu(table, molecules, template, models, addPiechart, addForm);
}

//description: Predicts ADMET properties for a molecule column.
//input: dataframe table { description: Target table for results. }
//input: column molecules { semType: Molecule; description: Molecule column. }
//input: list<string> props { optional: true; description: Properties to compute. All if omitted. }
//output: dataframe result { action: join(table) }
//meta.vectorFunc: true
export async function getAdmeProperties(table: DG.DataFrame, molecules: DG.Column, props?: string[]) : Promise<any> {
  return await PackageFunctions.getAdmeProperties(table, molecules, props);
}

//description: Predicts ADMET properties for a given molecule.
//input: string molecule { semType: Molecule; description: Molecule (SMILES or molfile). }
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
