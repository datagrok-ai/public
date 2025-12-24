import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: info
export function info() : void {
  PackageFunctions.info();
}

//output: object result
export async function getAutoDockService() : Promise<any> {
  return await PackageFunctions.getAutoDockService();
}

//name: autoDockApp
export async function autoDockApp() : Promise<void> {
  await PackageFunctions.autoDockApp();
}

//output: list<string> result
export async function getConfigFiles() : Promise<string[]> {
  return await PackageFunctions.getConfigFiles();
}

//input: string jsonForm 
//input: string containerId 
//output: string result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 1 * *
export async function dockLigandCached(jsonForm: string, containerId: string) : Promise<string> {
  return await PackageFunctions.dockLigandCached(jsonForm, containerId);
}

//name: AutoDock
//description: Autodock plugin UI
//input: dataframe table { description: 'Input data table' }
//input: categorical ligands { type: categorical; semType: Molecule; description: 'Small molecules to dock' }
//input: string target { choices: Docking:getConfigFiles; description: 'Folder with config and macromolecule' }
//input: int poses = 10 { description: 'Number of output conformations for each small molecule' }
//meta.role: hitTriageFunction
//top-menu: Chem | Docking | AutoDock...
export async function runAutodock(table: DG.DataFrame, ligands: DG.Column, target: string, poses: number) : Promise<void> {
  await PackageFunctions.runAutodock(table, ligands, target, poses);
}

//input: string molecule 
//output: bool result
export function isApplicableAutodock(molecule: string) : boolean {
  return PackageFunctions.isApplicableAutodock(molecule);
}

//name: AutoDock
//input: semantic_value molecule { semType: Molecule3D }
//output: widget result
//meta.role: widgets,panel
//meta.domain: chem
//condition: Docking:isApplicableAutodock(molecule)
export async function autodockWidget(molecule: DG.SemanticValue) : Promise<any> {
  return await PackageFunctions.autodockWidget(molecule);
}

//input: semantic_value molecule 
//input: bool showProperties 
//input: dataframe table 
//output: widget result
export async function getAutodockSingle(molecule: DG.SemanticValue, showProperties: boolean, table?: DG.DataFrame) : Promise<any> {
  return await PackageFunctions.getAutodockSingle(molecule, showProperties, table);
}

//name: Demo Docking
//description: Small molecule docking to a macromolecule with pose visualization
//meta.demoPath: Bioinformatics | Docking
export async function demoDocking() : Promise<void> {
  await PackageFunctions.demoDocking();
}

//name: Biology | AutoDock
//input: semantic_value smiles { semType: Molecule }
//output: widget result
//meta.role: widgets,panel
export async function autodockPanel(smiles: DG.SemanticValue) : Promise<any> {
  return await PackageFunctions.autodockPanel(smiles);
}

//name: Docking
//input: string path { meta.url: true; optional: true }
//output: view result
//meta.icon: images/docking-icon.png
//meta.browsePath: Bio
//meta.role: app
export async function dockingView(path?: string) : Promise<any> {
  return await PackageFunctions.dockingView(path);
}
