import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: info
//output: dynamic result
export function info() {
  return PackageFunctions.info();
}

//name: getAutoDockService
//output: dynamic result
export async function getAutoDockService() {
  return PackageFunctions.getAutoDockService();
}

//name: autoDockApp
export async function autoDockApp() {
  return PackageFunctions.autoDockApp();
}

//name: getConfigFiles
//output: list<string> configFiles
export async function getConfigFiles() {
  return PackageFunctions.getConfigFiles();
}

//name: dockLigandCached
//input: string jsonForm 
//input: string containerId 
//output: string result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 1 * *
export async function dockLigandCached(jsonForm: string, containerId: string) {
  return PackageFunctions.dockLigandCached(jsonForm, containerId);
}

//name: Autodock
//description: Autodock plugin UI
//tags: HitTriageFunction
//input: dataframe table 
//input: column ligands { type: categorical; semType: Molecule }
//input: string target { choices: Docking }
//input: int poses { default: 10 }
//top-menu: Chem | Docking | AutoDock...
export async function runAutodock5(table: DG.DataFrame, ligands: DG.Column, target: string, poses: number) {
  return PackageFunctions.runAutodock5(table, ligands, target, poses);
}

//name: isApplicableAutodock
//input: string molecule 
//output: bool result
export function isApplicableAutodock(molecule: string) {
  return PackageFunctions.isApplicableAutodock(molecule);
}

//name: AutoDock
//tags: panel, chem, widgets
//input: semantic_value molecule { semType: Molecule3D }
//output: widget result
//condition: Docking
export async function autodockWidget(molecule: DG.SemanticValue) {
  return PackageFunctions.autodockWidget(molecule);
}

//name: getAutodockSingle
//input: semantic_value molecule 
//input: bool showProperties 
//input: dataframe table 
//output: widget result
export async function getAutodockSingle(molecule: DG.SemanticValue, showProperties: boolean, table: DG.DataFrame) {
  return PackageFunctions.getAutodockSingle(molecule, showProperties, table);
}

//name: Demo Docking
//description: Small molecule docking to a macromolecule with pose visualization
//meta.demoPath: Bioinformatics | Docking
export async function demoDocking() {
  return PackageFunctions.demoDocking();
}

//name: Biology | AutoDock
//tags: panel
//input: semantic_value smiles { semType: Molecule }
//output: widget result
export async function autodockPanel(smiles: DG.SemanticValue) {
  return PackageFunctions.autodockPanel(smiles);
}

//name: Docking
//tags: app
//input: string path { meta.url: true; optional: true }
//output: view result
//meta.icon: images/docking-icon.png
//meta.browsePath: Bio
export async function dockingView(path: string) {
  return PackageFunctions.dockingView(path);
}
