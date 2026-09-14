import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: info
export function info() : void {
  PackageFunctions.info();
}

//description: Returns the AutoDock service client used by other packages
//output: object result
export async function getAutoDockService() : Promise<any> {
  return await PackageFunctions.getAutoDockService();
}

//description: Lists target folders that contain an AutoDock grid parameter (.gpf) file
//output: list<string> result
export async function getConfigFiles() : Promise<string[]> {
  return await PackageFunctions.getConfigFiles();
}

//description: Docks one ligand in the AutoDock container; results are cached per input
//input: string jsonForm 
//input: string containerId 
//output: string result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 1 * *
export async function dockLigandCached(jsonForm: string, containerId: string) : Promise<string> {
  return await PackageFunctions.dockLigandCached(jsonForm, containerId);
}

//name: AutoDock Results
//description: Dock small molecules to a macromolecule target with AutoDock and return binding poses and energies.
//input: dataframe table 
//input: column ligands { semType: Molecule; description: Small molecules to dock }
//input: string target { choices: Docking:getConfigFiles; description: Target folder with the macromolecule and docking config }
//input: double poses { description: Number of output conformations (poses) per molecule }
//output: dataframe result { action: join(table) }
//meta.vectorFunc: true
export async function getAutodockResults(table: DG.DataFrame, ligands: DG.Column, target: string, poses: number) : Promise<any> {
  return await PackageFunctions.getAutodockResults(table, ligands, target, poses);
}

//name: AutoDock
//description: Autodock plugin UI
//input: dataframe table { description: 'Input data table' }
//input: column ligands { semType: Molecule; description: 'Small molecules to dock' }
//input: string target { choices: Docking:getConfigFiles; description: 'Folder with config and macromolecule' }
//input: int poses = 10 { description: 'Number of output conformations for each small molecule' }
//meta.role: hitTriageFunction
//top-menu: Chem | Docking | AutoDock...
export async function runAutodock(table: DG.DataFrame, ligands: DG.Column, target: string, poses: number) : Promise<void> {
  await PackageFunctions.runAutodock(table, ligands, target, poses);
}

//description: Checks whether a Molecule3D value is an AutoDock pose with a binding energy
//input: string molecule 
//output: bool result
export function isApplicableAutodock(molecule: string) : boolean {
  return PackageFunctions.isApplicableAutodock(molecule);
}

//name: AutoDock
//description: Shows the receptor with the docked pose and its AutoDock energy terms
//input: semantic_value molecule { semType: Molecule3D }
//output: widget result
//meta.role: widgets,panel
//meta.domain: chem
//condition: Docking:isApplicableAutodock(molecule)
export async function autodockWidget(molecule: DG.SemanticValue) : Promise<any> {
  return await PackageFunctions.autodockWidget(molecule);
}

//description: Builds the AutoDock pose widget for a single Molecule3D value
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
//description: Docks the current molecule against a chosen target and shows the best pose
//input: semantic_value smiles { semType: Molecule }
//output: widget result
//meta.role: widgets,panel
export async function autodockPanel(smiles: DG.SemanticValue) : Promise<any> {
  return await PackageFunctions.autodockPanel(smiles);
}

//name: Docking
//description: Opens the Docking app
//input: string path { meta.url: true; optional: true }
//output: view result
//meta.icon: images/docking-icon.png
//meta.browsePath: Bio
//meta.role: app
export async function dockingView(path?: string) : Promise<any> {
  return await PackageFunctions.dockingView(path);
}
