import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: info
export function info() : void {
  PackageFunctions.info();
}

//output: list<string> result
export async function getBoltzConfigFolders() : Promise<string[]> {
  return await PackageFunctions.getBoltzConfigFolders();
}

//input: string config 
//input: string msa 
//output: string result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 1 * *
export async function runBoltz(config: string, msa: string) : Promise<string> {
  return await PackageFunctions.runBoltz(config, msa);
}

//name: Boltz
//input: dataframe table 
//input: column sequences { semType: Macromolecule }
//output: dataframe result { action: join(table) }
//top-menu: Bio | Folding | Boltz...
export async function folding(table: DG.DataFrame, sequences: DG.Column) : Promise<any> {
  return await PackageFunctions.folding(table, sequences);
}

//name: Boltz
//input: dataframe table 
//input: column ligands { semType: Molecule }
//input: string config { choices: Boltz1:getBoltzConfigFolders; description: 'Folder with config files for docking' }
//output: dataframe result { action: join(table) }
//top-menu: Chem | Docking | Boltz...
export async function docking(table: DG.DataFrame, ligands: DG.Column, config: string) : Promise<any> {
  return await PackageFunctions.docking(table, ligands, config);
}

//name: Boltz-1
//input: semantic_value molecule { semType: Molecule3D }
//output: widget result
//meta.role: widgets,panel
//meta.domain: chem
//condition: Boltz1:isApplicableBoltz(molecule)
export async function boltzWidget(molecule: DG.SemanticValue) : Promise<any> {
  return await PackageFunctions.boltzWidget(molecule);
}

//input: string molecule 
//output: bool result
export function isApplicableBoltz(molecule: string) : boolean {
  return PackageFunctions.isApplicableBoltz(molecule);
}

//name: Boltz-1
//output: view result
//meta.role: app
//meta.browsePath: Bio
export async function boltz1App() : Promise<any> {
  return await PackageFunctions.boltz1App();
}
