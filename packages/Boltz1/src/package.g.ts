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
//output: dataframe result
//top-menu: Bio | Folding | Boltz...
export async function folding(table: DG.DataFrame, sequences: DG.Column) : Promise<any> {
  return await PackageFunctions.folding(table, sequences);
}

//name: Boltz-1
//input: dataframe df { caption: Table }
//input: column molecules { semType: Molecule }
//input: string config { caption: Config; choices: Boltz1:getBoltzConfigFolders; description: 'Folder with config files for docking' }
//output: dataframe result
//top-menu: Chem | Docking | Boltz-1...
export async function docking(df: DG.DataFrame, molecules: DG.Column, config: string) : Promise<any> {
  return await PackageFunctions.docking(df, molecules, config);
}

//name: Boltz-1
//tags: panel, chem, widgets
//input: semantic_value molecule { semType: Molecule3D }
//output: widget result
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
//tags: app
//output: view result
//meta.browsePath: Bio
export async function boltz1App() : Promise<any> {
  return await PackageFunctions.boltz1App();
}
