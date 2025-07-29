import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: info
//description: Show package web root
export function info() {
  return PackageFunctions.info();
}

//name: getBoltzConfigFolders
//output: list result
export async function getBoltzConfigFolders() {
  return PackageFunctions.getBoltzConfigFolders();
}

//name: runBoltz
//input: string config 
//input: string msa 
//output: string result
//meta.cache: all
//meta.cache.invalidateOn: 0 0 1 * *
export async function runBoltz(config: string, msa: string) {
  return PackageFunctions.runBoltz(config, msa);
}

//name: Folding
//input: dataframe df 
//input: column sequences { semType: Macromolecule }
//output: dataframe result
//top-menu: Bio | Folding | Boltz-1...
export async function folding(df: DG.DataFrame, sequences: DG.Column) {
  return PackageFunctions.folding(df, sequences);
}

//name: Docking
//input: dataframe df 
//input: column molecules { semType: Molecule }
//input: string config { choices: Boltz1: getBoltzConfigFolders; caption: Folder with config files for docking }
//output: dataframe result
//top-menu: Chem | Docking | Boltz-1...
export async function docking(df: DG.DataFrame, molecules: DG.Column, config: string) {
  return PackageFunctions.docking(df, molecules, config);
}

//name: Boltz-1
//tags: panel, chem, widgets
//input: semantic_value molecule { semType: Molecule3D }
//output: widget result
//condition: Boltz1:isApplicableBoltz(molecule)
export async function boltzWidget(molecule: DG.SemanticValue) {
  return PackageFunctions.boltzWidget(molecule);
}

//name: isApplicableBoltz
//input: string molecule 
//output: bool result
export function isApplicableBoltz(molecule: string) {
  return PackageFunctions.isApplicableBoltz(molecule);
}

//name: Boltz-1
//tags: app
//output: view result
//meta.browsePath: Bio
export async function boltz1App() {
  return PackageFunctions.boltz1App();
}
