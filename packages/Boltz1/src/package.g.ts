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
//meta.vectorFunc: true
//top-menu: Bio | Folding | Boltz...
export async function folding(table: DG.DataFrame, sequences: DG.Column) : Promise<any> {
  return await PackageFunctions.folding(table, sequences);
}

//name: Boltz
//input: dataframe table 
//input: column ligands { semType: Molecule }
//input: string config { choices: Boltz1:getBoltzConfigFolders; description: 'Folder with config files for docking' }
//output: dataframe result { action: join(table) }
//meta.vectorFunc: true
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

//output: list<string> result
export async function getBoltzStructureBindingConfigs() : Promise<string[]> {
  return await PackageFunctions.getBoltzStructureBindingConfigs();
}

//output: list<string> result
export async function getBoltzSmDesignConfigs() : Promise<string[]> {
  return await PackageFunctions.getBoltzSmDesignConfigs();
}

//output: list<string> result
export async function getBoltzSmScreenConfigs() : Promise<string[]> {
  return await PackageFunctions.getBoltzSmScreenConfigs();
}

//output: list<string> result
export async function getBoltzProteinDesignConfigs() : Promise<string[]> {
  return await PackageFunctions.getBoltzProteinDesignConfigs();
}

//output: list<string> result
export async function getBoltzProteinScreenConfigs() : Promise<string[]> {
  return await PackageFunctions.getBoltzProteinScreenConfigs();
}

//input: dataframe table 
//input: column ligands { semType: Molecule }
//input: string config { choices: Boltz1:getBoltzStructureBindingConfigs }
//output: dataframe result { action: join(table) }
//meta.vectorFunc: true
export async function boltzStructureAndBinding(table: DG.DataFrame, ligands: DG.Column, config: string) : Promise<any> {
  return await PackageFunctions.boltzStructureAndBinding(table, ligands, config);
}

//input: dataframe table 
//input: column molecules { semType: Molecule }
//output: dataframe result { action: join(table) }
//meta.vectorFunc: true
export async function boltzAdme(table: DG.DataFrame, molecules: DG.Column) : Promise<any> {
  return await PackageFunctions.boltzAdme(table, molecules);
}

//input: string config { choices: Boltz1:getBoltzSmDesignConfigs }
//input: double numMolecules 
//output: dataframe result
export async function boltzDesignSmallMolecules(config: string, numMolecules: number) : Promise<any> {
  return await PackageFunctions.boltzDesignSmallMolecules(config, numMolecules);
}

//input: dataframe table 
//input: column molecules { semType: Molecule }
//input: string config { choices: Boltz1:getBoltzSmScreenConfigs }
//output: dataframe result { action: join(table) }
//meta.vectorFunc: true
export async function boltzScreenSmallMolecules(table: DG.DataFrame, molecules: DG.Column, config: string) : Promise<any> {
  return await PackageFunctions.boltzScreenSmallMolecules(table, molecules, config);
}

//input: string config { choices: Boltz1:getBoltzProteinDesignConfigs }
//input: double numProteins 
//output: dataframe result
export async function boltzDesignProteins(config: string, numProteins: number) : Promise<any> {
  return await PackageFunctions.boltzDesignProteins(config, numProteins);
}

//input: dataframe table 
//input: column proteins { semType: Macromolecule }
//input: string config { choices: Boltz1:getBoltzProteinScreenConfigs }
//output: dataframe result { action: join(table) }
//meta.vectorFunc: true
export async function boltzScreenProteins(table: DG.DataFrame, proteins: DG.Column, config: string) : Promise<any> {
  return await PackageFunctions.boltzScreenProteins(table, proteins, config);
}
