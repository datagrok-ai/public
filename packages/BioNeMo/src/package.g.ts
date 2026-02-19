import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: MolMIMModel
//input: string algorithm = 'CMA-ES' 
//input: double num_molecules = 30 
//input: string property_name = 'QED' 
//input: bool minimize = false 
//input: double min_similarity = 0.3 
//input: double particles = 30 
//input: double iterations = 10 
//input: string smi = '[H][C@@]12Cc3c[nH]c4cccc(C1=C[C@H](NC(=O)N(CC)CC)CN2C)c34' { semType: Molecule }
export async function molMIMModel(algorithm: string, num_molecules: number, property_name: string, minimize: boolean, min_similarity: number, particles: number, iterations: number, smi: string) : Promise<void> {
  await PackageFunctions.molMIMModel(algorithm, num_molecules, property_name, minimize, min_similarity, particles, iterations, smi);
}

//name: EsmFold
//input: dataframe table 
//input: column sequences { semType: Macromolecule }
//output: dataframe result { action: join(table) }
//top-menu: Bio | Folding | EsmFold...
export async function esmFoldModel(table: DG.DataFrame, sequences: DG.Column) : Promise<any> {
  return await PackageFunctions.esmFoldModel(table, sequences);
}

//name: Bio | EsmFold
//input: semantic_value sequence { semType: Macromolecule }
//output: widget result
export async function esmFoldModelPanel(sequence: DG.SemanticValue) : Promise<any> {
  return await PackageFunctions.esmFoldModelPanel(sequence);
}

//output: list<string> result
export async function getTargetFiles() : Promise<string[]> {
  return await PackageFunctions.getTargetFiles();
}

//input: string ligand 
//input: string target 
//input: double poses 
//output: string result
//meta.cache: client
//meta.cache.invalidateOn: 0 * * * *
export async function diffDockModelScript(ligand: string, target: string, poses: number) : Promise<any> {
  return await PackageFunctions.diffDockModelScript(ligand, target, poses);
}

//name: DiffDock
//input: dataframe table 
//input: column ligands { semType: Molecule }
//input: string target { choices: Bionemo: getTargetFiles }
//input: double poses = 10 
//output: dataframe result { action: join(table) }
//top-menu: Chem | Docking | DiffDock...
export async function diffDockModel(table: DG.DataFrame, ligands: DG.Column, target: string, poses: number) : Promise<any> {
  return await PackageFunctions.diffDockModel(table, ligands, target, poses);
}

//name: Biology | DiffDock
//input: semantic_value smiles { semType: Molecule }
//output: widget result
//meta.role: widgets,panel
export async function diffDockPanel(smiles: DG.SemanticValue) : Promise<any> {
  return await PackageFunctions.diffDockPanel(smiles);
}
