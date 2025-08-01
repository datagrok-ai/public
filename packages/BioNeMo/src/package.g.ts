import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: MolMIMModel
//input: string algorithm { default: CMA-ES }
//input: double num_molecules { default: 30 }
//input: string property_name { default: QED }
//input: bool minimize { default: false }
//input: double min_similarity { default: 0.3 }
//input: double particles { default: 30 }
//input: double iterations { default: 10 }
//input: string smi { default: [H][C@@]12Cc3c[nH]c4cccc(C1=C[C@H](NC(=O)N(CC)CC)CN2C)c34; semType: Molecule }
export async function molMIMModel(algorithm: string, num_molecules: number, property_name: string, minimize: boolean, min_similarity: number, particles: number, iterations: number, smi: string) {
  return PackageFunctions.molMIMModel(algorithm, num_molecules, property_name, minimize, min_similarity, particles, iterations, smi);
}

//name: EsmFoldModel
//input: dataframe df 
//input: column sequences { semType: Macromolecule }
//top-menu: Bio | Folding | EsmFold...
export async function esmFoldModel(df: DG.DataFrame, sequences: DG.Column) {
  return PackageFunctions.esmFoldModel(df, sequences);
}

//name: Bio | EsmFold
//input: semantic_value sequence { semType: Macromolecule }
//output: widget result
export async function esmFoldModelPanel(sequence: DG.SemanticValue) {
  return PackageFunctions.esmFoldModelPanel(sequence);
}

//name: getTargetFiles
//output: list result
export async function getTargetFiles() {
  return PackageFunctions.getTargetFiles();
}

//name: diffDockModelScript
//input: string ligand 
//input: string target 
//input: double poses 
//output: dynamic result
//meta.cache: client
//meta.cache.invalidateOn: 0 * * * *
export async function diffDockModelScript(ligand: string, target: string, poses: number) {
  return PackageFunctions.diffDockModelScript(ligand, target, poses);
}

//name: DiffDockModel
//input: dataframe df 
//input: column ligands { semType: Molecule }
//input: string target { choices: Bionemo: getTargetFiles }
//input: double poses { default: 5 }
//top-menu: Chem | Docking | DiffDock...
export async function diffDockModel(df: DG.DataFrame, ligands: DG.Column, target: string, poses: number) {
  return PackageFunctions.diffDockModel(df, ligands, target, poses);
}

//name: Biology | DiffDock
//tags: panel, widgets
//input: semantic_value smiles { semType: Molecule }
//output: widget result
export async function diffDockPanel(smiles: DG.SemanticValue) {
  return PackageFunctions.diffDockPanel(smiles);
}
