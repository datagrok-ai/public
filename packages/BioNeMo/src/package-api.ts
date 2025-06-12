import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';



export namespace Scripts {
  //Predicts the 3D structure of how a molecule interacts with a protein
  export async function diffdock(protein: string, ligand: string, num_poses: number, api_key: string): Promise<any> {
    return await grok.functions.call('Bionemo:Diffdock', { protein, ligand, num_poses, api_key });
  }

  //Predicts the 3D structure of a protein from its amino acid sequence
  export async function esmfold(sequence: string, api_key: string): Promise<string> {
    return await grok.functions.call('Bionemo:Esmfold', { sequence, api_key });
  }

  //MolMIM performs controlled generation, finding molecules with the right properties
  export async function molMIMGenerate(algorithm: string, num_molecules: number, property_name: string, minimize: boolean, min_similarity: number, particles: number, iterations: number, smi: string, api_key: string): Promise<string> {
    return await grok.functions.call('Bionemo:MolMIMGenerate', { algorithm, num_molecules, property_name, minimize, min_similarity, particles, iterations, smi, api_key });
  }
}

export namespace Funcs {
  export async function info(): Promise<any> {
    return await grok.functions.call('Bionemo:Info', {});
  }

  export async function molMIMModel(algorithm: string, num_molecules: number, property_name: string, minimize: boolean, min_similarity: number, particles: number, iterations: number, smi: string): Promise<any> {
    return await grok.functions.call('Bionemo:MolMIMModel', { algorithm, num_molecules, property_name, minimize, min_similarity, particles, iterations, smi });
  }

  export async function esmFoldModel(df: DG.DataFrame, sequences: DG.Column): Promise<any> {
    return await grok.functions.call('Bionemo:EsmFoldModel', { df, sequences });
  }

  export async function esmFoldModelPanel(sequence: any): Promise<any> {
    return await grok.functions.call('Bionemo:EsmFoldModelPanel', { sequence });
  }

  export async function getTargetFiles(): Promise<any> {
    return await grok.functions.call('Bionemo:GetTargetFiles', {});
  }

  export async function diffDockModelScript(ligand: string, target: string, poses: number): Promise<any> {
    return await grok.functions.call('Bionemo:DiffDockModelScript', { ligand, target, poses });
  }

  export async function diffDockModel(df: DG.DataFrame, ligands: DG.Column, target: string, poses: number): Promise<any> {
    return await grok.functions.call('Bionemo:DiffDockModel', { df, ligands, target, poses });
  }

  export async function diffDockPanel(smiles: any): Promise<any> {
    return await grok.functions.call('Bionemo:DiffDockPanel', { smiles });
  }
}
