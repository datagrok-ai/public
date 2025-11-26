import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';


export namespace Chem {
  export async function calculateLogD(table: DG.DataFrame, molecules: DG.Column, pH: number): Promise<DG.DataFrame> {
    return await grok.functions.call('Chem:CalculateLogD', {table, molecules, pH});
  }

  export async function calculateLogP(table: DG.DataFrame, molecules: DG.Column): Promise<DG.DataFrame> {
    return await grok.functions.call('Chem:CalculateLogP', {table, molecules});
  }

  export async function calculateLogS(table: DG.DataFrame, molecules: DG.Column): Promise<void> {
    return await grok.functions.call('Chem:CalculateLogS', {table, molecules});
  }

  export async function calculatePI(table: DG.DataFrame, molecules: DG.Column, pI_mean: boolean, pI_IPC2_peptide: boolean, pI_IPC_peptide: boolean, pI_ProMoST: boolean, pI_Gauci: boolean, pI_Grimsley: boolean, pI_Thurlkill: boolean, pI_Lehninger: boolean, pI_Toseland: boolean): Promise<DG.DataFrame> {
    return await grok.functions.call('Chem:CalculatePI', {table, molecules, pI_mean, pI_IPC2_peptide, pI_IPC_peptide, pI_ProMoST, pI_Gauci, pI_Grimsley, pI_Thurlkill, pI_Lehninger, pI_Toseland});
  }

  export async function calculatePKa(table: DG.DataFrame, molecules: DG.Column, pKa_acidic_list: boolean, pKa_basic_list: boolean, pKa_strongest_acidic: boolean, pKa_strongest_basic: boolean): Promise<DG.DataFrame> {
    return await grok.functions.call('Chem:CalculatePKa', {table, molecules, pKa_acidic_list, pKa_basic_list, pKa_strongest_acidic, pKa_strongest_basic});
  }

  export async function curate(data: DG.DataFrame, molecules: DG.Column, kekulization: boolean, normalization: boolean, reionization: boolean, neutralization: boolean, tautomerization: boolean, mainFragment: boolean): Promise<DG.DataFrame> {
    return await grok.functions.call('Chem:Curate', {data, molecules, kekulization, normalization, reionization, neutralization, tautomerization, mainFragment});
  }

  export async function desc(smiles: string, df1: DG.DataFrame, selected: string, df2: DG.DataFrame): Promise<DG.DataFrame> {
    return await grok.functions.call('Chem:Desc', {smiles, df1, selected, df2});
  }

  export async function findMCS(molecules: string, df: DG.DataFrame, exactAtomSearch: boolean, exactBondSearch: boolean): Promise<string> {
    return await grok.functions.call('Chem:FindMCS', {molecules, df, exactAtomSearch, exactBondSearch});
  }

  export async function findRGroupsWithCore(molecules: string, df: DG.DataFrame, core: string, onlyMatchAtRGroups: boolean): Promise<DG.DataFrame> {
    return await grok.functions.call('Chem:FindRGroupsWithCore', {molecules, df, core, onlyMatchAtRGroups});
  }

  export async function findRGroups(molecules: string, df: DG.DataFrame, core: string, prefix: string): Promise<DG.DataFrame> {
    return await grok.functions.call('Chem:FindRGroups', {molecules, df, core, prefix});
  }

  export async function chemistryGasteigerPartialCharges(mol: string, contours: number): Promise<any> {
    return await grok.functions.call('Chem:ChemistryGasteigerPartialCharges', {mol, contours});
  }

  export async function generateConformers(molecule: string, num_conformers: number, optimize: boolean, rms_threshold: number, max_attempts: number, random_seed: number): Promise<DG.DataFrame> {
    return await grok.functions.call('Chem:GenerateConformers', {molecule, num_conformers, optimize, rms_threshold, max_attempts, random_seed});
  }

  export async function inchiToMol(id: string): Promise<string> {
    return await grok.functions.call('Chem:InchiToMol', {id});
  }

  export async function mutate(molecule: string, steps: number, randomize: boolean, maxRandomResults: number): Promise<DG.DataFrame> {
    return await grok.functions.call('Chem:Mutate', {molecule, steps, randomize, maxRandomResults});
  }

  export async function generateScaffoldTree(data: DG.DataFrame, smilesColumn: string, ringCutoff: number, dischargeAndDeradicalize: boolean): Promise<any> {
    return await grok.functions.call('Chem:GenerateScaffoldTree', {data, smilesColumn, ringCutoff, dischargeAndDeradicalize});
  }

  export async function smilesTo3DCoordinates(molecule: string): Promise<string> {
    return await grok.functions.call('Chem:SmilesTo3DCoordinates', {molecule});
  }

  export async function amideReaction(amines: DG.DataFrame, amine_molecules: DG.Column, acids: DG.DataFrame, acid_molecules: DG.Column): Promise<DG.DataFrame> {
    return await grok.functions.call('Chem:AmideReaction', {amines, amine_molecules, acids, acid_molecules});
  }

  export async function butinaMoleculesClustering(data: DG.DataFrame, molecules: DG.Column): Promise<DG.DataFrame> {
    return await grok.functions.call('Chem:ButinaMoleculesClustering', {data, molecules});
  }

  export async function usrcat(data: DG.DataFrame, smiles: DG.Column): Promise<any> {
    return await grok.functions.call('Chem:USRCAT', {data, smiles});
  }

  export async function filterByCatalogs(data: DG.DataFrame, smiles: DG.Column, catalog: string): Promise<DG.DataFrame> {
    return await grok.functions.call('Chem:FilterByCatalogs', {data, smiles, catalog});
  }

  export async function murckoScaffolds(data: DG.DataFrame, smiles: DG.Column): Promise<DG.DataFrame> {
    return await grok.functions.call('Chem:MurckoScaffolds', {data, smiles});
  }

  export async function similarityMapsUsingFingerprints(mol: string, refmol: string, radius: number): Promise<number> {
    return await grok.functions.call('Chem:SimilarityMapsUsingFingerprints', {mol, refmol, radius});
  }

  export async function chemicalSpaceUsingTSNE(data: DG.DataFrame, smiles: DG.Column, components: number, minClusterSize: number): Promise<any> {
    return await grok.functions.call('Chem:ChemicalSpaceUsingTSNE', {data, smiles, components, minClusterSize});
  }

  export async function twoComponentReaction(data1: DG.DataFrame, reactants1: DG.Column, data2: DG.DataFrame, reactants2: DG.Column, reaction: string, matrixExpansion: boolean, randomize: boolean, seed: number, maxRandomReactions: number): Promise<DG.DataFrame> {
    return await grok.functions.call('Chem:TwoComponentReaction', {data1, reactants1, data2, reactants2, reaction, matrixExpansion, randomize, seed, maxRandomReactions});
  }

  export async function chemicalSpaceUsingUMAP(data: DG.DataFrame, smiles: DG.Column, neighbors: number, minClusterSize: number): Promise<any> {
    return await grok.functions.call('Chem:ChemicalSpaceUsingUMAP', {data, smiles, neighbors, minClusterSize});
  }

  export async function testPythonRunning(x: number, y: number): Promise<number> {
    return await grok.functions.call('Chem:TestPythonRunning', {x, y});
  }

  export async function getDescriptorsPy(
    smiles: string,
    df1: DG.DataFrame,
    selected: string,
    df2: DG.DataFrame): Promise<any> {
    return await grok.functions.call('Chem:Desc', {smiles, df1, selected, df2});
  }
}
