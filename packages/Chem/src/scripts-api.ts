import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';


export async function findMCS(molecules: string, df: DG.DataFrame, smarts?: boolean): Promise<string> {
  const returnSmarts = !!smarts;
  return await grok.functions.call('Chem:FindMCS', {molecules, df, returnSmarts});
}

export async function findRGroups(
  molecules: string,
  df: DG.DataFrame,
  core: string,
  prefix: string): Promise<DG.DataFrame> {
  return await grok.functions.call('Chem:FindRGroups', {molecules, df, core, prefix});
}

export async function smilesTo3DCoordinates(molecule: string): Promise<string> {
  return await grok.functions.call('Chem:SmilesTo3DCoordinates', {molecule});
}

export async function getDescriptorsTree(): Promise<{[key: string]: any}> {
  return JSON.parse((await grok.functions.call('Chem:DescTree')).replaceAll('\\"', '\'').replaceAll('\\', ''));
}

export async function getDescriptorsPy(
  smiles: string,
  df1: DG.DataFrame,
  selected: string,
  df2: DG.DataFrame): Promise<any> {
  return await grok.functions.call('Chem:Desc', {smiles, df1, selected, df2});
}

export async function generateScaffoldTree(
  data: DG.DataFrame, 
  smiles: string) : Promise<string> {
    return await grok.functions.call('Chem: GenerateScaffoldTree', {data, smiles});
}

export async function setupScaffold() : Promise<void> {
  return await grok.functions.call('Chem:SetupScaffold');
}
