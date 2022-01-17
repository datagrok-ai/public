import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';


export async function findMCS(smiles: string, df: DG.DataFrame): Promise<string> {
  return await grok.functions.call('Chem:FindMCS', {smiles, df});
}

export async function findRGroups(smiles: string, df: DG.DataFrame, core: string, prefix: string): Promise<DG.DataFrame> {
  return await grok.functions.call('Chem:FindRGroups', {smiles, df, core, prefix});
}

export async function smilesTo3DCoordinates(smiles: string): Promise<string> {
  return await grok.functions.call('Chem:SmilesTo3DCoordinates', {smiles});
}

export async function getDescriptorsTree(): Promise<any> {
  return JSON.parse((await grok.functions.call('Chem:DescTree')).replaceAll('\\"', '\'').replaceAll('\\', ''));
}

export async function getDescriptorsPy(smiles: string, df1: DG.DataFrame, selected: string, df2: DG.DataFrame): Promise<any> {
  return await grok.functions.call('Chem:Desc', {smiles, df1, selected, df2});
}
