import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';


export async function findMCS(molecules: string, df: DG.DataFrame,
  exactAtomSearch: boolean, exactBondSearch: boolean): Promise<string> {
  return await grok.functions.call('Chem:FindMCS', {molecules, df, exactAtomSearch, exactBondSearch});
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

export async function generateScaffoldTree(
  data: DG.DataFrame,
  smilesColumn: string,
  ringCutoff: number = 0,
  dischargeAndDeradicalize: boolean = false) : Promise<string> {
  return await grok.functions.call('Chem: GenerateScaffoldTree', {
    data, smilesColumn, ringCutoff, dischargeAndDeradicalize,
  });
}
