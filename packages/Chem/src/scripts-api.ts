import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as api from './package-api';

export async function findMCS(molecules: string, df: DG.DataFrame,
  exactAtomSearch: boolean, exactBondSearch: boolean): Promise<string> {
  return await api.scripts.findMCS(molecules, df, exactAtomSearch, exactBondSearch);
}

export async function findRGroups(
  molecules: string,
  df: DG.DataFrame,
  core: string,
  prefix: string): Promise<DG.DataFrame> {
  return await api.scripts.findRGroups(molecules, df, core, prefix);
}

export async function findRGroupsWithCore(
  molecules: string,
  df: DG.DataFrame,
  core: string,
  onlyMatchAtRGroups: boolean): Promise<DG.DataFrame> {
  return await api.scripts.findRGroupsWithCore(molecules, df, core, onlyMatchAtRGroups);
}

export async function smilesTo3DCoordinates(molecule: string): Promise<string> {
  return await api.scripts.smilesTo3DCoordinates(molecule);
}

export async function getDescriptorsPy(
  smiles: string,
  df1: DG.DataFrame,
  selected: string,
  df2: DG.DataFrame): Promise<any> {
  return await api.scripts.desc(smiles, df1, selected, df2);
}

export async function generateScaffoldTree(
  data: DG.DataFrame,
  smilesColumn: string,
  ringCutoff: number = 0,
  dischargeAndDeradicalize: boolean = false) : Promise<DG.FileInfo> {
  return await api.scripts.generateScaffoldTree(data, smilesColumn, ringCutoff, dischargeAndDeradicalize);
}
