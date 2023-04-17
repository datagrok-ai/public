import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { getRdKitModule } from './utils/chem-common-rdkit';
import { RDMol } from '@datagrok-libraries/chem-meta/src/rdkit-api';


export async function findMCS(molecules: string, df: DG.DataFrame): Promise<string> {
  // const returnSmarts = !!smarts;
  // return await grok.functions.call('Chem:FindMCS', {molecules, df, returnSmarts});

  let rdkit = getRdKitModule();
  //let test = rdkit.get_test_string();
  //console.log(test);

  let molCol = df.columns.byName(molecules);
  let n = molCol.length;

  let mols: RDMol[] = [];
  let arr = new Uint32Array(n);
  for(let i = 0; i < molCol.length; i++) {
    mols.push(rdkit.get_mol(molCol.get(i)));
    //@ts-ignore
    arr[i] = mols[i].$$.ptr;
  }
  // let mol1 = rdkit.get_mol("CCCCCC");
  // let mol2 = rdkit.get_mol("CCCC");

  //@ts-ignore
  let buff = rdkit.asm.Zb(n*4);

  // >> 2 is the reduction of element number of 32 bit vs 8 bit for offset
  //@ts-ignore
  rdkit.HEAPU32.set(arr, buff >> 2);

  let smarts: string = rdkit.get_mcs(buff, n);

  //@ts-ignore
  //rdkit.asm.gb("get_mcs", "void", ["number", "number"], [2, buff]);

  //@ts-ignore
  rdkit.asm.$b(buff);

  for(let i = 0; i < molCol.length; i++)
    mols[i].delete();

  return smarts;
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

export async function getDescriptorsPy(
  smiles: string,
  df1: DG.DataFrame,
  selected: string,
  df2: DG.DataFrame): Promise<any> {
  return await grok.functions.call('Chem:Desc', {smiles, df1, selected, df2});
}

export async function generateScaffoldTree(
  data: DG.DataFrame,
  smilesColumn: string,
  ringCutoff: number = 0,
  dischargeAndDeradicalize: boolean = false) : Promise<string> {
    return await grok.functions.call('Chem: GenerateScaffoldTree', {
      data, smilesColumn, ringCutoff, dischargeAndDeradicalize
    });
}
