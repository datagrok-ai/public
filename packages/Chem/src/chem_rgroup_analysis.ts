/* Do not change these import lines to match external modules in webpack configuration */
import * as DG from 'datagrok-api/dg';
import { rgroupGetter } from './scripts-api';

export let _package = new DG.Package();


export async function getRGroups(smiles: DG.Column, core: string, prefix: string) {
  let result: DG.DataFrame = await rgroupGetter(smiles.name, smiles.dataFrame, core, prefix);

  // for (let col of result.columns) {
  //   for (let i = 0; i < col.length; i++) {
  //     col.set(i, convertToRDKit(col.get(i)));
  //   }
  //   col.compact();
  // }
  return result;
}

export function convertToRDKit(smiles: string | null): string | null {
  if (smiles !== null) {
    let regexConv: RegExp = /(\[)(R)(\d+)(\])/g;
    let match = regexConv.exec(smiles);
    if (match !== null) {
      smiles = smiles.replace(regexConv, `${match[1]}*:${match[3]}${match[4]}`);
    }
  }
  return smiles;
}
