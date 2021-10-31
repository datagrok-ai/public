/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Column, DataFrame, TableView} from 'datagrok-api/dg';

export let _package = new DG.Package();


export async function getRGroups(smiles: Column, core: string) {

  let result: DataFrame = await grok.functions.call(
    "Chem:RGroupGetter", {
      'smiles': smiles.name,
      'df1': smiles.dataFrame,
      'core': core
    });

  let regexConv: RegExp = /(\[)(R)(\d+)(\])/

  function convert(smiles: string | null): string | null {
    if (smiles != null) {
      let match = regexConv.exec(smiles);
      console.log(smiles);
      if (match != null) {
        smiles = smiles.replace(regexConv, `${match[1]}*:${match[3]}${match[4]}`);
        console.log(smiles);
      }
    }
    return smiles;
  }

  for (let col of result.columns) {
    for (let i = 0; i < col.length; i++) {
      col.set(i, convert(col.get(i)));
    }
  }
  return result;

}

export async function getMCS(smiles: Column) {
  let res: string;
  res = await grok.functions.call(
    "Chem:MCSGetter", {
      'smiles': smiles.name,
      'df1': smiles.dataFrame
    });
  return res;
}