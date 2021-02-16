/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export let _package = new DG.Package();
let packageName = 'Chemblbrowser'

//name: test
//input: string s
export function test(s) {
  grok.shell.info(_package.webRoot);
}
//name: Browser
//tags: app
export function Browser() {
  getAllChemblStructures();
}



//name: getAllChemblStructures
//output: dataframe df
export async function getAllChemblStructures() {
  let queryName = 'allChemblStructures'
  return await grok.data.query(`${packageName}:${queryName}`);
  }