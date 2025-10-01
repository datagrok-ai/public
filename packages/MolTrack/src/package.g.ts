import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: CDD Vault
//tags: app
//output: view result
//meta.icon: images/cdd-icon-small.png
//meta.browsePath: Chem
export async function molTrackApp() : Promise<any> {
  return await PackageFunctions.molTrackApp();
}

//input: dynamic treeNode 
export async function cddVaultAppTreeBrowser(appNode: any) : Promise<void> {
  await PackageFunctions.cddVaultAppTreeBrowser(appNode);
}
