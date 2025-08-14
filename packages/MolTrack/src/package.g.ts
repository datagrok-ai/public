import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: CDD Vault
//tags: app
//output: view result
//meta.icon: images/cdd-icon-small.png
//meta.browsePath: Chem
export async function molTrackApp() {
  return PackageFunctions.molTrackApp();
}

//name: cddVaultAppTreeBrowser
//input: dynamic treeNode 
//input: view browseView 
export async function cddVaultAppTreeBrowser(appNode: any, browseView: any) {
  return PackageFunctions.cddVaultAppTreeBrowser(appNode, browseView);
}
