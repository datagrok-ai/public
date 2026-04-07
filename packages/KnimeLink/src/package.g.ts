import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//description: KnimeLink function registration
//meta.role: autostart
export async function knimeLinkAutostart() : Promise<void> {
  await PackageFunctions.knimeLinkAutostart();
}

//name: KNIME
//output: view result
//meta.role: app
//meta.browsePath: Compute
export async function knimeLinkApp() : Promise<any> {
  return await PackageFunctions.knimeLinkApp();
}

//input: dynamic treeNode 
//meta.role: appTreeBrowser
//meta.app: KNIME
export async function knimeLinkAppTreeBrowser(treeNode: any) : Promise<void> {
  await PackageFunctions.knimeLinkAppTreeBrowser(treeNode);
}
