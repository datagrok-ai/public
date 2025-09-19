import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Signals
//tags: app
//output: view result
//meta.icon: images/signals-icon.png
//meta.browsePath: Chem
export async function signalsApp() : Promise<any> {
  return await PackageFunctions.signalsApp();
}

//input: dynamic treeNode 
export async function signalsAppTreeBrowser(treeNode: any) : Promise<void> {
  await PackageFunctions.signalsAppTreeBrowser(treeNode);
}

//name: Signals
//input: string data 
export async function saveToSignalsEln(data: string) : Promise<void> {
  await PackageFunctions.saveToSignalsEln(data);
}
