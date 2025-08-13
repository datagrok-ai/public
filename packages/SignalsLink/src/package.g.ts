import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Signals
//tags: app
//output: view result
//meta.icon: images/signals-icon.png
//meta.browsePath: Chem
export async function signalsApp() {
  return PackageFunctions.signalsApp();
}

//name: signalsAppTreeBrowser
//input: dynamic treeNode 
//input: view browseView 
export async function signalsAppTreeBrowser(treeNode: any, browseView: any) {
  return PackageFunctions.signalsAppTreeBrowser(treeNode, browseView);
}

//name: Signals
//input: string data 
export async function saveToSignalsEln(data: string) {
  return PackageFunctions.saveToSignalsEln(data);
}
