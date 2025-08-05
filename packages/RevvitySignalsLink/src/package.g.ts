import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Revvity Signals
//tags: app
//output: view result
//meta.browsePath: Chem
export async function revvitySignalsLinkApp() {
  return PackageFunctions.revvitySignalsLinkApp();
}

//name: revvitySignalsLinkAppTreeBrowser
//input: dynamic treeNode 
//input: view browseView 
//output: dynamic result
export async function revvitySignalsLinkAppTreeBrowser(treeNode: any, browseView: DG.View) {
  return PackageFunctions.revvitySignalsLinkAppTreeBrowser(treeNode, browseView);
}

//name: Search Entities
//input: string query 
//input: string params 
//output: dataframe result
export async function searchEntities(query: string, params: string) {
  return PackageFunctions.searchEntities(query, params);
}

//name: Search Entities With Structures
//input: string query 
//input: string params 
//output: dataframe result
export async function searchEntitiesWithStructures(query: string, params: string) {
  return PackageFunctions.searchEntitiesWithStructures(query, params);
}

//name: Get Users
//output: string result
export async function getUsers() {
  return PackageFunctions.getUsers();
}

//name: Revvity Signals
//tags: panel, widgets
//input: string id { semType: RevvitySignalsId }
//output: widget result
export async function entityTreeWidget(id: string) {
  return PackageFunctions.entityTreeWidget(id);
}
