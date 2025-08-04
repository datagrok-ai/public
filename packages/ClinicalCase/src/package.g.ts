import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Clinical Case
//tags: app
//output: view result
//meta.browsePath: Clinical
export async function clinicalCaseApp() {
  return PackageFunctions.clinicalCaseApp();
}

//name: clinicalCaseAppTreeBrowser
//input: dynamic treeNode 
//input: view browseView 
//output: dynamic result
export async function clinicalCaseAppTreeBrowser(treeNode: any, browseView: DG.View) {
  return PackageFunctions.clinicalCaseAppTreeBrowser(treeNode, browseView);
}

//name: clinicalCaseFolderLauncher
//tags: folderViewer
//input: file folder 
//input: list<file> files 
//output: dynamic result
export async function clinicalCaseFolderLauncher(folder: DG.FileInfo, files: DG.FileInfo[]) {
  return PackageFunctions.clinicalCaseFolderLauncher(folder, files);
}
