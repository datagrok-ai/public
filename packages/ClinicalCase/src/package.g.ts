import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Clinical Case
//tags: app
//output: view result
//meta.browsePath: Clinical
export async function clinicalCaseApp() : Promise<any> {
  return PackageFunctions.clinicalCaseApp();
}

//name: clinicalCaseAppTreeBrowser
//input: dynamic treeNode 
//meta.role: appTreeBrowser
export async function clinicalCaseAppTreeBrowser(treeNode: any) : Promise<void> {
  PackageFunctions.clinicalCaseAppTreeBrowser(treeNode);
}

//name: clinicalCaseFolderLauncher
//tags: folderViewer
//input: file folder 
//input: list<file> files 
//output: widget result
export async function clinicalCaseFolderLauncher(folder: DG.FileInfo, files: DG.FileInfo[]) : Promise<any> {
  return PackageFunctions.clinicalCaseFolderLauncher(folder, files);
}
