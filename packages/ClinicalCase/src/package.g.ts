import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Clinical Case
//tags: app
//output: view result
//meta.browsePath: Clinical
export async function clinicalCaseApp() : Promise<any> {
  return await PackageFunctions.clinicalCaseApp();
}

//input: dynamic treeNode 
export async function clinicalCaseAppTreeBrowser(treeNode: any) : Promise<void> {
  await PackageFunctions.clinicalCaseAppTreeBrowser(treeNode);
}

//tags: folderViewer
//input: file folder 
//input: list<file> files 
//output: widget result
export async function clinicalCaseFolderLauncher(folder: DG.FileInfo, files: DG.FileInfo[]) : Promise<any> {
  return await PackageFunctions.clinicalCaseFolderLauncher(folder, files);
}

//tags: file-handler
//input: list file 
//output: list<dataframe> result
//meta.ext: xpt
export async function xptFileHandler(file: DG.FileInfo) : Promise<any> {
  return await PackageFunctions.xptFileHandler(file);
}
