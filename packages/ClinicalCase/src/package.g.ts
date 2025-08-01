import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Clinical Case
//tags: app
//output: view v
//meta.browsePath: Clinical
export async function clinicalCaseApp() {
  return PackageFunctions.clinicalCaseApp();
}

//name: clinicalCaseFolderLauncher
//tags: folderViewer
//input: file folder 
//input: list<file> files 
//output: widget res
export async function clinicalCaseFolderLauncher(folder: DG.FileInfo, files: DG.FileInfo[]) {
  return PackageFunctions.clinicalCaseFolderLauncher(folder, files);
}
