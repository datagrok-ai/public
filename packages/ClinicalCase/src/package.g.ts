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

//name: Get list of studies
//description: Return list of clinical and preclinical studies loaded into Clinical Case application
//input: string name { optional: true }
//input: string description { optional: true; description: More detailed study information including species, drug, dosing }
//input: int numSubjects { optional: true }
//input: string numSubjectsOperator { optional: true; description: >, <, = }
//input: datetime startDate { optional: true }
//input: string startDateOperator { optional: true; description: >, <, = }
//input: datetime endDate { optional: true }
//input: string endDateOperator { optional: true; description: >, <, = }
//input: string standard { optional: true; description: CDISC data format, either SDTM or SEND }
//output: widget result
export async function getListOfStudies(name?: string, description?: string, numSubjects?: number, numSubjectsOperator?: string, startDate?: any, startDateOperator?: string, endDate?: any, endDateOperator?: string, standard?: any) : Promise<any> {
  return await PackageFunctions.getListOfStudies(name, description, numSubjects, numSubjectsOperator, startDate, startDateOperator, endDate, endDateOperator, standard);
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
