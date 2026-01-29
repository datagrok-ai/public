import {SdiscRuleViolationCellRenderer} from './utils/rule-violation-cell-renderer';
import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Clinical Case
//output: view result
//meta.role: app
//meta.icon: /img/clin_case_icon.png
export async function clinicalCaseApp() : Promise<any> {
  return await PackageFunctions.clinicalCaseApp();
}

//name: Preclinical Case
//output: view result
//meta.role: app
//meta.icon: /img/preclinical_case_icon.png
export async function PreclinicalCaseApp() : Promise<any> {
  return await PackageFunctions.PreclinicalCaseApp();
}

//input: dynamic treeNode 
export async function clinicalCaseAppTreeBrowser(treeNode: any) : Promise<void> {
  await PackageFunctions.clinicalCaseAppTreeBrowser(treeNode);
}

//input: dynamic treeNode 
export async function preclinicalCaseAppTreeBrowser(treeNode: any) : Promise<void> {
  await PackageFunctions.preclinicalCaseAppTreeBrowser(treeNode);
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
//input: bool ongoing { optional: true }
//input: string standard { optional: true; description: CDISC data format, either SDTM or SEND }
//output: widget result
export async function getListOfStudies(name?: string, description?: string, numSubjects?: number, numSubjectsOperator?: string, startDate?: any, startDateOperator?: string, endDate?: any, endDateOperator?: string, ongoing?: boolean, standard?: any) : Promise<any> {
  return await PackageFunctions.getListOfStudies(name, description, numSubjects, numSubjectsOperator, startDate, startDateOperator, endDate, endDateOperator, ongoing, standard);
}

//input: file folder 
//input: list<file> files 
//output: widget result
//meta.role: folderViewer
export async function clinicalCaseFolderLauncher(folder: DG.FileInfo, files: DG.FileInfo[]) : Promise<any> {
  return await PackageFunctions.clinicalCaseFolderLauncher(folder, files);
}

//input: list file 
//output: list<dataframe> result
//meta.role: fileHandler
//meta.ext: xpt
export async function xptFileHandler(file: DG.FileInfo) : Promise<any> {
  return await PackageFunctions.xptFileHandler(file);
}

//name: Run CDISC CORE Validation
//description: Run CDISC CORE validation on datasets
//input: string standard 
//input: string dataPath 
//input: string version 
//input: string outputFormat 
//input: dynamic options 
//output: string result
export async function runCoreValidate(standard: string, dataPath: string, version?: string, outputFormat?: string, options?: any) : Promise<string> {
  return await PackageFunctions.runCoreValidate(standard, dataPath, version, outputFormat, options);
}
//name: sdiscRuleViolationRenderer
//tags: cellRenderer
//output: grid_cell_renderer renderer
//meta.cellType: sdisc-rule-violation
export function _SdiscRuleViolationCellRenderer() {
  return new SdiscRuleViolationCellRenderer();
}

