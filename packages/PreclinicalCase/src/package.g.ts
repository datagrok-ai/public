import {YAxisCellRenderer} from './utils/y-axis-cell-renderer';
import {SdiscRuleViolationCellRenderer} from './utils/rule-violation-cell-renderer';
import {CombinedMeasurementsCellRenderer} from './utils/combined-measurements-cell-renderer';
import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Preclinical Case
//output: view result
//meta.role: app
//meta.icon: /img/preclinical_case_icon.png
export async function PreclinicalCaseApp() : Promise<any> {
  return await PackageFunctions.PreclinicalCaseApp();
}

//input: dynamic treeNode 
//meta.role: appTreeBrowser
//meta.app: Preclinical Case
export async function preclinicalCaseAppTreeBrowser(treeNode: any) : Promise<void> {
  await PackageFunctions.preclinicalCaseAppTreeBrowser(treeNode);
}

//input: list file 
//output: list<dataframe> result
//meta.role: fileHandler
//meta.ext: xpt
export async function xptFileHandler(file: DG.FileInfo) : Promise<any> {
  return await PackageFunctions.xptFileHandler(file);
}

//name: Run CDISC CORE Validation
//description: Run CDISC CORE validation on SEND datasets
//input: string standard 
//input: string dataPath 
//input: string version 
//input: string outputFormat 
//input: dynamic options 
//output: string result
export async function runCoreValidate(standard: string, dataPath: string, version?: string, outputFormat?: string, options?: any) : Promise<string> {
  return await PackageFunctions.runCoreValidate(standard, dataPath, version, outputFormat, options);
}
//name: combinedMeasurementsRenderer
//output: grid_cell_renderer renderer
//meta.role: cellRenderer
//meta.cellType: combined-measurements
export function _CombinedMeasurementsCellRenderer() {
  return new CombinedMeasurementsCellRenderer();
}

//name: sdiscRuleViolationRenderer
//output: grid_cell_renderer renderer
//meta.role: cellRenderer
//meta.cellType: sdisc-rule-violation
export function _SdiscRuleViolationCellRenderer() {
  return new SdiscRuleViolationCellRenderer();
}

//name: yAxisRenderer
//output: grid_cell_renderer renderer
//meta.role: cellRenderer
//meta.cellType: y-axis
export function _YAxisCellRenderer() {
  return new YAxisCellRenderer();
}

