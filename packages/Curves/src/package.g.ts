import {RawPNGRenderer} from './pngRenderer';
import {PlateGridCellRenderer} from './plate/plate-cell-renderer';
import {PackageFunctions} from './package';
import {MultiCurveViewer} from './fit/multi-curve-viewer';
import {FitChartCellRenderer} from './fit/fit-renderer';
import * as DG from 'datagrok-api/dg';
//name: Fit
//tags: cellRenderer
//output: grid_cell_renderer renderer
//meta.cellType: fit
//meta.virtual: true
export function _FitChartCellRenderer() {
  return new FitChartCellRenderer();
}

//name: MultiCurveViewer
//description: A viewer that superimposes multiple in-cell curves on one chart
//tags: viewer
//output: viewer result
//meta.icon: icons/multi-curve-viewer.png
//meta.trellisable: true
export function _MultiCurveViewer() {
  return new MultiCurveViewer();
}


//name: Curve fitting
//description: Curve fitting is the process of constructing a curve, or mathematical function, that has the best fit to a series of data points
//meta.demoPath: Curves | Curve Fitting
//test: curveFitDemo() //wait: 2000 
export async function curveFitDemo() : Promise<void> {
  await PackageFunctions.curveFitDemo();
}

//tags: init
export function _initCurves() : void {
  PackageFunctions._initCurves();
}

//input: dataframe df 
//input: column concentrationCol 
//input: column readoutCol 
//input: column batchIDCol 
//input: column assayCol 
//input: column runIDCol 
//input: column compoundIDCol 
//input: column targetEntityCol 
//input: column excludeOutliersCol { nullable: true }
//input: dataframe parentTable { nullable: true }
//input: list<string> fitParamColumns { nullable: true }
//input: string reportedIC50Column { nullable: true }
//input: string reportedQualifiedIC50Column { nullable: true }
//input: string experimentIDColumn { nullable: true }
//input: string qualifierColumn { nullable: true }
//input: list<string> additionalColumns { nullable: true }
//input: string wellLevelJoinCol { nullable: true }
//input: string parentLevelJoinCol { nullable: true }
//output: dataframe result
export async function dataToCurves(df: DG.DataFrame, concentrationCol: DG.Column, readoutCol: DG.Column, batchIDCol: DG.Column, assayCol: DG.Column, runIDCol: DG.Column, compoundIDCol: DG.Column, targetEntityCol: DG.Column, excludeOutliersCol?: DG.Column, parentTable?: DG.DataFrame, fitParamColumns?: string[], reportedIC50Column?: string, reportedQualifiedIC50Column?: string, experimentIDColumn?: string, qualifierColumn?: string, additionalColumns?: string[], wellLevelJoinCol?: string, parentLevelJoinCol?: string) : Promise<any> {
  return await PackageFunctions.dataToCurves(df, concentrationCol, readoutCol, batchIDCol, assayCol, runIDCol, compoundIDCol, targetEntityCol, excludeOutliersCol, parentTable, fitParamColumns, reportedIC50Column, reportedQualifiedIC50Column, experimentIDColumn, qualifierColumn, additionalColumns, wellLevelJoinCol, parentLevelJoinCol);
}

//top-menu: Data | Curves | Data to Curves
export async function dataToCurvesTopMenu() : Promise<void> {
  await PackageFunctions.dataToCurvesTopMenu();
}

//tags: Transform
//input: dataframe table 
//input: string colName 
//input: string propName 
//input: int seriesNumber 
//output: column result
//meta.vectorFunc: true
export function addStatisticsColumn(table: DG.DataFrame, colName: string, propName: string, seriesNumber: number) : any {
  return PackageFunctions.addStatisticsColumn(table, colName, propName, seriesNumber);
}

//tags: Transform
//input: dataframe table 
//input: string colName 
//input: string propName 
//input: string aggrType 
//output: column result
//meta.vectorFunc: true
export function addAggrStatisticsColumn(table: DG.DataFrame, colName: string, propName: string, aggrType: string) : any {
  return PackageFunctions.addAggrStatisticsColumn(table, colName, propName, aggrType);
}

//tags: folderViewer
//input: file folder 
//input: list<file> files 
export async function platesFolderPreview(folder: DG.FileInfo, files: DG.FileInfo[]) : Promise<void> {
  await PackageFunctions.platesFolderPreview(folder, files);
}

//tags: fileViewer
//input: file file 
//output: view result
//meta.fileViewer: txt
//meta.fileViewerCheck: Curves:checkFileIsPlate
export function previewPlate(file: DG.FileInfo) : any {
  return PackageFunctions.previewPlate(file);
}

//tags: file-handler
//input: string fileContent 
//output: list<dataframe> result
//meta.ext: txt
//meta.fileViewerCheck: Curves:checkFileIsPlate
export async function importPlate(fileContent: string) : Promise<any> {
  return await PackageFunctions.importPlate(fileContent);
}

//description: Checks if a CSV file can be parsed as a plate.
//input: file file 
//output: bool result
export async function checkCsvIsPlate(file: DG.FileInfo) : Promise<boolean> {
  return await PackageFunctions.checkCsvIsPlate(file);
}

//name: Browse
//tags: app
//output: view result
//meta.browsePath: Plates
export function platesApp() : any {
  return PackageFunctions.platesApp();
}

//input: string barcode 
//output: dynamic result
export async function getPlateByBarcode(barcode: string) : Promise<any> {
  return await PackageFunctions.getPlateByBarcode(barcode);
}

//name: createDummyPlateData
export async function createDummyPlateData() : Promise<void> {
  await PackageFunctions.createDummyPlateData();
}
//name: PlateGridCellRenderer
//tags: cellRenderer
//output: grid_cell_renderer renderer
//meta.cellType: Plate
export function _PlateGridCellRenderer() {
  return new PlateGridCellRenderer();
}

//name: rawPng
//tags: cellRenderer
//output: grid_cell_renderer renderer
//meta.cellType: rawPng
export function _RawPNGRenderer() {
  return new RawPNGRenderer();
}

