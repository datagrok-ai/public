import {PlateGridCellRenderer} from './plate/plate-cell-renderer';
import {PackageFunctions} from './package';
import {MultiCurveViewer} from './fit/multi-curve-viewer';
import {FitChartCellRenderer} from './fit/fit-renderer';

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
export function _MultiCurveViewer() {
  return new MultiCurveViewer();
}

//name: Curve fitting
//description: Curve fitting is the process of constructing a curve, or mathematical function, that has the best fit to a series of data points
//meta.demoPath: Curves | Curve Fitting
//test: curveFitDemo() //wait: 2000 
export async function curveFitDemo() {
  return PackageFunctions.curveFitDemo();
}

//name: Assay Plates
//description: Assasy plates with concentration, layout and readout data
//meta.demoPath: Curves | Assay Plates
export async function assayPlatesDemo() {
  return PackageFunctions.assayPlatesDemo();
}

//name: _initCurves
//tags: init
export function _initCurves() {
  return PackageFunctions._initCurves();
}

//name: addStatisticsColumn
//tags: Transform
//input: dataframe df 
//input: string colName 
//input: string propName 
//input: string seriesName 
//input: double seriesNumber 
//input: string newColName 
export function addStatisticsColumn(df: any, colName: string, propName: string, seriesName: string, seriesNumber: number, newColName: string) {
  return PackageFunctions.addStatisticsColumn(df, colName, propName, seriesName, seriesNumber, newColName);
}

//name: addAggrStatisticsColumn
//tags: Transform
//input: dataframe df 
//input: string colName 
//input: string propName 
//input: string aggrType 
export function addAggrStatisticsColumn(df: any, colName: string, propName: string, aggrType: string) {
  return PackageFunctions.addAggrStatisticsColumn(df, colName, propName, aggrType);
}

//name: platesFolderPreview
//tags: folderViewer
//input: file folder 
//input: list<file> files 
//output: widget result
export async function platesFolderPreview(folder: any, files: any) {
  return PackageFunctions.platesFolderPreview(folder, files);
}

//name: previewPlate
//tags: fileViewer
//input: file file 
//output: view result
//meta.fileViewer: txt
//meta.fileViewerCheck: Curves:checkFileIsPlate
export function previewPlate(file: any) {
  return PackageFunctions.previewPlate(file);
}

//name: importPlate
//tags: file-handler
//input: string fileContent 
//output: list<dataframe> result
//meta.ext: txt
//meta.fileViewerCheck: Curves:checkFileIsPlate
export async function importPlate(fileContent: string) {
  return PackageFunctions.importPlate(fileContent);
}

//name: importPlateXlsx
//tags: file-handler
//input: blob fileContent 
//meta.ext: xlsx
//meta.fileViewerCheck: Curves:checkExcelIsPlate
export async function importPlateXlsx(fileContent: any) {
  return PackageFunctions.importPlateXlsx(fileContent);
}

//name: viewPlateXlsx
//tags: fileViewer
//input: file file 
//output: view result
//meta.fileViewer: xlsx
//meta.fileViewerCheck: Curves:checkExcelIsPlate
export async function previewPlateXlsx(file: any) {
  return PackageFunctions.previewPlateXlsx(file);
}

//name: checkExcelIsPlate
//input: blob content 
//output: bool result
export async function checkExcelIsPlate(content: any) {
  return PackageFunctions.checkExcelIsPlate(content);
}

//name: checkFileIsPlate
//input: string content 
//output: bool result
export function checkFileIsPlate(content: string) {
  return PackageFunctions.checkFileIsPlate(content);
}

//name: Browse
//tags: app
//output: view result
//meta.browsePath: Plates
export function platesApp() {
  return PackageFunctions.platesApp();
}

//name: platesAppTreeBrowser
//input: dynamic treeNode 
export async function platesAppTreeBrowser(treeNode: any) {
  return PackageFunctions.platesAppTreeBrowser(treeNode);
}

//name: getPlateByBarcode
//input: string barcode 
//output: dynamic result
export async function getPlateByBarcode(barcode: string) {
  return PackageFunctions.getPlateByBarcode(barcode);
}

//name: createDummyPlateData
export async function createDummyPlateData() {
  return PackageFunctions.createDummyPlateData();
}

//name: PlateGridCellRenderer
//tags: cellRenderer
//output: grid_cell_renderer renderer
//meta.cellType: Plate
export function _PlateGridCellRenderer() {
  return new PlateGridCellRenderer();
}

