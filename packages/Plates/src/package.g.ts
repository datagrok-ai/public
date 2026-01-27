import {PlateGridCellRenderer} from './plate/plate-cell-renderer';
import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Assay Plates
//description: Assasy plates with concentration, layout and readout data
//meta.demoPath: Plates | Assay Plates
export async function assayPlatesDemo() : Promise<void> {
  await PackageFunctions.assayPlatesDemo();
}

//meta.role: init
export async function _initPlates() : Promise<void> {
  await PackageFunctions._initPlates();
}

//input: file folder 
//input: list<file> files 
//output: dynamic result
//meta.role: folderViewer
export async function platesFolderPreview(folder: DG.FileInfo, files: DG.FileInfo[]) {
  return await PackageFunctions.platesFolderPreview(folder, files);
}

//input: file file 
//output: view result
//meta.role: fileViewer
//meta.fileViewer: txt
//meta.fileViewerCheck: Plates:checkFileIsPlate
export function previewPlate(file: DG.FileInfo) : any {
  return PackageFunctions.previewPlate(file);
}

//input: string fileContent 
//output: list<dataframe> result
//meta.role: fileHandler
//meta.ext: txt
//meta.fileViewerCheck: Plates:checkFileIsPlate
export async function importPlate(fileContent: string) : Promise<any> {
  return await PackageFunctions.importPlate(fileContent);
}

//description: Checks if an Excel file contains plate data.
//input: blob content 
//output: bool result
export async function checkExcelIsPlate(content: Uint8Array) : Promise<boolean> {
  return await PackageFunctions.checkExcelIsPlate(content);
}

//input: blob fileContent 
//meta.role: fileHandler
//meta.ext: xlsx
//meta.fileViewerCheck: Plates:checkExcelIsPlate
export async function importPlateXlsx(fileContent: Uint8Array) : Promise<any> {
  return await PackageFunctions.importPlateXlsx(fileContent);
}

//name: viewPlateXlsx
//input: file file 
//output: view result
//meta.role: fileViewer
//meta.fileViewer: xlsx
//meta.fileViewerCheck: Plates:checkExcelIsPlate
export async function previewPlateXlsx(file: DG.FileInfo) : Promise<any> {
  return await PackageFunctions.previewPlateXlsx(file);
}

//description: Checks if a CSV file can be parsed as a plate.
//input: file file 
//output: bool result
export async function checkCsvIsPlate(file: DG.FileInfo) : Promise<boolean> {
  return await PackageFunctions.checkCsvIsPlate(file);
}

//input: string content 
//output: bool result
export function checkFileIsPlate(content: string) : boolean {
  return PackageFunctions.checkFileIsPlate(content);
}

//name: Plates
//output: view result
//meta.role: app
export function platesApp() : any {
  return PackageFunctions.platesApp();
}

//input: dynamic treeNode 
//meta.role: appTreeBrowser
//meta.app: Plates
export async function platesAppTreeBrowser(treeNode: any) : Promise<void> {
  await PackageFunctions.platesAppTreeBrowser(treeNode);
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
//output: grid_cell_renderer renderer
//meta.role: cellRenderer
//meta.cellType: Plate
export function _PlateGridCellRenderer() {
  return new PlateGridCellRenderer();
}

