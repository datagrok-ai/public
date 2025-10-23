import {PlateGridCellRenderer} from './plate/plate-cell-renderer';
import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//tags: init
export function _initPlates() : void {
  PackageFunctions._initPlates();
}

//tags: folderViewer
//input: file folder 
//input: list<file> files 
//output: dynamic result
export async function platesFolderPreview(folder: DG.FileInfo, files: DG.FileInfo[]) {
  return await PackageFunctions.platesFolderPreview(folder, files);
}

//tags: fileViewer
//input: file file 
//output: view result
//meta.fileViewer: txt
//meta.fileViewerCheck: Plates:checkFileIsPlate
export function previewPlate(file: DG.FileInfo) : any {
  return PackageFunctions.previewPlate(file);
}

//tags: file-handler
//input: string fileContent 
//output: list<dataframe> result
//meta.ext: txt
//meta.fileViewerCheck: Plates:checkFileIsPlate
export async function importPlate(fileContent: string) : Promise<any> {
  return await PackageFunctions.importPlate(fileContent);
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

//name: Browse
//tags: app
//output: view result
//meta.browsePath: Plates
export function platesApp() : any {
  return PackageFunctions.platesApp();
}

//input: dynamic treeNode 
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
//tags: cellRenderer
//output: grid_cell_renderer renderer
//meta.cellType: Plate
export function _PlateGridCellRenderer() {
  return new PlateGridCellRenderer();
}

