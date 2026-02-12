import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: compareColumns
//top-menu: Data | Compare Columns...
export function _compareColumns() : void {
  PackageFunctions._compareColumns();
}

//name: welcomeView
//output: view home
//meta.autostartImmediate: true
export function _welcomeView() : any {
  return PackageFunctions._welcomeView();
}

//name: Activity dashboard
//output: widget result
//meta.showName: false
//meta.role: dashboard
//meta.order: -1
export function activityDashboardWidget() : any {
  return PackageFunctions.activityDashboardWidget();
}

//name: Community
//output: widget result
//meta.role: dashboard
//meta.order: 6
export function communityWidget() : any {
  return PackageFunctions.communityWidget();
}

//output: widget result
export function webWidget() : any {
  return PackageFunctions.webWidget();
}

//output: widget result
export function htmlWidget() : any {
  return PackageFunctions.htmlWidget();
}

//output: widget result
export function kpiWidget() : any {
  return PackageFunctions.kpiWidget();
}

//output: object result
//meta.propertyType: string
//meta.semType: cron
//meta.role: valueEditor
export function cronInput() : any {
  return PackageFunctions.cronInput();
}

//input: column col 
//output: bool result
export function isFormulaColumn(col: DG.Column) : boolean {
  return PackageFunctions.isFormulaColumn(col);
}

//name: Formula
//input: column col 
//output: widget result
//meta.role: panel
//condition: PowerPack:isFormulaColumn(col)
export function formulaWidget(col: DG.Column) : any {
  return PackageFunctions.formulaWidget(col);
}

//input: dynamic func 
//input: dynamic inputParams 
//output: widget result
export function getFuncTableViewWidget(func: any, inputParams: any) : any {
  return PackageFunctions.getFuncTableViewWidget(func, inputParams);
}

//output: dynamic result
//meta.role: searchProvider
export function powerPackSearchProvider() : any {
  return PackageFunctions.powerPackSearchProvider();
}

//name: formulaLinesEditor
//input: dataframe src { optional: true }
//input: int currentIndexToSet { optional: true }
//input: bool isDataFrameValue { optional: true }
//input: bool isAnnotationArea { optional: true }
export function formulaLinesDialog(src: any, currentIndexToSet?: number, isDataFrameValue?: boolean, isAnnotationArea?: boolean) : void {
  PackageFunctions.formulaLinesDialog(src, currentIndexToSet, isDataFrameValue, isAnnotationArea);
}

//meta.role: init
export async function powerPackInit() : Promise<void> {
  await PackageFunctions.powerPackInit();
}

//description: Windows Manager
//meta.role: autostart
export function windowsManager() : void {
  PackageFunctions.windowsManager();
}

//description: Open 'Viewer Gallery' dialog
//input: dynamic tv 
export function viewerDialog(tv: any) : void {
  PackageFunctions.viewerDialog(tv);
}

//description: ViewerGallery
//meta.role: autostart
export function viewerGallery() : void {
  PackageFunctions.viewerGallery();
}

//input: file file 
//output: view result
//meta.role: fileViewer
//meta.fileViewer: md,mdx
export async function markdownFileViewer(file: DG.FileInfo) : Promise<any> {
  return await PackageFunctions.markdownFileViewer(file);
}

//description: Opens Excel file
//input: list bytes 
//input: string sheetName { optional: true }
//output: list<dataframe> result
//meta.role: fileHandler
//meta.ext: xlsx
export async function xlsxFileHandler(bytes: Uint8Array, sheetName?: string) : Promise<any> {
  return await PackageFunctions.xlsxFileHandler(bytes, sheetName);
}

//input: dynamic conn 
//input: string schema 
//input: string table 
//input: string column 
//input: string name 
//input: dataframe df 
//input: string db 
//meta.role: transform
export async function runEnrichment(conn: any, schema: string, table: string, column: string, name: string, df: DG.DataFrame, db: string) : Promise<void> {
  await PackageFunctions.runEnrichment(conn, schema, table, column, name, df, db);
}
