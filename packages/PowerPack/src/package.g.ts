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
//tags: dashboard
//output: widget result
//meta.showName: false
//meta.order: 1
export function activityDashboardWidget() : any {
  return PackageFunctions.activityDashboardWidget();
}

//name: Recent projects
//tags: dashboard
//output: widget result
//meta.order: 2
export function recentProjectsWidget() : any {
  return PackageFunctions.recentProjectsWidget();
}

//name: Community
//tags: dashboard
//output: widget result
//meta.order: 6
export function communityWidget() : any {
  return PackageFunctions.communityWidget();
}

//name: webWidget
//output: widget result
export function webWidget() : any {
  return PackageFunctions.webWidget();
}

//name: htmlWidget
//output: widget result
export function htmlWidget() : any {
  return PackageFunctions.htmlWidget();
}

//name: Learn
//tags: dashboard
//output: widget result
//meta.order: 5
export function learnWidget() : any {
  return PackageFunctions.learnWidget();
}

//name: kpiWidget
//output: widget result
export function kpiWidget() : any {
  return PackageFunctions.kpiWidget();
}

//name: isFormulaColumn
//input: column col 
//output: bool result
export function isFormulaColumn(col: DG.Column) : boolean {
  return PackageFunctions.isFormulaColumn(col);
}

//name: Formula
//tags: panel
//input: column col 
//output: widget result
//condition: PowerPack:isFormulaColumn(col)
export function formulaWidget(col: DG.Column) : any {
  return PackageFunctions.formulaWidget(col);
}

//name: powerPackSearchProvider
//tags: searchProvider
//output: dynamic result
export function powerPackSearchProvider() : any {
  return PackageFunctions.powerPackSearchProvider();
}

//name: formulaLinesEditor
//input: dataframe src { optional: true }
export function formulaLinesDialog(src: any) : void {
  PackageFunctions.formulaLinesDialog(src);
}

//name: powerPackInit
//tags: init
export async function powerPackInit() : Promise<void> {
  PackageFunctions.powerPackInit();
}

//name: windowsManager
//description: Windows Manager
//tags: autostart
export function windowsManager() : void {
  PackageFunctions.windowsManager();
}

//name: viewerDialog
//description: Open 'Viewer Gallery' dialog
//input: dynamic tv 
export function viewerDialog(tv: any) : void {
  PackageFunctions.viewerDialog(tv);
}

//name: viewerGallery
//description: ViewerGallery
//tags: autostart
export function viewerGallery() : void {
  PackageFunctions.viewerGallery();
}

//name: markdownFileViewer
//tags: fileViewer
//input: file file 
//output: view result
//meta.fileViewer: md,mdx
export async function markdownFileViewer(file: DG.FileInfo) : Promise<any> {
  return PackageFunctions.markdownFileViewer(file);
}

//name: xlsxFileHandler
//description: Opens Excel file
//tags: file-handler
//input: list bytes 
//input: string sheetName { optional: true }
//output: list<dataframe> result
//meta.ext: xlsx
export async function xlsxFileHandler(bytes: Uint8Array, sheetName?: string) : Promise<any> {
  return PackageFunctions.xlsxFileHandler(bytes, sheetName);
}
