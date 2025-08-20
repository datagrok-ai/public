import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: compareColumns
//top-menu: Data | Compare Columns...
export function _compareColumns() {
  return PackageFunctions._compareColumns();
}

//name: welcomeView
//output: view result
//meta.autostartImmediate: true
export function _welcomeView() {
  return PackageFunctions._welcomeView();
}

//name: Activity dashboard
//tags: dashboard
//output: widget result
//meta.showName: false
//meta.order: 1
export function activityDashboardWidget() {
  return PackageFunctions.activityDashboardWidget();
}

//name: Recent projects
//tags: dashboard
//output: widget result
//meta.order: 2
export function recentProjectsWidget() {
  return PackageFunctions.recentProjectsWidget();
}

//name: Community
//tags: dashboard
//output: widget result
//meta.order: 6
export function communityWidget() {
  return PackageFunctions.communityWidget();
}

//name: webWidget
//output: widget result
export function webWidget() {
  return PackageFunctions.webWidget();
}

//name: htmlWidget
//output: widget result
export function htmlWidget() {
  return PackageFunctions.htmlWidget();
}

//name: Learn
//tags: dashboard
//output: widget result
//meta.order: 5
export function learnWidget() {
  return PackageFunctions.learnWidget();
}

//name: kpiWidget
//output: widget result
export function kpiWidget() {
  return PackageFunctions.kpiWidget();
}

//name: isFormulaColumn
//input: column col 
//output: bool result
export function isFormulaColumn(col: DG.Column) {
  return PackageFunctions.isFormulaColumn(col);
}

//name: Formula
//tags: panel
//input: column col 
//output: widget result
//condition: PowerPack:isFormulaColumn(col)
export function formulaWidget(col: DG.Column) {
  return PackageFunctions.formulaWidget(col);
}

//name: powerPackSearchProvider
//tags: searchProvider
//output: dynamic result
export function powerPackSearchProvider() {
  return PackageFunctions.powerPackSearchProvider();
}

//name: formulaLinesEditor
//input: dataframe src { optional: true }
//top-menu: Data | Formula Lines...
export function formulaLinesDialog(src: any) {
  return PackageFunctions.formulaLinesDialog(src);
}

//name: powerPackInit
//tags: init
export async function powerPackInit() {
  return PackageFunctions.powerPackInit();
}

//name: windowsManager
//description: Windows Manager
//tags: autostart
export function windowsManager() {
  return PackageFunctions.windowsManager();
}

//name: viewerDialog
//description: Open 'Viewer Gallery' dialog
//input: dynamic tv 
export function viewerDialog(tv: any) {
  return PackageFunctions.viewerDialog(tv);
}

//name: viewerGallery
//description: ViewerGallery
//tags: autostart
export function viewerGallery() {
  return PackageFunctions.viewerGallery();
}

//name: markdownFileViewer
//tags: fileViewer
//input: file file 
//output: view result
//meta.fileViewer: md,mdx
export async function markdownFileViewer(file: DG.FileInfo) {
  return PackageFunctions.markdownFileViewer(file);
}

//name: xlsxFileHandler
//description: Opens Excel file
//tags: file-handler
//input: list bytes 
//input: string sheetName { optional: true }
//output: list<dataframe> result
//meta.ext: xlsx
export async function xlsxFileHandler(bytes: Uint8Array, sheetName?: string) {
  return PackageFunctions.xlsxFileHandler(bytes, sheetName);
}
