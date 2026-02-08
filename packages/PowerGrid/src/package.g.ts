import {ScatterPlotCellRenderer} from './sparklines/scatter-plot';
import {RawPNGRenderer} from './pngRenderer';
import {PackageFunctions} from './package';
import {HtmlTestCellRenderer} from './cell-types/test-cell-renderer';
import {StarsCellRenderer} from './cell-types/stars-cell-renderer';
import {MultiChoiceCellRenderer} from './cell-types/multi-choice-cell-renderer';
import {ImageCellRenderer} from './cell-types/image-cell-renderer';
import {HyperlinkCellRenderer} from './cell-types/hyperlink-cell-renderer';
import {BinaryImageCellRenderer} from './cell-types/binary-image-cell-renderer';
import * as DG from 'datagrok-api/dg';
//name: binaryImageCellRenderer
//output: grid_cell_renderer renderer
//meta.role: cellRenderer
//meta.cellType: BinaryImage
export function _BinaryImageCellRenderer() {
  return new BinaryImageCellRenderer();
}

//name: hyperlinkCellRenderer
//output: grid_cell_renderer renderer
//meta.role: cellRenderer
//meta.cellType: Hyperlink
export function _HyperlinkCellRenderer() {
  return new HyperlinkCellRenderer();
}

//name: imageUrlCellRenderer
//output: grid_cell_renderer renderer
//meta.role: cellRenderer
//meta.cellType: ImageUrl
export function _ImageCellRenderer() {
  return new ImageCellRenderer();
}

//name: Multi Choice
//output: grid_cell_renderer renderer
//meta.role: cellRenderer
//meta.cellType: MultiChoice
export function _MultiChoiceCellRenderer() {
  return new MultiChoiceCellRenderer();
}

//name: Stars
//output: grid_cell_renderer renderer
//meta.role: cellRenderer
//meta.cellType: Stars
export function _StarsCellRenderer() {
  return new StarsCellRenderer();
}

//name: htestCellRenderer
//output: grid_cell_renderer renderer
//meta.role: cellRenderer
//meta.cellType: htest
export function _HtmlTestCellRenderer() {
  return new HtmlTestCellRenderer();
}


//output: grid_cell_renderer result
//meta.gridChart: true
//meta.cellType: bar
//meta.role: cellRenderer
export function barCellRenderer() {
  return PackageFunctions.barCellRenderer();
}

//name: Sparklines
//output: grid_cell_renderer result
//meta.cellType: sparkline
//meta.gridChart: true
//meta.virtual: true
//meta.role: cellRenderer
export function sparklineCellRenderer() {
  return PackageFunctions.sparklineCellRenderer();
}

//name: Bar Chart
//output: grid_cell_renderer result
//meta.cellType: barchart
//meta.gridChart: true
//meta.virtual: true
//meta.role: cellRenderer
export function barchartCellRenderer() {
  return PackageFunctions.barchartCellRenderer();
}

//name: Pie Chart
//output: grid_cell_renderer result
//meta.cellType: piechart
//meta.gridChart: true
//meta.virtual: true
//meta.role: cellRenderer
export function piechartCellRenderer() {
  return PackageFunctions.piechartCellRenderer();
}

//name: Radar
//output: grid_cell_renderer result
//meta.cellType: radar
//meta.gridChart: true
//meta.virtual: true
//meta.role: cellRenderer
export function radarCellRenderer() {
  return PackageFunctions.radarCellRenderer();
}

//name: Smart Form
//output: grid_cell_renderer result
//meta.cellType: smartform
//meta.gridChart: true
//meta.virtual: true
//meta.role: cellRenderer
export function smartFormCellRenderer() {
  return PackageFunctions.smartFormCellRenderer();
}

//name: Tags
//output: grid_cell_renderer result
//meta.cellType: Tags
//meta.gridChart: true
//meta.virtual: true
//meta.role: cellRenderer
export function tagsCellRenderer() {
  return PackageFunctions.tagsCellRenderer();
}

//name: Confidence Interval
//output: grid_cell_renderer result
//meta.cellType: ConfidenceInterval
//meta.gridChart: true
//meta.virtual: true
//meta.role: cellRenderer
export function confidenceIntervalCellRenderer() {
  return PackageFunctions.confidenceIntervalCellRenderer();
}

//description: Adds a sparkline column for the selected columns
//input: list columns { type: numerical }
//meta.action: Sparklines...
export function summarizeColumns(columns: DG.Column[]) : void {
  PackageFunctions.summarizeColumns(columns);
}

//description: Adds a 'form' column for the selected columns
//input: list columns 
//meta.action: Smart form...
export function addFormColumn(columns: DG.Column[]) : void {
  PackageFunctions.addFormColumn(columns);
}

//output: grid_cell_renderer result
//meta.cellType: testUnitsKg
//meta.columnTags: foo=bar,units=kg
//meta.role: cellRenderer
export function testUnitsKgCellRenderer() {
  return PackageFunctions.testUnitsKgCellRenderer();
}

//output: grid_cell_renderer result
//meta.cellType: testUnitsTon
//meta.columnTags: foo=bar,units=ton
//meta.role: cellRenderer
export function testUnitsTonCellRenderer() {
  return PackageFunctions.testUnitsTonCellRenderer();
}

//input: object gridCol 
//output: object result
export function addPinnedColumn(gridCol: any) : any {
  return PackageFunctions.addPinnedColumn(gridCol);
}

//name: demoTestUnitsCellRenderer
export function demoTestUnitsCellRenderer() : void {
  PackageFunctions.demoTestUnitsCellRenderer();
}

//meta.role: autostart
export async function _autoPowerGrid() : Promise<void> {
  await PackageFunctions._autoPowerGrid();
}

//name: Forms
//description: Forms viewer
//output: viewer result
//meta.icon: files/icons/formviewer.svg
//meta.viewerPosition: bottom
//meta.toolbox: true
//meta.role: viewer
export function formsViewer() {
  return PackageFunctions.formsViewer();
}

//name: Content
//description: Image content
//input: string imageUrl { semType: ImageUrl }
//output: widget result
//meta.role: widgets,panel
export function imgContent(imageUrl: string) : any {
  return PackageFunctions.imgContent(imageUrl);
}

//name: demoCellTypes
export function demoCellTypes() : void {
  PackageFunctions.demoCellTypes();
}

//input: dynamic sc 
//input: bool show 
export async function _scWebGPURender(sc: any, show: boolean) : Promise<void> {
  await PackageFunctions._scWebGPURender(sc, show);
}

//input: dynamic sc 
//input: dynamic pt 
//output: int result
export async function _scWebGPUPointHitTest(sc: any, pt: any) {
  return await PackageFunctions._scWebGPUPointHitTest(sc, pt);
}

//output: bool result
export function isWebGPUAvailable() : boolean {
  return PackageFunctions.isWebGPUAvailable();
}

//input: dynamic sc 
//output: bool result
export function isWebGPURenderValid(sc: any) : boolean {
  return PackageFunctions.isWebGPURenderValid(sc);
}
//name: rawPng
//output: grid_cell_renderer renderer
//meta.role: cellRenderer
//meta.cellType: rawPng
export function _RawPNGRenderer() {
  return new RawPNGRenderer();
}

//name: Scatter Plot
//output: grid_cell_renderer renderer
//meta.role: cellRenderer
//meta.cellType: scatterplot
//meta.virtual: true
export function _ScatterPlotCellRenderer() {
  return new ScatterPlotCellRenderer();
}

