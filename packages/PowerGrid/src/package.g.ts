import {ScatterPlotCellRenderer} from './sparklines/scatter-plot';
import {RawPNGRenderer} from './pngRenderer';
import {PackageFunctions} from './package';
import {HtmlTestCellRenderer} from './cell-types/test-cell-renderer';
import {TagsCellRenderer} from './cell-types/tags-cell-renderer';
import {MultiChoiceCellRenderer} from './cell-types/multi-choice-cell-renderer';
import {ImageCellRenderer} from './cell-types/image-cell-renderer';
import {HyperlinkCellRenderer} from './cell-types/hyperlink-cell-renderer';
import {BinaryImageCellRenderer} from './cell-types/binary-image-cell-renderer';
import * as DG from 'datagrok-api/dg';
//name: binaryImageCellRenderer
//tags: cellRenderer
//output: grid_cell_renderer renderer
//meta.cellType: BinaryImage
export function _BinaryImageCellRenderer() {
  return new BinaryImageCellRenderer();
}

//name: hyperlinkCellRenderer
//tags: cellRenderer
//output: grid_cell_renderer renderer
//meta.cellType: Hyperlink
export function _HyperlinkCellRenderer() {
  return new HyperlinkCellRenderer();
}

//name: imageUrlCellRenderer
//tags: cellRenderer
//output: grid_cell_renderer renderer
//meta.cellType: ImageUrl
export function _ImageCellRenderer() {
  return new ImageCellRenderer();
}

//name: Multi Choice
//tags: cellRenderer
//output: grid_cell_renderer renderer
//meta.cellType: MultiChoice
export function _MultiChoiceCellRenderer() {
  return new MultiChoiceCellRenderer();
}

//name: Tags
//tags: cellRenderer
//output: grid_cell_renderer renderer
//meta.cellType: Tags
export function _TagsCellRenderer() {
  return new TagsCellRenderer();
}

//name: htestCellRenderer
//tags: cellRenderer
//output: grid_cell_renderer renderer
//meta.cellType: htest
export function _HtmlTestCellRenderer() {
  return new HtmlTestCellRenderer();
}


//tags: cellRenderer
//output: grid_cell_renderer result
//meta.gridChart: true
//meta.cellType: bar
export function barCellRenderer() {
  return PackageFunctions.barCellRenderer();
}

//name: Sparklines
//tags: cellRenderer
//output: grid_cell_renderer result
//meta.cellType: sparkline
//meta.gridChart: true
//meta.virtual: true
export function sparklineCellRenderer() {
  return PackageFunctions.sparklineCellRenderer();
}

//name: Bar Chart
//tags: cellRenderer
//output: grid_cell_renderer result
//meta.cellType: barchart
//meta.gridChart: true
//meta.virtual: true
export function barchartCellRenderer() {
  return PackageFunctions.barchartCellRenderer();
}

//name: Pie Chart
//tags: cellRenderer
//output: grid_cell_renderer result
//meta.cellType: piechart
//meta.gridChart: true
//meta.virtual: true
export function piechartCellRenderer() {
  return PackageFunctions.piechartCellRenderer();
}

//name: Radar
//tags: cellRenderer
//output: grid_cell_renderer result
//meta.cellType: radar
//meta.gridChart: true
//meta.virtual: true
export function radarCellRenderer() {
  return PackageFunctions.radarCellRenderer();
}

//name: Smart Form
//tags: cellRenderer
//output: grid_cell_renderer result
//meta.cellType: smartform
//meta.gridChart: true
//meta.virtual: true
export function smartFormCellRenderer() {
  return PackageFunctions.smartFormCellRenderer();
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

//tags: cellRenderer
//output: grid_cell_renderer result
//meta.cellType: testUnitsKg
//meta.columnTags: foo=bar,units=kg
export function testUnitsKgCellRenderer() {
  return PackageFunctions.testUnitsKgCellRenderer();
}

//tags: cellRenderer
//output: grid_cell_renderer result
//meta.cellType: testUnitsTon
//meta.columnTags: foo=bar,units=ton
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

//tags: autostart
export async function _autoPowerGrid() : Promise<void> {
  await PackageFunctions._autoPowerGrid();
}

//name: Forms
//description: Forms viewer
//tags: viewer
//output: viewer result
//meta.icon: files/icons/formviewer.svg
//meta.viewerPosition: bottom
//meta.toolbox: true
export function formsViewer() {
  return PackageFunctions.formsViewer();
}

//name: Content
//description: Image content
//tags: panel, powergrid, widgets
//input: string imageUrl { semType: ImageUrl }
//output: widget result
export function imgContent(imageUrl: string) : any {
  return PackageFunctions.imgContent(imageUrl);
}

//name: demoCellTypes
export function demoCellTypes() : void {
  PackageFunctions.demoCellTypes();
}

//tags: scWebGPURender
//input: dynamic sc 
//input: bool show 
export async function _scWebGPURender(sc: any, show: boolean) : Promise<void> {
  await PackageFunctions._scWebGPURender(sc, show);
}

//tags: scWebGPUPointHitTest
//input: dynamic sc 
//input: dynamic pt 
//output: int result
export async function _scWebGPUPointHitTest(sc: any, pt: any) {
  return await PackageFunctions._scWebGPUPointHitTest(sc, pt);
}

//tags: isWebGPUAvailable
//output: bool result
export function isWebGPUAvailable() : boolean {
  return PackageFunctions.isWebGPUAvailable();
}

//tags: isWebGPURenderValid
//input: dynamic sc 
//output: bool result
export function isWebGPURenderValid(sc: any) : boolean {
  return PackageFunctions.isWebGPURenderValid(sc);
}
//name: rawPng
//tags: cellRenderer
//output: grid_cell_renderer renderer
//meta.cellType: rawPng
export function _RawPNGRenderer() {
  return new RawPNGRenderer();
}

//name: Scatter Plot
//tags: cellRenderer
//output: grid_cell_renderer renderer
//meta.cellType: scatterplot
//meta.virtual: true
export function _ScatterPlotCellRenderer() {
  return new ScatterPlotCellRenderer();
}

