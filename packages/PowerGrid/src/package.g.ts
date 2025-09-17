import {ScatterPlotCellRenderer} from './sparklines/scatter-plot';
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


//name: barCellRenderer
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

//name: summarizeColumns
//description: Adds a sparkline column for the selected columns
//input: list columns { type: numerical }
//meta.action: Sparklines...
export function summarizeColumns(columns: DG.Column[]) {
  return PackageFunctions.summarizeColumns(columns);
}

//name: addFormColumn
//description: Adds a 'form' column for the selected columns
//input: list columns 
//meta.action: Smart form...
export function addFormColumn(columns: DG.Column[]) {
  return PackageFunctions.addFormColumn(columns);
}

//name: testUnitsKgCellRenderer
//tags: cellRenderer
//output: grid_cell_renderer result
//meta.cellType: testUnitsKg
//meta.columnTags: foo=bar,units=kg
export function testUnitsKgCellRenderer() {
  return PackageFunctions.testUnitsKgCellRenderer();
}

//name: testUnitsTonCellRenderer
//tags: cellRenderer
//output: grid_cell_renderer result
//meta.cellType: testUnitsTon
//meta.columnTags: foo=bar,units=ton
export function testUnitsTonCellRenderer() {
  return PackageFunctions.testUnitsTonCellRenderer();
}

//name: addPinnedColumn
//input: object gridCol 
//output: object result
export function addPinnedColumn(gridCol: any) {
  return PackageFunctions.addPinnedColumn(gridCol);
}

//name: demoTestUnitsCellRenderer
export function demoTestUnitsCellRenderer() {
  return PackageFunctions.demoTestUnitsCellRenderer();
}

//name: _autoPowerGrid
//tags: autostart
export async function _autoPowerGrid() {
  return PackageFunctions._autoPowerGrid();
}

//name: Forms
//description: Forms viewer
//tags: viewer
//output: viewer result
//meta.icon: files/icons/formviewer.svg
//meta.viewerPosition: bottom
export function formsViewer() {
  return PackageFunctions.formsViewer();
}

//name: Content
//description: Image content
//tags: panel, powergrid, widgets
//input: string imageUrl { semType: ImageUrl }
//output: widget result
export function imgContent(imageUrl: string) {
  return PackageFunctions.imgContent(imageUrl);
}

//name: demoCellTypes
export function demoCellTypes() {
  return PackageFunctions.demoCellTypes();
}

//name: _scWebGPURender
//tags: scWebGPURender
//input: dynamic sc 
//input: bool show 
export async function _scWebGPURender(sc: any, show: boolean) {
  return PackageFunctions._scWebGPURender(sc, show);
}

//name: _scWebGPUPointHitTest
//tags: scWebGPUPointHitTest
//input: dynamic sc 
//input: dynamic pt 
//output: int result
export async function _scWebGPUPointHitTest(sc: any, pt: any) {
  return PackageFunctions._scWebGPUPointHitTest(sc, pt);
}

//name: isWebGPUAvailable
//tags: isWebGPUAvailable
//output: bool result
export function isWebGPUAvailable() {
  return PackageFunctions.isWebGPUAvailable();
}

//name: isWebGPURenderValid
//tags: isWebGPURenderValid
//input: dynamic sc 
//output: bool result
export function isWebGPURenderValid(sc: any) {
  return PackageFunctions.isWebGPURenderValid(sc);
}
//name: Scatter Plot
//tags: cellRenderer
//output: grid_cell_renderer renderer
//meta.cellType: scatterplot
//meta.virtual: true
export function _ScatterPlotCellRenderer() {
  return new ScatterPlotCellRenderer();
}

