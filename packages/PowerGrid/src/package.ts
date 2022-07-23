/* Do not change these import lines to match external modules in webpack configuration */
// import * as grok from 'datagrok-api/grok';
// import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {ImageCellRenderer} from './cell-types/image-cell-renderer';
import {HyperlinkCellRenderer} from './cell-types/hyperlink-cell-renderer';
import {HtmlTestCellRenderer, TestCellRenderer} from './cell-types/test-cell-renderer';
import {BarCellRenderer} from './cell-types/bar-cell-renderer';
import {BinaryImageCellRenderer} from './cell-types/binary-image-cell-renderer';

import {SparklineCellRenderer} from './sparklines/sparklines-lines';
import {BarChartCellRenderer} from './sparklines/bar-chart';
import {PieChartCellRenderer} from './sparklines/piechart';
import {RadarChartCellRender} from './sparklines/radar-chart';

export const _package = new DG.Package();

//name: imageUrlCellRenderer
//tags: cellRenderer
//meta.cellType: ImageUrl
//output: grid_cell_renderer result
export function imageUrlCellRenderer() {
  return new ImageCellRenderer();
}

//name: binaryImageCellRenderer
//tags: cellRenderer
//meta.cellType: BinaryImage
//output: grid_cell_renderer result
export function binaryImageCellRenderer() {
  return new BinaryImageCellRenderer();
}

//name: hyperlinkCellRenderer
//tags: cellRenderer
//meta.cellType: Hyperlink
//output: grid_cell_renderer result
export function hyperlinkCellRenderer() {
  return new HyperlinkCellRenderer();
}

//name: testCellRenderer
//tags: cellRenderer
//meta.cellType: test
//meta.columnTag.units: bar
//output: grid_cell_renderer result
export function testCellRenderer() {
  return new TestCellRenderer();
}

//name: htestCellRenderer
//tags: cellRenderer
//meta.cellType: htest
//output: grid_cell_renderer result
export function htestCellRenderer() {
  return new HtmlTestCellRenderer();
}

//name: barCellRenderer
//tags: cellRenderer
//meta.cellType: bar
//output: grid_cell_renderer result
export function barCellRenderer() {
  return new BarCellRenderer();
}

//name: Sparklines
//tags: cellRenderer
//meta.cellType: sparkline
//meta.virtual: true
//output: grid_cell_renderer result
export function sparklineCellRenderer() {
  return new SparklineCellRenderer();
}

//name: Bar Chart
//tags: cellRenderer
//meta.cellType: barchart
//meta.virtual: true
//output: grid_cell_renderer result
export function barchartCellRenderer() {
  return new BarChartCellRenderer();
}

//name: Pie Chart
//tags: cellRenderer
//meta.cellType: piechart
//meta.virtual: true
//output: grid_cell_renderer result
export function piechartCellRenderer() {
  return new PieChartCellRenderer();
}

//name: Radar
//tags: cellRenderer
//meta.cellType: radar
//meta.virtual: true
//output: grid_cell_renderer result
export function radarCellRenderer() {
  return new RadarChartCellRender();
}
