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
import {BarChartCellRenderer} from './sparklines/barchart';
import {PieChartCellRenderer} from './sparklines/piechart';
import {RadarChartCellRender} from './sparklines/radarchart';

export const _package = new DG.Package();

//name: imageUrlCellRenderer
//tags: cellRenderer, cellRenderer-ImageUrl
//meta.cellType: ImageUrl
//output: grid_cell_renderer result
export function imageUrlCellRenderer() {
  return new ImageCellRenderer();
}

//name: binaryImageCellRenderer
//tags: cellRenderer, cellRenderer-BinaryImage
//meta.cellType: BinaryImage
//output: grid_cell_renderer result
export function binaryImageCellRenderer() {
  return new BinaryImageCellRenderer();
}

//name: hyperlinkCellRenderer
//tags: cellRenderer, cellRenderer-Hyperlink
//meta.cellType: Hyperlink
//output: grid_cell_renderer result
export function hyperlinkCellRenderer() {
  return new HyperlinkCellRenderer();
}

//name: testCellRenderer
//tags: cellRenderer, cellRenderer-test
//meta.cellType: test
//output: grid_cell_renderer result
export function testCellRenderer() {
  return new TestCellRenderer();
}

//name: htestCellRenderer
//tags: cellRenderer, cellRenderer-htest
//meta.cellType: htest
//output: grid_cell_renderer result
export function htestCellRenderer() {
  return new HtmlTestCellRenderer();
}

//name: barCellRenderer
//tags: cellRenderer, cellRenderer-bar
//meta.cellType: bar
//output: grid_cell_renderer result
export function barCellRenderer() {
  return new BarCellRenderer();
}

//name: sparklineCellRenderer
//tags: cellRenderer, cellRenderer-sparkline
//meta.cellType: sparkline_ts
//meta.virtual: true
//output: grid_cell_renderer result
export function sparklineCellRenderer() {
  return new SparklineCellRenderer();
}

//name: barchartCellRenderer
//tags: cellRenderer, cellRenderer-barchart
//meta.cellType: barchart_ts
//meta.virtual: true
//output: grid_cell_renderer result
export function barchartCellRenderer() {
  return new BarChartCellRenderer();
}

//name: piechartCellRender
//tags: cellRenderer, cellRenderer-piechart
//meta.cellType: piechart_ts
//meta.virtual: true
//output: grid_cell_renderer result
export function piechartCellRenderer() {
  return new PieChartCellRenderer();
}

//name: radarCellRenderer
//tags: cellRenderer, cellRenderer-radar
//meta.cellType: radarchart_ts
//meta.virtual: true
//output: grid_cell_renderer result
export function radarCellRenderer() {
  return new RadarChartCellRender();
}
