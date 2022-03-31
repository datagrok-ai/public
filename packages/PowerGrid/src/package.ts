/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ImageCellRenderer} from "./cell-types/image-cell-renderer";
import {HyperlinkCellRenderer} from "./cell-types/hyperlink-cell-renderer";
import {HtmlTestCellRenderer, TestCellRenderer} from "./cell-types/test-cell-renderer";
import {SparklineCellRenderer} from "./sparklines/sparklines-lines";

export const _package = new DG.Package();

//name: imageUrlCellRenderer
//tags: cellRenderer, cellRenderer-ImageUrl
//meta.cellType: ImageUrl
//output: grid_cell_renderer result
export function imageUrlCellRenderer() {
  return new ImageCellRenderer();
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

//name: sparklineCellRenderer
//tags: cellRenderer, cellRenderer-sparkline
//meta.cellType: sparkline
//meta.virtual: true
//output: grid_cell_renderer result
export function sparklineCellRenderer() {
  return new SparklineCellRenderer();
}

//name: htestCellRenderer
//tags: cellRenderer, cellRenderer-htest
//meta.cellType: htest
//output: grid_cell_renderer result
export function htestCellRenderer() {
  return new HtmlTestCellRenderer();
}