import {BinaryImageCellRenderer} from './cell-types/binary-image-cell-renderer';
import {HyperlinkCellRenderer} from './cell-types/hyperlink-cell-renderer';
import {ImageCellRenderer} from './cell-types/image-cell-renderer';
import {MultiChoiceCellRenderer} from './cell-types/multi-choice-cell-renderer';
import {TagsCellRenderer} from './cell-types/tags-cell-renderer';
import {HtmlTestCellRenderer} from './cell-types/test-cell-renderer';
import {ScatterPlotCellRenderer} from './sparklines/scatter-plot';

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

//name: Scatter Plot
//tags: cellRenderer
//output: grid_cell_renderer renderer
//meta.cellType: scatterplot
//meta.virtual: true
export function _ScatterPlotCellRenderer() {
  return new ScatterPlotCellRenderer();
}

