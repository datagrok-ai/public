/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ImageCellRenderer} from "./cell-types/image-cell-renderer";
import {HyperlinkCellRenderer} from "./cell-types/hyperlink-cell-renderer";

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