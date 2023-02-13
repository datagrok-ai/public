import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';

import {MarkupNodeType} from './markup';
import {RectangleTreePlacer} from './rectangle-tree-placer';


export class GridTreePlacer<TNode extends MarkupNodeType> extends RectangleTreePlacer<TNode> {
  private readonly grid: DG.Grid;

  constructor(grid: DG.Grid, totalLength: number) {
    const top: number = Math.floor(grid.vertScroll.min) - 0.5;
    const rowsGridHeight: number = grid.root.clientHeight - grid.colHeaderHeight;
    const bottom = top + (rowsGridHeight / grid.props.rowHeight);

    super(top, bottom, totalLength);

    this.grid = grid;

    this.grid.onBeforeDrawContent.subscribe(this.gridOnChanged.bind(this));
    ui.onSizeChanged(this.grid.root).subscribe(this.gridOnChanged.bind(this));
  }

  // getNode(node: MarkupNodeType, point: DG.Point): RectangleTreeHoverType | null {
  //   // TODO: implement getNode
  //   return null;
  // }

  private gridOnChanged() {
    // const firstRowIndex: number = Math.floor(this.grid.vertScroll.min);
    // const rowsGridHeight: number = this.grid.root.clientHeight - this.grid.colHeaderHeight;
    // const lastRowIndex: number = firstRowIndex + Math.ceil(rowsGridHeight / this.grid.props.rowHeight);

    const top: number = Math.floor(this.grid.vertScroll.min) - 0.5;
    const rowsGridHeight: number = this.grid.root.clientHeight - this.grid.colHeaderHeight;
    const bottom = top + (rowsGridHeight / this.grid.props.rowHeight);

    this.update({top, bottom});
  }
}
