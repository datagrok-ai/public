import * as DG from 'datagrok-api/dg';
import {PinnedColumn} from "../pinned/PinnedColumn";

export class GridCellRendererEx extends DG.GridCellRenderer { // temporary to address a bug of importing during tests | extends DG.GridCellRenderer {
  onMouseDownEx(cellGrid : DG.GridCell, e : MouseEvent, nXOnCell : number, nYOnCell : number) : void {
    //this.onMouseDown(cellGrid, e);
  }

  onMouseUpEx(cellGrid : DG.GridCell, e : MouseEvent, nXOnCell : number, nYOnCell : number) : void {
    //this.onMouseUp(cellGrid, e);
  }

  onClickEx(cellGrid : DG.GridCell, e : MouseEvent, nXOnCell : number, nYOnCell : number) : void {
    //this.onClick(cellGrid, e);
  }

  onMouseMoveEx(cellGrid : DG.GridCell, e : MouseEvent, nXOnCell : number, nYOnCell : number) : void {
    //this.onMouseMove(cellGrid, e);
  }

  onMouseEnterEx(cellGrid : DG.GridCell, e : MouseEvent, nXOnCell : number, nYOnCell : number) : void {
    //this.onMouseEnter(cellGrid, e);
  }

  onMouseLeaveEx(cellGrid : DG.GridCell, e : MouseEvent, nXOnCell : number, nYOnCell : number) : void {
    //this.onMouseLeave(cellGrid, e);
  }

  onResizeWidth(colGrid : DG.GridColumn | PinnedColumn, grid : DG.Grid, nWCol : number, bAdjusting : boolean) : void  {

  }
  onResizeHeight(colGrid : DG.GridColumn | PinnedColumn, grid : DG.Grid, nHRow : number, bAdjusting : boolean) : void {

  }

  render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, value: any, context: any): void {

  }
}
