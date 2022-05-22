import * as DG from 'datagrok-api/dg';

export class GridCellRendererEx extends DG.GridCellRenderer {
  isClickable(cellGrid : DG.GridCell, nXOnCell : number, nYOnCell : number) {
    return false;
  }

  onMouseDownEx(cellGrid : DG.GridCell, e : MouseEvent, nXOnCell : number, nYOnCell : number) : void {
    this.onMouseDown(cellGrid, e);
  }

  onMouseUpEx(cellGrid : DG.GridCell, e : MouseEvent, nXOnCell : number, nYOnCell : number) : void {
    this.onMouseUp(cellGrid, e);
  }

  onClickEx(cellGrid : DG.GridCell, e : MouseEvent, nXOnCell : number, nYOnCell : number) : void {
    this.onClick(cellGrid, e);
  }

  onMouseMoveEx(cellGrid : DG.GridCell, e : MouseEvent, nXOnCell : number, nYOnCell : number) : void {
    this.onMouseMove(cellGrid, e);
  }

  onMouseEnterEx(cellGrid : DG.GridCell, e : MouseEvent, nXOnCell : number, nYOnCell : number) : void {
    this.onMouseEnter(cellGrid, e);
  }

  onMouseLeaveEx(cellGrid : DG.GridCell, e : MouseEvent, nXOnCell : number, nYOnCell : number) : void {
    this.onMouseLeave(cellGrid, e);
  }
}
