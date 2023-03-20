import * as DG from 'datagrok-api/dg';
import {PinnedColumn} from "../pinned/PinnedColumn";
import * as GridUtils from '../utils/GridUtils';
import {RendererUIManager} from '../renderer/RendererUIManager';
import {TooltipManager} from "../tooltip/TooltipManager";

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

  render(g: CanvasRenderingContext2D, nX: number, nY: number, nW: number, nH: number, cellGrid: DG.GridCell,  style: DG.GridCellStyle): void {
    const r = GridUtils.getGridColumnRenderer(cellGrid.gridColumn);
    if (r === null)
      GridUtils.setGridColumnRenderer(cellGrid.gridColumn, this);

    if (!RendererUIManager.isRegistered(cellGrid.grid))
      RendererUIManager.register(cellGrid.grid);

    if (!TooltipManager.isRegisted(cellGrid.grid))
      TooltipManager.register(cellGrid.grid)
  }

  tooltip(cell: DG.GridCell) : HTMLElement | null{
    return null;
  }
}
