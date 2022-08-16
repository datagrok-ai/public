import * as DG from 'datagrok-api/dg';
import { PinnedColumn } from "../pinned/PinnedColumn";
export declare class GridCellRendererEx {
    onMouseDownEx(cellGrid: DG.GridCell, e: MouseEvent, nXOnCell: number, nYOnCell: number): void;
    onMouseUpEx(cellGrid: DG.GridCell, e: MouseEvent, nXOnCell: number, nYOnCell: number): void;
    onClickEx(cellGrid: DG.GridCell, e: MouseEvent, nXOnCell: number, nYOnCell: number): void;
    onMouseMoveEx(cellGrid: DG.GridCell, e: MouseEvent, nXOnCell: number, nYOnCell: number): void;
    onMouseEnterEx(cellGrid: DG.GridCell, e: MouseEvent, nXOnCell: number, nYOnCell: number): void;
    onMouseLeaveEx(cellGrid: DG.GridCell, e: MouseEvent, nXOnCell: number, nYOnCell: number): void;
    onResizeWidth(colGrid: DG.GridColumn | PinnedColumn, grid: DG.Grid, nWCol: number, bAdjusting: boolean): void;
    onResizeHeight(colGrid: DG.GridColumn | PinnedColumn, grid: DG.Grid, nHRow: number, bAdjusting: boolean): void;
    render(g: CanvasRenderingContext2D, x: number, y: number, w: number, h: number, value: any, context: any): void;
}
//# sourceMappingURL=GridCellRendererEx.d.ts.map