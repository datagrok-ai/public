import * as DG from 'datagrok-api/dg';
import { GridCellRendererEx } from "./GridCellRendererEx";
import { PinnedColumn } from "../pinned/PinnedColumn";
export declare class LoadableImageRenderer extends GridCellRendererEx {
    constructor();
    valueToImageId(obValue: any): string | null;
    createImage(strImageId: string, nWImage: number, nHImage: number, fnImageRadyCallback: any): void;
    onResizeWidth(colGridOrPinned: DG.GridColumn | PinnedColumn, grid: DG.Grid, nWCol: number, bAdjusting: boolean): void;
    onResizeHeight(colGridOrPinned: DG.GridColumn | PinnedColumn, grid: DG.Grid, nHRow: number, bAdjusting: boolean): void;
    private onResize;
    render(g: CanvasRenderingContext2D, nX: number, nY: number, nW: number, nH: number, value: any, context: any): void;
    private m_mapImages;
    private m_bCellsAdjusting;
}
//# sourceMappingURL=LoadableImageRenderer.d.ts.map