import * as DG from 'datagrok-api/dg';
import { GridCellRendererEx } from "../renderer/GridCellRendererEx";
export declare function isRowHeader(colGrid: DG.GridColumn): boolean;
export declare function getInstalledGridForColumn(colGrid: DG.GridColumn): DG.Grid | null;
export declare function installGridForColumn(grid: DG.Grid, colGrid: DG.GridColumn): boolean;
export declare function setGridColumnRenderer(colGrid: DG.GridColumn, renderer: GridCellRendererEx): void;
export declare function getGridColumnRenderer(colGrid: DG.GridColumn): GridCellRendererEx | null;
export declare function getGridColumnHeaderHeight(grid: DG.Grid): number;
export declare function getGridRowHeight(grid: DG.Grid): number;
export declare function getGridVisibleRowCount(grid: DG.Grid): number;
export declare function fillVisibleViewportRows(arMinMaxRowIdxs: Array<number>, grid: DG.Grid): void;
export declare function fillVisibleViewportGridCells(arColRowIdxs: Array<number>, grid: DG.Grid): void;
export declare function scaleFont(font: string, fFactor: number): string;
export declare function paintColHeaderCell(g: CanvasRenderingContext2D | null, nX: number, nY: number, nW: number, nH: number, colGrid: DG.GridColumn): void;
//# sourceMappingURL=GridUtils.d.ts.map