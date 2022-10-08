import * as DG from 'datagrok-api/dg';
import { PinnedColumn } from "./PinnedColumn";
export declare function getTotalPinnedRowsHeight(grid: DG.Grid): number;
export declare function getPinnedColumnLeft(colPinned: PinnedColumn): number;
export declare function getTotalPinnedColsWidth(grid: DG.Grid): number;
export declare function findPinnedColumnByRoot(eCanvas: HTMLCanvasElement, grid: DG.Grid): PinnedColumn | null;
export declare function getPinnedColumnCount(grid: DG.Grid): number;
export declare function getPinnedColumn(nIdx: number, grid: DG.Grid): PinnedColumn;
export declare function addPinnedColumn(colGrid: DG.GridColumn): PinnedColumn;
export declare function closeAllPinnedColumns(grid: DG.Grid): void;
export declare function installPinnedColumns(grid: DG.Grid): void;
export declare function isPinnedColumn(colGrid: DG.GridColumn): boolean;
export declare function isPinnableColumn(colGrid: DG.GridColumn): boolean;
export declare function registerPinnedColumns(): void;
export declare function handleContextMenu(args: any, fnMenuCallback: Function): void;
//# sourceMappingURL=PinnedUtils.d.ts.map