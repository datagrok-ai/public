import * as DG from 'datagrok-api/dg';
export declare class PinnedColumn {
    private static MIN_ROW_HEIGHT;
    private static MAX_ROW_HEIGHT;
    private static SELECTION_COLOR;
    private static ACTIVE_CELL_COLOR;
    private static Y_RESIZE_SENSITIVITY;
    private m_fDevicePixelRatio;
    private m_colGrid;
    private m_root;
    private m_nWidthBug;
    private m_observerResizeGrid;
    private m_handlerVScroll;
    private m_handlerRowsFiltering;
    private m_handlerCurrRow;
    private m_handlerSel;
    private m_handlerRowsResized;
    private m_handlerRowsSorted;
    private m_handlerMouseDown;
    private m_handlerMouseUp;
    private m_handlerMouseLeave;
    private m_handlerMouseMove;
    private m_handlerMouseWheel;
    constructor(colGrid: DG.GridColumn);
    isPinned(): boolean;
    getGridColumn(): DG.GridColumn | null;
    getWidth(): number;
    getRoot(): HTMLCanvasElement | null;
    close(): void;
    private paint;
    private static hitTestRows;
}
//# sourceMappingURL=PinnedColumn.d.ts.map