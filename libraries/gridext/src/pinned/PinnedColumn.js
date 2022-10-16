import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as GridUtils from '../utils/GridUtils';
import * as TextUtils from '../utils/TextUtils';
import { ColorUtils } from '../utils/ColorUtils';
import * as rxjs from 'rxjs';
import { GridCellRendererEx } from "../renderer/GridCellRendererEx";
import * as PinnedUtils from "./PinnedUtils";
import { MouseDispatcher } from "../ui/MouseDispatcher";
import { toDart } from "datagrok-api/dg";
//import {TableView} from "datagrok-api/dg";
/*
const hSubscriber  = grok.events.onViewLayoutApplied.subscribe((layout : DG.ViewLayout) => {
  const view : DG.TableView = layout.view as TableView;
  const itViewers = view.viewers;
  const arViewers = Array.from(itViewers);

  let viewer = null;
  const nViewerCount = arViewers.length;
  for (let n = 0; n < nViewerCount; ++n) {
    viewer = arViewers[n];
    if (viewer.type !== "Grid")
      continue;

    PinnedUtils.installPinnedColumns(viewer as DG.Grid);
  }
});
*/
function getRenderer(cell) {
    const colGrid = cell.gridColumn;
    if (colGrid === null || colGrid === undefined) {
        throw new Error('Grid cell is detached from the Grid column');
    }
    let renderer = GridUtils.getGridColumnRenderer(colGrid);
    if (renderer instanceof GridCellRendererEx) {
        return renderer;
    }
    return cell.renderer;
}
function getGrid(colGrid) {
    let grid = colGrid.grid;
    if (grid === null) {
        grid = GridUtils.getInstalledGridForColumn(colGrid);
        if (grid instanceof DG.Grid)
            return grid;
    }
    return grid;
}
function notifyAllColsRowsResized(grid, nHRows, bAdjusting) {
    let renderer = null;
    let colGrid = null;
    const lstColsGrid = grid.columns;
    const nColCount = lstColsGrid.length;
    for (let nCol = 0; nCol < nColCount; ++nCol) {
        colGrid = lstColsGrid.byIndex(nCol);
        if (colGrid === null || !colGrid.visible) {
            continue;
        }
        renderer = GridUtils.getGridColumnRenderer(colGrid);
        if (renderer instanceof GridCellRendererEx) {
            renderer.onResizeHeight(colGrid, grid, nHRows, bAdjusting);
        }
    }
}
function notifyAllPinnedColsRowsResized(colPinnedSource, nHRows, bAdjusting) {
    const colGridSource = colPinnedSource.getGridColumn();
    if (colGridSource === null) {
        return;
    }
    const grid = getGrid(colGridSource);
    const dart = DG.toDart(grid);
    if (dart.m_arPinnedCols === undefined) {
        throw new Error('Pinned Columns are not installed.');
    }
    let renderer = null;
    let colPinned = null;
    let colGrid = null;
    const nPinnedColCount = dart.m_arPinnedCols.length;
    for (let nColPin = 0; nColPin < nPinnedColCount; ++nColPin) {
        colPinned = dart.m_arPinnedCols[nColPin];
        colGrid = colPinned.m_colGrid;
        if (colGrid === null) {
            throw new Error('Pinned Column is detached.');
        }
        renderer = GridUtils.getGridColumnRenderer(colGrid);
        if (renderer instanceof GridCellRendererEx && colPinned.m_root !== null && grid !== null) {
            renderer.onResizeHeight(colPinned, grid, nHRows, bAdjusting);
        }
    }
}
const DEBUG = false;
export class PinnedColumn {
    constructor(colGrid) {
        var _a;
        this.m_nHResizeRowsBeforeDrag = -1;
        this.m_nResizeRowGridDragging = -1;
        this.m_nYResizeDraggingAnchor = -1;
        this.m_nResizeRowGridMoving = -1;
        this.m_nYDraggingAnchor = -1;
        this.m_nRowGridDragging = -1;
        this.m_nWheelCount = 0;
        this.m_arXYMouseOnCellDown = [-2, -2];
        this.m_arXYMouseOnCellUp = [-1, -1];
        this.m_bSortedAscending = null;
        this.m_cellCurrent = null;
        MouseDispatcher.create();
        const grid = getGrid(colGrid);
        if (grid === null) {
            throw new Error("Column '" + colGrid.name + "' is not attached to the grid.");
        }
        if (!PinnedUtils.isPinnableColumn(colGrid)) {
            throw new Error("Column '" + colGrid.name + "' cannot be pinned. It either pinned or HTML.");
        }
        //let nRowMin = grid.minVisibleRow;
        //let nRowMax = grid.maxVisibleRow;
        //let nColMin = grid.minVisibleColumn;
        //let nColMax = grid.maxVisibleColumn;
        //const it = grid.pinnedRows;
        //const ar = Array.from(it);
        this.m_fDevicePixelRatio = window.devicePixelRatio;
        const dart = DG.toDart(grid);
        if (dart.m_arPinnedCols === undefined)
            dart.m_arPinnedCols = [];
        if (dart.m_arPinnedCols.length === 0 && !GridUtils.isRowHeader(colGrid)) {
            const colGrid0 = grid.columns.byIndex(0);
            if (colGrid0 !== null && colGrid0 !== undefined)
                new PinnedColumn(colGrid0);
        }
        const nWTotalPinnedCols = PinnedUtils.getTotalPinnedColsWidth(grid);
        dart.m_arPinnedCols.push(this);
        const viewTable = grid.view;
        const dframe = grid.dataFrame;
        const nW = colGrid.width;
        this.m_colGrid = colGrid;
        this.m_nWidthBug = -1;
        try {
            colGrid.visible = false;
        }
        catch (e) {
            //DG bug
            console.error("ERROR: Couldn't hide column '" + colGrid.name + "' due to a DG bug. Attempt to set the width to 0");
            try {
                this.m_nWidthBug = colGrid.width;
                colGrid.width = 0;
            }
            catch (e) {
                //DG bug
                console.error("ERROR: Couldn't set the width to 0 for column '" + colGrid.name + "' due to a DG bug. This could be ignored if the column visually looks ok.");
            }
        }
        if (!GridUtils.isRowHeader(colGrid)) {
            if (colGrid.settings === null || colGrid.settings === undefined)
                colGrid.settings = {};
            colGrid.settings.isPinned = true; //this will be saved with the layout
            colGrid.settings.idxPinned = dart.m_arPinnedCols.length - 1;
        }
        grid.canvas.style.left = (grid.canvas.offsetLeft + nW).toString() + "px";
        grid.overlay.style.left = (grid.overlay.offsetLeft + nW).toString() + "px";
        grid.canvas.style.width = (grid.canvas.offsetWidth - nW).toString() + "px";
        grid.overlay.style.width = (grid.overlay.offsetWidth - nW).toString() + "px";
        const nHeight = grid.canvas.height; //canvas pixel height
        const eCanvasThis = ui.canvas(nW * window.devicePixelRatio, nHeight);
        const tabIndex = grid.canvas.getAttribute("tabIndex");
        if (tabIndex !== null)
            eCanvasThis.setAttribute("tabIndex", tabIndex);
        eCanvasThis.style.position = "absolute";
        eCanvasThis.style.left = nWTotalPinnedCols + "px";
        eCanvasThis.style.top = grid.canvas.offsetTop + "px";
        eCanvasThis.style.width = nW + "px";
        eCanvasThis.style.height = Math.round(nHeight / window.devicePixelRatio) + "px";
        //console.log("h " + grid.canvas.height + " offset " + grid.canvas.offsetHeight);
        if (grid.canvas.parentNode === null)
            throw new Error("Parent node for canvas cannot be null.");
        grid.canvas.parentNode.insertBefore(eCanvasThis, grid.canvas);
        this.m_root = eCanvasThis;
        const colGrid0 = grid.columns.byIndex(0);
        if (colGrid0 !== null && colGrid0 !== undefined) { //DG Bug from reading layout
            try {
                colGrid0.visible = false;
            }
            catch (e) {
                console.error("ERROR: Couldn't hide row header.");
            }
        }
        //OnResize Row header
        const headerThis = this; /*
        this.m_observerResize = new ResizeObserver(entries => {
          const g = headerThis.m_root.getContext('2d');
          for (let entry of entries) {
            headerThis.paint(g, grid);
          }
        });
        this.m_observerResize.observe(headerThis.m_root);*/
        //OnResize Grid
        this.m_observerResizeGrid = new ResizeObserver(function (entries) {
            const bCurrent = DG.toDart(grok.shell.v) === DG.toDart(viewTable);
            if (!bCurrent)
                return;
            if (headerThis.m_fDevicePixelRatio !== window.devicePixelRatio || grid.canvas.height !== eCanvasThis.height) {
                eCanvasThis.width = nW * window.devicePixelRatio;
                eCanvasThis.height = grid.canvas.height;
                eCanvasThis.style.top = grid.canvas.offsetTop + "px";
                eCanvasThis.style.width = nW + "px";
                eCanvasThis.style.height = Math.round(grid.canvas.height / window.devicePixelRatio) + "px";
                headerThis.m_fDevicePixelRatio = window.devicePixelRatio;
            }
            //console.log("Grid Resize: " + grid.canvas.height + " " + window.devicePixelRatio);
            //eCanvasThis.style.height = grid.root.style.height;
            /*
                  const eCanvasNew = ui.canvas(nW, grid.root.offsetHeight);
                  if(headerThis.m_root.parentNode !== null) {
                    headerThis.m_root.parentNode.replaceChild(eCanvasNew, headerThis.m_root);
                    headerThis.m_root = eCanvasNew;
                  }*/
            //headerThis.m_root.height = grid.root.offsetHeight;
            const g = eCanvasThis.getContext('2d');
            for (let entry of entries) {
                setTimeout(() => { headerThis.paint(g, grid); }, 100);
            }
        });
        (_a = this.m_observerResizeGrid) === null || _a === void 0 ? void 0 : _a.observe(grid.canvas);
        this.m_handlerKeyDown = rxjs.fromEvent(eCanvasThis, 'keydown').subscribe((e) => {
            //alert('up');
            setTimeout(() => {
                const ee = new KeyboardEvent(e.type, e);
                try {
                    grid.overlay.dispatchEvent(ee);
                }
                catch (ex) {
                    //console.error(ex.message);
                }
            }, 1);
        });
        const scrollVert = grid.vertScroll;
        this.m_handlerVScroll = scrollVert.onValuesChanged.subscribe(() => {
            const g = eCanvasThis.getContext('2d');
            headerThis.paint(g, grid);
        });
        this.m_handlerRowsFiltering = dframe.onRowsFiltering.subscribe(() => {
            setTimeout(() => {
                const g = eCanvasThis.getContext('2d');
                headerThis.paint(g, grid);
            }, 100);
        });
        this.m_handlerCurrRow = dframe.onCurrentRowChanged.subscribe(() => {
            const g = eCanvasThis.getContext('2d');
            headerThis.paint(g, grid);
        });
        this.m_handlerSel = dframe.onSelectionChanged.subscribe((e) => {
            const g = eCanvasThis.getContext('2d');
            headerThis.paint(g, grid);
        });
        this.m_handlerColsRemoved = dframe.onColumnsRemoved.subscribe((e) => {
            if (headerThis.m_colGrid === null)
                return;
            for (let nC = 0; nC < e.columns.length; ++nC) {
                if (e.columns[nC].name === headerThis.m_colGrid.name)
                    headerThis.close();
            }
        });
        this.m_handlerColNameChanged = dframe.onColumnNameChanged.subscribe((e) => {
            var _a;
            const dart = toDart(e);
            const strColNameOld = dart.newName;
            if (strColNameOld === ((_a = headerThis.m_colGrid) === null || _a === void 0 ? void 0 : _a.name)) {
                const g = eCanvasThis.getContext('2d');
                headerThis.paint(g, grid);
            }
        });
        /*
            this.m_handlerFilter = dframe.onRowsFiltered.subscribe((e : any) => {
                const g = eCanvasThis.getContext('2d');
                headerThis.paint(g, grid);
              }
            );
        */
        this.m_handlerRowsResized = grid.onRowsResized.subscribe((e) => {
            const g = eCanvasThis.getContext('2d');
            headerThis.paint(g, grid);
        });
        this.m_handlerRowsSorted = grid.onRowsSorted.subscribe((e) => {
            const g = eCanvasThis.getContext('2d');
            headerThis.paint(g, grid);
        });
    }
    isPinned() {
        return this.m_colGrid !== null;
    }
    getGridColumn() {
        return this.m_colGrid;
    }
    getWidth() {
        return this.m_root === null ? -1 : this.m_root.offsetWidth;
    }
    getRoot() {
        return this.m_root;
    }
    close() {
        var _a, _b, _c, _d, _e, _f, _g, _h, _j;
        if (this.m_colGrid === null) {
            throw new Error("Column has already been unpinned");
        }
        if (this.m_observerResizeGrid !== null) {
            this.m_observerResizeGrid.disconnect();
            this.m_observerResizeGrid = null;
        }
        /*my changes
            if(this.m_observerResize !== null) {
              this.m_observerResize.disconnect();
              this.m_observerResize = null;
            }
            */
        (_a = this.m_handlerKeyDown) === null || _a === void 0 ? void 0 : _a.unsubscribe();
        this.m_handlerKeyDown = null;
        (_b = this.m_handlerColsRemoved) === null || _b === void 0 ? void 0 : _b.unsubscribe();
        this.m_handlerColsRemoved = null;
        (_c = this.m_handlerColNameChanged) === null || _c === void 0 ? void 0 : _c.unsubscribe();
        this.m_handlerColNameChanged = null;
        (_d = this.m_handlerVScroll) === null || _d === void 0 ? void 0 : _d.unsubscribe();
        this.m_handlerVScroll = null;
        (_e = this.m_handlerRowsResized) === null || _e === void 0 ? void 0 : _e.unsubscribe();
        this.m_handlerRowsResized = null;
        (_f = this.m_handlerRowsSorted) === null || _f === void 0 ? void 0 : _f.unsubscribe();
        this.m_handlerRowsSorted = null;
        (_g = this.m_handlerRowsFiltering) === null || _g === void 0 ? void 0 : _g.unsubscribe();
        this.m_handlerRowsFiltering = null;
        (_h = this.m_handlerCurrRow) === null || _h === void 0 ? void 0 : _h.unsubscribe();
        this.m_handlerCurrRow = null;
        (_j = this.m_handlerSel) === null || _j === void 0 ? void 0 : _j.unsubscribe();
        this.m_handlerSel = null;
        const grid = getGrid(this.m_colGrid);
        if (grid === null) {
            throw new Error("Column '" + this.m_colGrid.name + "' is disconnected from grid.");
        }
        const dart = DG.toDart(grid);
        const ar = dart.m_arPinnedCols;
        const nIdx = ar.indexOf(this);
        ar.splice(nIdx, 1);
        if (this.m_root === null)
            throw new Error('Root cannot be null');
        let nIdxPinned = -1;
        let colGridTmp = null;
        for (let n = nIdx; n < ar.length; ++n) {
            colGridTmp = ar[n];
            colGridTmp.m_root.style.left = (colGridTmp.m_root.offsetLeft - this.m_root.offsetWidth).toString() + "px";
            nIdxPinned = colGridTmp.m_colGrid.settings.idxPinned;
            colGridTmp.m_colGrid.settings.idxPinned = n;
        }
        if (!GridUtils.isRowHeader(this.m_colGrid)) {
            this.m_colGrid.settings.idxPinned = -1;
            this.m_colGrid.settings.isPinned = false;
        }
        if (this.m_nWidthBug >= 0) {
            try {
                this.m_colGrid.width = this.m_nWidthBug;
            }
            catch (e) {
                //DG bug
                console.error("ERROR: Couldn't set the width to " + this.m_nWidthBug + " for column '" + this.m_colGrid.name + "' due to a DG bug. This could be ignored if the column visually looks ok.");
            }
        }
        try {
            this.m_colGrid.width = this.m_root.offsetWidth;
            this.m_colGrid.visible = true;
        }
        catch (e) {
            //DG bug
            console.error("ERROR: Couldn't show column '" + this.m_colGrid.name + "' due to a DG bug. This could be ignored if the column visually looks ok.");
        }
        grid.canvas.style.left = (grid.canvas.offsetLeft - this.m_root.offsetWidth).toString() + "px";
        grid.overlay.style.left = (grid.overlay.offsetLeft - this.m_root.offsetWidth).toString() + "px";
        grid.canvas.style.width = (grid.canvas.offsetWidth + this.m_root.offsetWidth).toString() + "px";
        grid.overlay.style.width = (grid.overlay.offsetWidth + this.m_root.offsetWidth).toString() + "px";
        if (this.m_root.parentNode !== null)
            this.m_root.parentNode.removeChild(this.m_root);
        this.m_root = null;
        if (dart.m_arPinnedCols.length === 1 && dart.m_arPinnedCols[0].m_colGrid.idx === 0 && this.m_colGrid.idx !== 0) {
            // try{colGrid0.visible = true;}
            try {
                dart.m_arPinnedCols[0].close();
            }
            catch (e) {
                console.error("ERROR: Couldn't close pinned column '" + dart.m_arPinnedCols[0].m_colGrid.name + "' ");
            }
        }
        this.m_colGrid = null;
    }
    onMouseEnter(e) {
        var _a;
        if (DEBUG)
            console.log('Mouse Enter Pinned Column: ' + ((_a = this.getGridColumn()) === null || _a === void 0 ? void 0 : _a.name));
    }
    onMouseMove(e) {
        var _a;
        if (DEBUG)
            console.log('Mouse Move Pinned Column: ' + ((_a = this.getGridColumn()) === null || _a === void 0 ? void 0 : _a.name));
        if (this.m_colGrid === null || this.m_root === null)
            return;
        const grid = this.m_colGrid.grid;
        const viewTable = grid.view;
        if (DG.toDart(grok.shell.v) !== DG.toDart(viewTable)) {
            return;
        }
        const arXYOnCell = [-1, -1];
        let nRowGrid = PinnedColumn.hitTestRows(this.m_root, grid, e, false, arXYOnCell);
        if (nRowGrid >= 0) {
            const cell = grid.cell(this.m_colGrid.name, nRowGrid);
            const renderer = getRenderer(cell);
            if (renderer instanceof GridCellRendererEx) {
                if (this.m_cellCurrent === null) {
                    renderer.onMouseEnterEx(cell, e, arXYOnCell[0], arXYOnCell[1]);
                }
                if (this.m_cellCurrent !== null && nRowGrid !== this.m_cellCurrent.gridRow) {
                    renderer.onMouseLeaveEx(this.m_cellCurrent, e, -1, -1);
                    renderer.onMouseEnterEx(cell, e, arXYOnCell[0], arXYOnCell[1]);
                }
                renderer.onMouseMoveEx(cell, e, arXYOnCell[0], arXYOnCell[1]);
            }
            this.m_cellCurrent = cell;
        }
        else if (this.m_cellCurrent !== null) {
            const renderer = getRenderer(this.m_cellCurrent);
            if (renderer instanceof GridCellRendererEx) {
                renderer.onMouseLeaveEx(this.m_cellCurrent, e, -1, -1);
            }
            this.m_cellCurrent = null;
        }
        nRowGrid = PinnedColumn.hitTestRows(this.m_root, grid, e, true, undefined);
        if (nRowGrid >= 0) {
            this.m_nResizeRowGridMoving = nRowGrid;
            document.body.style.cursor = "row-resize";
            return;
        }
        if (this.m_nResizeRowGridMoving >= 0) {
            this.m_nResizeRowGridMoving = -1;
            document.body.style.cursor = "auto";
        }
        //Hamburger Menu
        const colGrid = this.getGridColumn();
        if (colGrid === null || colGrid.name === '')
            return;
        const eDivHamb = GridUtils.getToolIconDiv(colGrid.grid);
        const nHColHeader = GridUtils.getGridColumnHeaderHeight(colGrid.grid);
        if (0 <= e.offsetY && e.offsetY < nHColHeader) {
            eDivHamb === null || eDivHamb === void 0 ? void 0 : eDivHamb.style.removeProperty('visibility');
            eDivHamb === null || eDivHamb === void 0 ? void 0 : eDivHamb.setAttribute('column_name', colGrid.name);
            //console.log('ToolsIcon for column ' + colGrid.name);
            // @ts-ignore
            eDivHamb === null || eDivHamb === void 0 ? void 0 : eDivHamb.style.left = (PinnedUtils.getPinnedColumnLeft(this) + this.getWidth() - 18) + 'px';
            // @ts-ignore
            eDivHamb === null || eDivHamb === void 0 ? void 0 : eDivHamb.style.top = (GridUtils.getGridColumnHeaderHeight(colGrid.grid) - 16) + "px";
        }
        else {
            const colGrid = this.getGridColumn();
            if (colGrid != null) {
                eDivHamb === null || eDivHamb === void 0 ? void 0 : eDivHamb.setAttribute('column_name', '');
                // @ts-ignore
                eDivHamb === null || eDivHamb === void 0 ? void 0 : eDivHamb.style.visibility = 'hidden';
            }
        }
    }
    onMouseDrag(e) {
        var _a;
        if (DEBUG)
            console.log('Mouse Drag Pinned Column: ' + ((_a = this.getGridColumn()) === null || _a === void 0 ? void 0 : _a.name));
        if (this.m_colGrid === null || this.m_root === null)
            return;
        const grid = this.m_colGrid.grid;
        const viewTable = grid.view;
        if (DG.toDart(grok.shell.v) !== DG.toDart(viewTable)) {
            return;
        }
        const bResizing = this.m_nResizeRowGridDragging >= 0;
        if (bResizing) {
            //console.log("Dragging : " + headerThis.m_strColName);
            const nYDiff = e.clientY - this.m_nYResizeDraggingAnchor;
            let nHRowGrid = this.m_nHResizeRowsBeforeDrag + nYDiff;
            if (nHRowGrid < PinnedColumn.MIN_ROW_HEIGHT)
                nHRowGrid = PinnedColumn.MIN_ROW_HEIGHT;
            else if (nHRowGrid > PinnedColumn.MAX_ROW_HEIGHT)
                nHRowGrid = PinnedColumn.MAX_ROW_HEIGHT;
            const eCanvasThis = this.m_root;
            let g = eCanvasThis.getContext('2d');
            if (g === null)
                return;
            g.fillStyle = "white";
            const nHHeaderCols = GridUtils.getGridColumnHeaderHeight(grid);
            g.fillRect(0, nHHeaderCols, eCanvasThis.offsetWidth, eCanvasThis.offsetHeight);
            grid.setOptions({
                rowHeight: nHRowGrid //this won't trigger onRowsRezized event, which is a DG bug
            });
            notifyAllPinnedColsRowsResized(this, nHRowGrid, true);
            notifyAllColsRowsResized(grid, nHRowGrid, true);
            let header = null;
            const ar = grid.dart.m_arPinnedCols;
            for (let n = 0; n < ar.length; ++n) {
                header = ar[n];
                g = header.m_root.getContext('2d');
                header.paint(g, grid);
            }
            try {
                const colGrid0 = grid.columns.byIndex(0);
                if (colGrid0 !== null)
                    colGrid0.visible = false; //temporary addressed the DG bug
            }
            catch (e) {
                //DG bug
            }
            return;
        }
    }
    onMouseLeave(e, bOverlap) {
        var _a;
        if (DEBUG)
            console.log('Mouse Left Pinned Column: ' + ((_a = this.getGridColumn()) === null || _a === void 0 ? void 0 : _a.name) + '  overlap: ' + bOverlap);
        if (this.m_nResizeRowGridMoving >= 0) {
            this.m_nResizeRowGridMoving = -1;
            document.body.style.cursor = "auto";
        }
        if (this.m_cellCurrent !== null) {
            const renderer = getRenderer(this.m_cellCurrent);
            if (renderer instanceof GridCellRendererEx) {
                const eMouse = e;
                renderer.onMouseLeaveEx(this.m_cellCurrent, eMouse, -1, -1);
            }
            this.m_cellCurrent = null;
        }
        const colGrid = this.getGridColumn();
        if (colGrid != null && !bOverlap) {
            const eDivHamb = GridUtils.getToolIconDiv(colGrid.grid);
            eDivHamb === null || eDivHamb === void 0 ? void 0 : eDivHamb.setAttribute('column_name', '');
            // @ts-ignore
            eDivHamb === null || eDivHamb === void 0 ? void 0 : eDivHamb.style.visibility = 'hidden';
        }
    }
    onMouseDblClick(e) {
        var _a, _b, _c;
        if (DEBUG)
            console.log('Mouse Dbl Clicked Pinned Column: ' + ((_a = this.getGridColumn()) === null || _a === void 0 ? void 0 : _a.name));
        if (this.m_colGrid === null || this.m_root === null)
            return;
        const grid = this.m_colGrid.grid;
        const viewTable = grid === null || grid === void 0 ? void 0 : grid.view;
        if (DG.toDart(grok.shell.v) !== DG.toDart(viewTable)) {
            return;
        }
        if (((_b = this.m_colGrid) === null || _b === void 0 ? void 0 : _b.name) === '')
            return;
        if (this.m_bSortedAscending == null)
            this.m_bSortedAscending = true;
        else if (this.m_bSortedAscending)
            this.m_bSortedAscending = false;
        else
            this.m_bSortedAscending = true;
        const nHHeaderCols = GridUtils.getGridColumnHeaderHeight(grid);
        if (0 <= e.offsetX && e.offsetX <= this.m_root.offsetWidth &&
            0 <= e.offsetY && e.offsetY <= nHHeaderCols) //on the rows header
         {
            grid === null || grid === void 0 ? void 0 : grid.sort([(_c = this.m_colGrid) === null || _c === void 0 ? void 0 : _c.name], [this.m_bSortedAscending]);
        }
    }
    onMouseDown(e) {
        var _a, _b;
        if (DEBUG)
            console.log('Mouse Down Pinned Column: ' + ((_a = this.getGridColumn()) === null || _a === void 0 ? void 0 : _a.name));
        /*
            if(e.view != null) {
              const ee = document.createEvent( "MouseEvent" );
              ee.initMouseEvent(e.type, e.bubbles, e.cancelable, e.view, e.detail, e.screenX + 100, e.screenY, e.clientX + 100, e.clientY, e.ctrlKey, e.altKey, e.shiftKey, e.metaKey, e.button, e.relatedTarget);
              this.m_colGrid?.grid.root.dispatchEvent(ee);
              return;
            }
        */
        if (this.m_colGrid === null)
            return;
        const grid = (_b = this.m_colGrid) === null || _b === void 0 ? void 0 : _b.grid;
        const viewTable = grid === null || grid === void 0 ? void 0 : grid.view;
        if (DG.toDart(grok.shell.v) !== DG.toDart(viewTable))
            return;
        if (e.buttons !== 1)
            return;
        let eCanvasThis = this.m_root;
        if (eCanvasThis === null)
            return;
        this.m_nResizeRowGridMoving = -1;
        const bAddToSel = e.ctrlKey || e.shiftKey;
        let nRowGrid = bAddToSel ? -1 : PinnedColumn.hitTestRows(eCanvasThis, grid, e, true, undefined);
        if (nRowGrid >= 0) {
            const nHRows = GridUtils.getGridRowHeight(grid);
            this.m_nResizeRowGridDragging = nRowGrid;
            this.m_nYResizeDraggingAnchor = e.clientY;
            this.m_nHResizeRowsBeforeDrag = nHRows;
        }
        else {
            nRowGrid = PinnedColumn.hitTestRows(eCanvasThis, grid, e, false, this.m_arXYMouseOnCellDown);
            this.m_nRowGridDragging = nRowGrid;
            this.m_nYDraggingAnchor = e.clientY;
            const cell = grid.cell(this.m_colGrid.name, nRowGrid);
            const renderer = getRenderer(cell);
            if (renderer instanceof GridCellRendererEx) {
                renderer.onMouseDownEx(cell, e, this.m_arXYMouseOnCellDown[0], this.m_arXYMouseOnCellDown[1]);
            }
        }
    }
    onMouseUp(e) {
        var _a, _b;
        if (DEBUG)
            console.log('Mouse Up Pinned Column: ' + ((_a = this.getGridColumn()) === null || _a === void 0 ? void 0 : _a.name));
        /*
            if(e.view != null) {
              const ee = document.createEvent( "MouseEvent" );
              ee.initMouseEvent(e.type, e.bubbles, e.cancelable, e.view, e.detail, e.screenX + 100, e.screenY, e.clientX + 100, e.clientY, e.ctrlKey, e.altKey, e.shiftKey, e.metaKey, e.button, e.relatedTarget);
              this.m_colGrid?.grid.root.dispatchEvent(ee);
              return;
            }
        */
        if (this.m_colGrid === null || this.m_root == null)
            return;
        const grid = (_b = this.m_colGrid) === null || _b === void 0 ? void 0 : _b.grid;
        const viewTable = grid === null || grid === void 0 ? void 0 : grid.view;
        if (DG.toDart(grok.shell.v) !== DG.toDart(viewTable)) {
            return;
        }
        if (this.m_nResizeRowGridDragging >= 0) {
            const nHRow = GridUtils.getGridRowHeight(grid);
            notifyAllPinnedColsRowsResized(this, nHRow, false);
            notifyAllColsRowsResized(grid, nHRow, false);
        }
        this.m_nHResizeRowsBeforeDrag = -1;
        this.m_nResizeRowGridDragging = -1;
        this.m_nYResizeDraggingAnchor = -1;
        this.m_nResizeRowGridMoving = -1;
        document.body.style.cursor = "auto";
        if (this.m_nRowGridDragging >= 0) {
            const dframe = grid.dataFrame;
            const bCtrl = e.ctrlKey;
            const bRangeSel = e.shiftKey;
            let bSel = true;
            const nRowGrid = PinnedColumn.hitTestRows(this.m_root, grid, e, false, this.m_arXYMouseOnCellUp);
            if (!bCtrl && !bRangeSel && nRowGrid === this.m_nRowGridDragging) { //click on the same row which will become active
                let cellRH = null;
                try {
                    cellRH = grid.cell("", nRowGrid);
                }
                catch (e) {
                    let colG = null;
                    const lstCols = grid.columns;
                    for (let nC = 1; nC < lstCols.length; ++nC) {
                        colG = lstCols.byIndex(nC);
                        cellRH = colG === null ? null : grid.cell(colG.name, nRowGrid);
                        if (cellRH !== null)
                            break;
                    }
                }
                if (cellRH !== null) {
                    const nRowTable = cellRH.tableRowIndex;
                    if (nRowTable !== null)
                        dframe.currentRow = nRowTable;
                }
            }
            else {
                const bitsetSel = dframe.selection;
                let nRowMin = this.m_nRowGridDragging < nRowGrid ? this.m_nRowGridDragging : nRowGrid;
                let nRowMax = this.m_nRowGridDragging > nRowGrid ? this.m_nRowGridDragging : nRowGrid;
                if (bCtrl) {
                    nRowMin = nRowGrid;
                    nRowMax = nRowGrid;
                    let bCurSel = bitsetSel.get(nRowGrid);
                    bSel = !bCurSel;
                }
                else if (bRangeSel) {
                    let nRowGridActive = GridUtils.getActiveGridRow(grid);
                    if (nRowGridActive === null)
                        nRowGridActive = 0;
                    if (nRowMin == nRowMax) {
                        bitsetSel.setAll(false, true);
                        nRowMin = nRowGridActive < nRowGrid ? nRowGridActive : nRowGrid;
                        nRowMax = nRowGridActive > nRowGrid ? nRowGridActive : nRowGrid;
                    }
                }
                //if(!bCtrl || bRangeSel)
                //bitsetSel.setAll(false, true);
                let cellRH = null;
                let nRowTable = -1;
                for (let nRow = nRowMin; nRow <= nRowMax; ++nRow) {
                    try {
                        cellRH = grid.cell("", nRow);
                    }
                    catch (e) {
                        let colG = null;
                        const lstCols = grid.columns;
                        for (let nC = 1; nC < lstCols.length; ++nC) {
                            colG = lstCols.byIndex(nC);
                            cellRH = colG === null ? null : grid.cell(colG.name, nRowGrid);
                            if (cellRH !== null)
                                break;
                        }
                    }
                    if (cellRH !== null && cellRH.tableRowIndex !== null) {
                        nRowTable = cellRH.tableRowIndex;
                        bitsetSel.set(nRowTable, bSel, true);
                    }
                }
            }
            const cell = grid.cell(this.m_colGrid.name, nRowGrid);
            const renderer = getRenderer(cell);
            if (renderer instanceof GridCellRendererEx) {
                renderer.onMouseUpEx(cell, e, this.m_arXYMouseOnCellUp[0], this.m_arXYMouseOnCellUp[1]);
            }
            if (this.m_arXYMouseOnCellUp[0] === this.m_arXYMouseOnCellDown[0] && this.m_arXYMouseOnCellDown[1] === this.m_arXYMouseOnCellUp[1]) {
                if (renderer instanceof GridCellRendererEx) {
                    renderer.onClickEx(cell, e, this.m_arXYMouseOnCellUp[0], this.m_arXYMouseOnCellUp[1]);
                }
            }
            this.m_nRowGridDragging = -1;
            this.m_nYDraggingAnchor = -1;
            this.m_arXYMouseOnCellDown[0] = -2;
            this.m_arXYMouseOnCellDown[1] = -2;
            this.m_arXYMouseOnCellUp[0] = -1;
            this.m_arXYMouseOnCellUp[1] = -1;
        }
    }
    onContextMenu(e) {
        var _a;
        if (DEBUG)
            console.log('Context menu Pinned Column: ' + ((_a = this.getGridColumn()) === null || _a === void 0 ? void 0 : _a.name));
    }
    onMouseWheel(e) {
        var _a;
        if (this.m_colGrid === null || this.m_root == null)
            return;
        const grid = (_a = this.m_colGrid) === null || _a === void 0 ? void 0 : _a.grid;
        const viewTable = grid === null || grid === void 0 ? void 0 : grid.view;
        if (DG.toDart(grok.shell.v) !== DG.toDart(viewTable)) {
            return;
        }
        if (e.deltaX !== 0 || e.deltaZ !== 0) {
            return;
        }
        setTimeout(() => {
            const ee = new WheelEvent(e.type, e);
            try {
                grid.overlay.dispatchEvent(ee);
            }
            catch (ex) {
                //console.error(ex.message);
            }
        }, 1);
        if (true)
            return;
        //e.clientX = 5;
        if (this.m_nWheelCount === 1) {
            //scroll +
            const nRowCount = GridUtils.getGridVisibleRowCount(grid);
            const scrollY = grid.vertScroll;
            if (nRowCount - 1 > scrollY.max) {
                scrollY.setValues(scrollY.minRange, scrollY.maxRange, scrollY.min + 1, scrollY.max + 1);
            }
            this.m_nWheelCount = 0;
        }
        else if (this.m_nWheelCount === -1) {
            //scroll -
            const scrollY = grid.vertScroll;
            if (scrollY.min >= 1) {
                scrollY.setValues(scrollY.minRange, scrollY.maxRange, scrollY.min - 1, scrollY.max - 1);
            }
            this.m_nWheelCount = 0;
        }
        else {
            this.m_nWheelCount = e.deltaY > 0 ? 1 : -1;
        }
    }
    paint(g, grid) {
        //const nWDiv = entry.contentBoxSize ? entry.contentBoxSize[0].inlineSize : entry.contentRect.width;
        if (g === null) {
            return;
        }
        if (this.m_root === null) {
            throw new Error('Root cannot be null.');
        }
        if (this.m_colGrid === null) {
            throw new Error('Column grid cannot be null.');
        }
        const dframe = grid.dataFrame;
        const nW = this.m_root.offsetWidth;
        const nH = this.m_root.offsetHeight;
        g.fillStyle = "white";
        g.fillRect(0, 0, nW * window.devicePixelRatio, nH * window.devicePixelRatio);
        if (this.m_colGrid.name === null)
            return;
        const bitsetFilter = dframe.filter;
        if (bitsetFilter.falseCount === dframe.rowCount)
            return; //everything is filtered
        //column Header
        const options = grid.getOptions(true);
        const fontCellDefault = options.look.defaultCellFont;
        let font = options.look.colHeaderFont == null || options.look.colHeaderFont === undefined ? "bold 14px Volta Text, Arial" : options.look.colHeaderFont;
        let fontScaled = GridUtils.scaleFont(font, window.devicePixelRatio);
        g.font = fontScaled;
        let str = TextUtils.trimText(this.m_colGrid.name, g, nW);
        const tm = g.measureText(str);
        const nWLabel = tm.width;
        const nAscent = Math.abs(tm.actualBoundingBoxAscent);
        const nDescent = tm.actualBoundingBoxDescent;
        const nHFont = nAscent + nDescent; // + 2*nYInset;
        //let cellCH = grid.cell(this.m_colGrid.name, -1);
        //let renderer = cellCH.renderer;
        let nX = 0;
        let nY = 0;
        const nHCH = GridUtils.getGridColumnHeaderHeight(grid) * window.devicePixelRatio;
        g.textAlign = 'start';
        g.fillStyle = "Black";
        let nYOffset = Math.floor((nHCH - nHFont) / 2);
        const nXX = nX + ((nW * window.devicePixelRatio - nWLabel) >> 1);
        let nYY = (nY + nHCH - Math.ceil(3 * window.devicePixelRatio)); //-2*window.devicePixelRatio);
        //onsole.log("nXX " + nXX + " nYY = " + nYY + " CHH " + nHCH);
        g.fillText(str, nXX, nYY);
        //Regular cells
        const nRowCurrent = dframe.currentRow.idx;
        const bitsetSel = dframe.selection;
        const arRowsMinMax = [-1, -1];
        GridUtils.fillVisibleViewportRows(arRowsMinMax, grid);
        const nRowMin = arRowsMinMax[0];
        const nRowMax = arRowsMinMax[1];
        //console.log(nRowMin + " " + nRowMax);
        const nHRow = GridUtils.getGridRowHeight(grid);
        nYOffset = nHCH;
        const nHRowGrid = nHRow * window.devicePixelRatio;
        let cellRH = null;
        let nWW = nW * window.devicePixelRatio;
        //const nHH = nHRowGrid;
        //my changes const nPinnedRowCount = Array.from(grid.pinnedRows).length;
        const arTableRows = new Array(nRowMax - nRowMin + 1);
        let nRowTable = -1;
        let bSel = false;
        for (let nRG = nRowMin; nRG <= nRowMax; ++nRG) {
            try {
                cellRH = grid.cell(this.m_colGrid.name, nRG);
            }
            catch (e) //to address DG bug when everything is filtered
             {
                continue;
            }
            if (cellRH.tableRowIndex === undefined) //DG bug
                continue;
            //my changes if(this.m_colGrid.name == '')
            //my changescellRH.customText = nRG - nRowMin < nPinnedRowCount ? '' : (nRG - nPinnedRowCount +1).toString();
            nRowTable = cellRH.tableRowIndex === null ? -1 : cellRH.tableRowIndex;
            arTableRows[nRG - nRowMin] = nRowTable;
            nYY = nYOffset + (nRG - nRowMin) * nHRowGrid;
            let renderer = GridUtils.getGridColumnRenderer(cellRH.gridColumn);
            if (renderer === null) {
                try {
                    renderer = cellRH.renderer;
                }
                catch (e) {
                    console.error("Could not obtain renderer for DG cell. DG bug " + this.m_colGrid.name + " row " + nRG);
                    continue;
                }
            }
            if (renderer === null || renderer === undefined) {
                console.error("Couldn't find renderer for pinned column " + this.m_colGrid.name + " row " + nRG);
                continue;
            }
            //let nYY = nY;//*window.devicePixelRatio;
            font = cellRH.style.font;
            fontScaled = GridUtils.scaleFont(font, window.devicePixelRatio);
            if (fontScaled !== null) {
                cellRH.style.font = fontScaled;
            }
            if (nW > 0 && nHRowGrid > 0) { //to address a bug caused either DG or client app
                try {
                    if (renderer.name === 'Molecule') {
                        renderer.render(g, 0, nYY / window.devicePixelRatio, nWW / window.devicePixelRatio, nHRowGrid / window.devicePixelRatio, cellRH, cellRH.style);
                    }
                    else
                        renderer.render(g, 0, nYY, nWW, nHRowGrid, cellRH, cellRH.style);
                }
                catch (e) {
                    console.error("Could not paint cell for pinned column " + this.m_colGrid.name + " row " + nRG);
                    continue;
                    //throw e;
                }
            }
        }
        //Paint Grid
        g.strokeStyle = "Gainsboro";
        g.beginPath();
        g.moveTo(0, nY * window.devicePixelRatio);
        g.lineTo(0, (nY + nHCH - 1 * window.devicePixelRatio));
        g.stroke();
        g.beginPath();
        g.moveTo(0, nYOffset + 1);
        g.lineTo(nWW, nYOffset + 1);
        g.stroke();
        const nPinnedColCount = PinnedUtils.getPinnedColumnCount(grid);
        const colPinned = PinnedUtils.getPinnedColumn(nPinnedColCount - 1, grid);
        const bLast = this === colPinned;
        for (let nRG = nRowMin; nRG <= nRowMax; ++nRG) {
            nYY = nYOffset + (nRG - nRowMin) * nHRowGrid;
            //if(options.look.showRowGridlines) {
            g.strokeStyle = "Gainsboro";
            g.beginPath();
            g.moveTo(0, nYY + nHRowGrid + 1);
            g.lineTo(nWW, nYY + nHRowGrid + 1);
            g.stroke();
            g.beginPath();
            g.moveTo(0, nYY);
            g.lineTo(0, nYY + nHRowGrid + 1);
            g.stroke();
            if (bLast) { //my changes  && (nRG - nRowMin) >= nPinnedRowCount) {
                g.strokeStyle = "black";
                g.beginPath();
                g.moveTo(nWW - 1, nYY);
                g.lineTo(nWW - 1, nYY + nHRowGrid + 1);
                g.stroke();
            }
            //}
            nRowTable = arTableRows[nRG - nRowMin];
            try {
                bSel = nRowTable === undefined || nRowTable < 0 ? false : bitsetSel.get(nRowTable);
            }
            catch (e) {
                console.error('PaintError: row_min: ' + nRowMin + ' row_max: ' + nRowMax + ' nR ' + nRG + ' ' + nRowTable);
                throw e;
            }
            if (bSel) {
                g.globalAlpha = 0.2;
                g.fillStyle = PinnedColumn.SELECTION_COLOR;
                g.fillRect(0, nYY, nWW, nHRowGrid);
                g.globalAlpha = 1;
            }
            if (nRowCurrent === nRowTable) {
                g.globalAlpha = 0.2;
                g.fillStyle = PinnedColumn.ACTIVE_CELL_COLOR;
                g.fillRect(0, nYY, nWW, nHRowGrid);
                g.globalAlpha = 1;
            }
        } //for
    }
    static hitTestRows(eCanvasPinned, grid, e, bBorder, arXYOnCell) {
        const rect = eCanvasPinned.getBoundingClientRect();
        const scrollLeft = window.pageXOffset || document.documentElement.scrollLeft;
        const scrollTop = window.pageYOffset || document.documentElement.scrollTop;
        const nY = rect.top + scrollTop;
        const nX = rect.left + scrollLeft;
        if (nX <= e.clientX && e.clientX <= nX + eCanvasPinned.offsetWidth) //on the rows header
         {
            const nHHeaderCols = GridUtils.getGridColumnHeaderHeight(grid);
            const nHRowGrid = GridUtils.getGridRowHeight(grid);
            const arMinMaxRows = [-1, -1];
            GridUtils.fillVisibleViewportRows(arMinMaxRows, grid);
            const nRowMin = arMinMaxRows[0];
            const nRowMax = arMinMaxRows[1];
            const nYMouseOnHeader = e.clientY - nY;
            let nYBorder = -1;
            let nYDiff = -1;
            for (let nRow = nRowMin; nRow <= nRowMax; ++nRow) {
                nYBorder = nHHeaderCols + (nRow - nRowMin + 1) * nHRowGrid;
                nYDiff = nYMouseOnHeader - nYBorder;
                if (bBorder && Math.abs(nYDiff) <= PinnedColumn.Y_RESIZE_SENSITIVITY) {
                    return nRow;
                }
                if (!bBorder && nYBorder - nHRowGrid <= nYMouseOnHeader && nYMouseOnHeader <= nYBorder) {
                    if (arXYOnCell !== undefined) {
                        arXYOnCell[0] = e.clientX - nX;
                        arXYOnCell[1] = nYMouseOnHeader - nYBorder + nHRowGrid;
                    }
                    return nRow;
                }
            }
        }
        return -1;
    }
}
PinnedColumn.MIN_ROW_HEIGHT = 20;
PinnedColumn.MAX_ROW_HEIGHT = 500;
PinnedColumn.SELECTION_COLOR = ColorUtils.toRgb(ColorUtils.colSelection); //"rgba(237, 220, 88, 0.15)";
PinnedColumn.ACTIVE_CELL_COLOR = ColorUtils.toRgb(ColorUtils.currentRow); //"rgba(153, 237, 82, 0.25)";
PinnedColumn.Y_RESIZE_SENSITIVITY = 2;
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJmaWxlIjoiUGlubmVkQ29sdW1uLmpzIiwic291cmNlUm9vdCI6IiIsInNvdXJjZXMiOlsiUGlubmVkQ29sdW1uLnRzIl0sIm5hbWVzIjpbXSwibWFwcGluZ3MiOiJBQUFBLE9BQU8sS0FBSyxJQUFJLE1BQU0sbUJBQW1CLENBQUM7QUFDMUMsT0FBTyxLQUFLLEVBQUUsTUFBTSxpQkFBaUIsQ0FBQztBQUN0QyxPQUFPLEtBQUssRUFBRSxNQUFNLGlCQUFpQixDQUFDO0FBQ3RDLE9BQU8sS0FBSyxTQUFTLE1BQU0sb0JBQW9CLENBQUM7QUFDaEQsT0FBTyxLQUFLLFNBQVMsTUFBTSxvQkFBb0IsQ0FBQztBQUNoRCxPQUFPLEVBQUMsVUFBVSxFQUFDLE1BQU0scUJBQXFCLENBQUM7QUFDL0MsT0FBTyxLQUFLLElBQUksTUFBTSxNQUFNLENBQUM7QUFDN0IsT0FBTyxFQUFFLGtCQUFrQixFQUFDLE1BQU0sZ0NBQWdDLENBQUM7QUFDbkUsT0FBTyxLQUFLLFdBQVcsTUFBTSxlQUFlLENBQUM7QUFFN0MsT0FBTyxFQUFDLGVBQWUsRUFBQyxNQUFNLHVCQUF1QixDQUFDO0FBQ3RELE9BQU8sRUFBYyxNQUFNLEVBQUMsTUFBTSxpQkFBaUIsQ0FBQztBQUNwRCw0Q0FBNEM7QUFHNUM7Ozs7Ozs7Ozs7Ozs7Ozs7RUFnQkU7QUFFRixTQUFTLFdBQVcsQ0FBQyxJQUFrQjtJQUNyQyxNQUFNLE9BQU8sR0FBRyxJQUFJLENBQUMsVUFBVSxDQUFDO0lBQ2hDLElBQUksT0FBTyxLQUFLLElBQUksSUFBSSxPQUFPLEtBQUssU0FBUyxFQUFFO1FBQzdDLE1BQU0sSUFBSSxLQUFLLENBQUMsNENBQTRDLENBQUMsQ0FBQztLQUMvRDtJQUVELElBQUksUUFBUSxHQUFHLFNBQVMsQ0FBQyxxQkFBcUIsQ0FBQyxPQUFPLENBQUMsQ0FBQztJQUN4RCxJQUFHLFFBQVEsWUFBWSxrQkFBa0IsRUFBRTtRQUN6QyxPQUFPLFFBQVEsQ0FBQztLQUNqQjtJQUVELE9BQU8sSUFBSSxDQUFDLFFBQVEsQ0FBQztBQUN2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QkFFRCxJQUFHLENBQUMsT0FBTyxJQUFJLFFBQVEsR0FBRyxTQUFTLElBQUksZUFBZSxJQUFJLGVBQWUsSUFBSSxRQUFRLEVBQUU7b0JBRXJGLElBQUcsVUFBVSxLQUFLLFNBQVMsRUFBRTt3QkFDM0IsVUFBVSxDQUFDLENBQUMsQ0FBQyxHQUFHLENBQUMsQ0FBQyxPQUFPLEdBQUcsRUFBRSxDQUFDO3dCQUMvQixVQUFVLENBQUMsQ0FBQyxDQUFDLEdBQUcsZUFBZSxHQUFHLFFBQVEsR0FBRyxTQUFTLENBQUM7cUJBQ3hEO29CQUVELE9BQU8sSUFBSSxDQUFDO2lCQUNiO2FBQ0Y7U0FDRjtRQUVELE9BQU8sQ0FBQyxDQUFDLENBQUM7SUFDWixDQUFDOztBQTNtQ2MsMkJBQWMsR0FBRyxFQUFFLENBQUM7QUFDcEIsMkJBQWMsR0FBRyxHQUFHLENBQUM7QUFDckIsNEJBQWUsR0FBRyxVQUFVLENBQUMsS0FBSyxDQUFDLFVBQVUsQ0FBQyxZQUFZLENBQUMsQ0FBQyxDQUFDLDZCQUE2QjtBQUMxRiw4QkFBaUIsR0FBRyxVQUFVLENBQUMsS0FBSyxDQUFDLFVBQVUsQ0FBQyxVQUFVLENBQUMsQ0FBQyxDQUFDLDZCQUE2QjtBQUMxRixpQ0FBb0IsR0FBRyxDQUFDLENBQUMiLCJzb3VyY2VzQ29udGVudCI6WyJpbXBvcnQgKiBhcyBncm9rIGZyb20gJ2RhdGFncm9rLWFwaS9ncm9rJztcclxuaW1wb3J0ICogYXMgREcgZnJvbSAnZGF0YWdyb2stYXBpL2RnJztcclxuaW1wb3J0ICogYXMgdWkgZnJvbSAnZGF0YWdyb2stYXBpL3VpJztcclxuaW1wb3J0ICogYXMgR3JpZFV0aWxzIGZyb20gJy4uL3V0aWxzL0dyaWRVdGlscyc7XHJcbmltcG9ydCAqIGFzIFRleHRVdGlscyBmcm9tICcuLi91dGlscy9UZXh0VXRpbHMnO1xyXG5pbXBvcnQge0NvbG9yVXRpbHN9IGZyb20gJy4uL3V0aWxzL0NvbG9yVXRpbHMnO1xyXG5pbXBvcnQgKiBhcyByeGpzIGZyb20gJ3J4anMnO1xyXG5pbXBvcnQgeyBHcmlkQ2VsbFJlbmRlcmVyRXh9IGZyb20gXCIuLi9yZW5kZXJlci9HcmlkQ2VsbFJlbmRlcmVyRXhcIjtcclxuaW1wb3J0ICogYXMgUGlubmVkVXRpbHMgZnJvbSBcIi4vUGlubmVkVXRpbHNcIjtcclxuaW1wb3J0IHtnZXRHcmlkRGFydFBvcHVwTWVudSwgaXNIaXRUZXN0T25FbGVtZW50fSBmcm9tIFwiLi4vdXRpbHMvR3JpZFV0aWxzXCI7XHJcbmltcG9ydCB7TW91c2VEaXNwYXRjaGVyfSBmcm9tIFwiLi4vdWkvTW91c2VEaXNwYXRjaGVyXCI7XHJcbmltcG9ydCB7Q29sdW1uc0FyZ3MsIHRvRGFydH0gZnJvbSBcImRhdGFncm9rLWFwaS9kZ1wiO1xyXG4vL2ltcG9ydCB7VGFibGVWaWV3fSBmcm9tIFwiZGF0YWdyb2stYXBpL2RnXCI7XHJcblxyXG5cclxuLypcclxuY29uc3QgaFN1YnNjcmliZXIgID0gZ3Jvay5ldmVudHMub25WaWV3TGF5b3V0QXBwbGllZC5zdWJzY3JpYmUoKGxheW91dCA6IERHLlZpZXdMYXlvdXQpID0+IHtcclxuICBjb25zdCB2aWV3IDogREcuVGFibGVWaWV3ID0gbGF5b3V0LnZpZXcgYXMgVGFibGVWaWV3O1xyXG4gIGNvbnN0IGl0Vmlld2VycyA9IHZpZXcudmlld2VycztcclxuICBjb25zdCBhclZpZXdlcnMgPSBBcnJheS5mcm9tKGl0Vmlld2Vycyk7XHJcblxyXG4gIGxldCB2aWV3ZXIgPSBudWxsO1xyXG4gIGNvbnN0IG5WaWV3ZXJDb3VudCA9IGFyVmlld2Vycy5sZW5ndGg7XHJcbiAgZm9yIChsZXQgbiA9IDA7IG4gPCBuVmlld2VyQ291bnQ7ICsrbikge1xyXG4gICAgdmlld2VyID0gYXJWaWV3ZXJzW25dO1xyXG4gICAgaWYgKHZpZXdlci50eXBlICE9PSBcIkdyaWRcIilcclxuICAgICAgY29udGludWU7XHJcblxyXG4gICAgUGlubmVkVXRpbHMuaW5zdGFsbFBpbm5lZENvbHVtbnModmlld2VyIGFzIERHLkdyaWQpO1xyXG4gIH1cclxufSk7XHJcbiovXHJcblxyXG5mdW5jdGlvbiBnZXRSZW5kZXJlcihjZWxsIDogREcuR3JpZENlbGwpIDogR3JpZENlbGxSZW5kZXJlckV4IHwgREcuR3JpZENlbGxSZW5kZXJlciB7XHJcbiAgY29uc3QgY29sR3JpZCA9IGNlbGwuZ3JpZENvbHVtbjtcclxuICBpZiAoY29sR3JpZCA9PT0gbnVsbCB8fCBjb2xHcmlkID09PSB1bmRlZmluZWQpIHtcclxuICAgIHRocm93IG5ldyBFcnJvcignR3JpZCBjZWxsIGlzIGRldGFjaGVkIGZyb20gdGhlIEdyaWQgY29sdW1uJyk7XHJcbiAgfVxyXG5cclxuICBsZXQgcmVuZGVyZXIgPSBHcmlkVXRpbHMuZ2V0R3JpZENvbHVtblJlbmRlcmVyKGNvbEdyaWQpO1xyXG4gIGlmKHJlbmRlcmVyIGluc3RhbmNlb2YgR3JpZENlbGxSZW5kZXJlckV4KSB7XHJcbiAgICByZXR1cm4gcmVuZGVyZXI7XHJcbiAgfVxyXG5cclxuICByZXR1cm4gY2VsbC5yZW5kZXJlcjtcclxufVxyXG5cclxuXHJcbmZ1bmN0aW9uIGdldEdyaWQoY29sR3JpZCA6IERHLkdyaWRDb2x1bW4pIDogREcuR3JpZCB8IG51bGwge1xyXG4gIGxldCBncmlkIDogREcuR3JpZCB8IG51bGwgPSBjb2xHcmlkLmdyaWQ7XHJcbiAgaWYoIGdyaWQgPT09IG51bGwpIHtcclxuICAgIGdyaWQgPSBHcmlkVXRpbHMuZ2V0SW5zdGFsbGVkR3JpZEZvckNvbHVtbihjb2xHcmlkKTtcclxuICAgIGlmKGdyaWQgaW5zdGFuY2VvZiBERy5HcmlkKVxyXG4gICAgICByZXR1cm4gZ3JpZDtcclxuICB9XHJcblxyXG4gIHJldHVybiBncmlkO1xyXG59XHJcblxyXG5cclxuZnVuY3Rpb24gbm90aWZ5QWxsQ29sc1Jvd3NSZXNpemVkKGdyaWQgOiBERy5HcmlkLCBuSFJvd3MgOiBudW1iZXIsIGJBZGp1c3RpbmcgOiBib29sZWFuKSA6IHZvaWQge1xyXG5cclxuICBsZXQgcmVuZGVyZXIgOiBHcmlkQ2VsbFJlbmRlcmVyRXggfCBudWxsID0gbnVsbFxyXG4gIGxldCBjb2xHcmlkID0gbnVsbDtcclxuICBjb25zdCBsc3RDb2xzR3JpZCA9IGdyaWQuY29sdW1ucztcclxuICBjb25zdCBuQ29sQ291bnQgPSBsc3RDb2xzR3JpZC5sZW5ndGg7XHJcbiAgZm9yKGxldCBuQ29sPTA7IG5Db2w8bkNvbENvdW50OyArK25Db2wpIHtcclxuICAgIGNvbEdyaWQgPSBsc3RDb2xzR3JpZC5ieUluZGV4KG5Db2wpO1xyXG4gICAgaWYoY29sR3JpZCA9PT0gbnVsbCB8fCAhY29sR3JpZC52aXNpYmxlKXtcclxuICAgICAgY29udGludWVcclxuICAgIH1cclxuXHJcbiAgICByZW5kZXJlciA9IEdyaWRVdGlscy5nZXRHcmlkQ29sdW1uUmVuZGVyZXIoY29sR3JpZCk7XHJcbiAgICBpZiAocmVuZGVyZXIgaW5zdGFuY2VvZiBHcmlkQ2VsbFJlbmRlcmVyRXgpIHtcclxuICAgICAgcmVuZGVyZXIub25SZXNpemVIZWlnaHQoY29sR3JpZCwgZ3JpZCwgbkhSb3dzLCBiQWRqdXN0aW5nKTtcclxuICAgIH1cclxuICB9XHJcbn1cclxuXHJcblxyXG5mdW5jdGlvbiBub3RpZnlBbGxQaW5uZWRDb2xzUm93c1Jlc2l6ZWQoY29sUGlubmVkU291cmNlIDogUGlubmVkQ29sdW1uLCBuSFJvd3MgOiBudW1iZXIsIGJBZGp1c3RpbmcgOiBib29sZWFuKSA6IHZvaWQge1xyXG5cclxuICBjb25zdCBjb2xHcmlkU291cmNlICA9IGNvbFBpbm5lZFNvdXJjZS5nZXRHcmlkQ29sdW1uKCk7XHJcbiAgaWYoY29sR3JpZFNvdXJjZSA9PT0gbnVsbCl7XHJcbiAgICByZXR1cm47XHJcbiAgfVxyXG5cclxuICBjb25zdCBncmlkID0gZ2V0R3JpZChjb2xHcmlkU291cmNlKTtcclxuICBjb25zdCBkYXJ0ID0gREcudG9EYXJ0KGdyaWQpO1xyXG4gIGlmKGRhcnQubV9hclBpbm5lZENvbHMgPT09IHVuZGVmaW5lZCkge1xyXG4gICAgdGhyb3cgbmV3IEVycm9yKCdQaW5uZWQgQ29sdW1ucyBhcmUgbm90IGluc3RhbGxlZC4nKTtcclxuICB9XHJcblxyXG4gIGxldCByZW5kZXJlciA6IEdyaWRDZWxsUmVuZGVyZXJFeCB8IG51bGwgPSBudWxsXHJcbiAgbGV0IGNvbFBpbm5lZCA9IG51bGw7XHJcbiAgbGV0IGNvbEdyaWQgPSBudWxsO1xyXG4gIGNvbnN0IG5QaW5uZWRDb2xDb3VudCA9IGRhcnQubV9hclBpbm5lZENvbHMubGVuZ3RoO1xyXG4gIGZvcihsZXQgbkNvbFBpbj0wOyBuQ29sUGluPG5QaW5uZWRDb2xDb3VudDsgKytuQ29sUGluKSB7XHJcbiAgICBjb2xQaW5uZWQgPSBkYXJ0Lm1fYXJQaW5uZWRDb2xzW25Db2xQaW5dO1xyXG4gICAgY29sR3JpZCA9IGNvbFBpbm5lZC5tX2NvbEdyaWQ7XHJcbiAgICBpZihjb2xHcmlkID09PSBudWxsKSB7XHJcbiAgICAgIHRocm93IG5ldyBFcnJvcignUGlubmVkIENvbHVtbiBpcyBkZXRhY2hlZC4nKTtcclxuICAgIH1cclxuXHJcbiAgICByZW5kZXJlciA9IEdyaWRVdGlscy5nZXRHcmlkQ29sdW1uUmVuZGVyZXIoY29sR3JpZCk7XHJcbiAgICBpZiAocmVuZGVyZXIgaW5zdGFuY2VvZiBHcmlkQ2VsbFJlbmRlcmVyRXggICYmIGNvbFBpbm5lZC5tX3Jvb3QgIT09IG51bGwgJiYgZ3JpZCAhPT0gbnVsbCkge1xyXG4gICAgICByZW5kZXJlci5vblJlc2l6ZUhlaWdodChjb2xQaW5uZWQsIGdyaWQsIG5IUm93cywgYkFkanVzdGluZyk7XHJcbiAgICB9XHJcbiAgfVxyXG59XHJcblxyXG5cclxuY29uc3QgREVCVUcgOiBib29sZWFuID0gZmFsc2U7XHJcblxyXG5cclxuZXhwb3J0IGNsYXNzIFBpbm5lZENvbHVtbiB7XHJcblxyXG4gIHByaXZhdGUgc3RhdGljIE1JTl9ST1dfSEVJR0hUID0gMjA7XHJcbiAgcHJpdmF0ZSBzdGF0aWMgTUFYX1JPV19IRUlHSFQgPSA1MDA7XHJcbiAgcHJpdmF0ZSBzdGF0aWMgU0VMRUNUSU9OX0NPTE9SID0gQ29sb3JVdGlscy50b1JnYihDb2xvclV0aWxzLmNvbFNlbGVjdGlvbik7IC8vXCJyZ2JhKDIzNywgMjIwLCA4OCwgMC4xNSlcIjtcclxuICBwcml2YXRlIHN0YXRpYyBBQ1RJVkVfQ0VMTF9DT0xPUiA9IENvbG9yVXRpbHMudG9SZ2IoQ29sb3JVdGlscy5jdXJyZW50Um93KTsgLy9cInJnYmEoMTUzLCAyMzcsIDgyLCAwLjI1KVwiO1xyXG4gIHByaXZhdGUgc3RhdGljIFlfUkVTSVpFX1NFTlNJVElWSVRZID0gMjtcclxuXHJcbiAgcHJpdmF0ZSBtX2ZEZXZpY2VQaXhlbFJhdGlvIDogbnVtYmVyO1xyXG4gIHByaXZhdGUgbV9jb2xHcmlkIDogREcuR3JpZENvbHVtbiB8IG51bGw7XHJcbiAgcHJpdmF0ZSBtX3Jvb3QgOiBIVE1MQ2FudmFzRWxlbWVudCB8IG51bGw7XHJcbiAgcHJpdmF0ZSBtX25XaWR0aEJ1ZyA6IG51bWJlcjtcclxuICAvL3ByaXZhdGUgbV9vYnNlcnZlclJlc2l6ZSA6IFJlc2l6ZU9ic2VydmVyIHwgbnVsbDtcclxuICBwcml2YXRlIG1fb2JzZXJ2ZXJSZXNpemVHcmlkIDogUmVzaXplT2JzZXJ2ZXIgfCBudWxsO1xyXG4gIHByaXZhdGUgbV9oYW5kbGVyS2V5RG93biA6IHJ4anMuU3Vic2NyaXB0aW9uIHwgbnVsbDtcclxuICBwcml2YXRlIG1faGFuZGxlckNvbHNSZW1vdmVkIDogcnhqcy5TdWJzY3JpcHRpb24gfCBudWxsO1xyXG4gIHByaXZhdGUgbV9oYW5kbGVyQ29sTmFtZUNoYW5nZWQgOiByeGpzLlN1YnNjcmlwdGlvbiB8IG51bGw7XHJcbiAgcHJpdmF0ZSBtX2hhbmRsZXJWU2Nyb2xsIDogcnhqcy5TdWJzY3JpcHRpb24gfCBudWxsO1xyXG4gIHByaXZhdGUgbV9oYW5kbGVyUm93c0ZpbHRlcmluZyA6IHJ4anMuU3Vic2NyaXB0aW9uIHwgbnVsbDtcclxuICBwcml2YXRlIG1faGFuZGxlckN1cnJSb3cgOiByeGpzLlN1YnNjcmlwdGlvbiB8IG51bGw7XHJcbiAgcHJpdmF0ZSBtX2hhbmRsZXJTZWwgOiByeGpzLlN1YnNjcmlwdGlvbiB8IG51bGw7XHJcbiAgLy9wcml2YXRlIG1faGFuZGxlckZpbHRlciA6IGFueTtcclxuICBwcml2YXRlIG1faGFuZGxlclJvd3NSZXNpemVkIDogcnhqcy5TdWJzY3JpcHRpb24gfCBudWxsO1xyXG4gIHByaXZhdGUgbV9oYW5kbGVyUm93c1NvcnRlZCA6IHJ4anMuU3Vic2NyaXB0aW9uIHwgbnVsbDtcclxuXHJcbiAgcHJpdmF0ZSBtX25IUmVzaXplUm93c0JlZm9yZURyYWcgPSAtMTtcclxuICBwcml2YXRlIG1fblJlc2l6ZVJvd0dyaWREcmFnZ2luZyA9IC0xO1xyXG4gIHByaXZhdGUgbV9uWVJlc2l6ZURyYWdnaW5nQW5jaG9yID0gLTE7XHJcbiAgcHJpdmF0ZSBtX25SZXNpemVSb3dHcmlkTW92aW5nID0gLTE7XHJcblxyXG4gIHByaXZhdGUgbV9uWURyYWdnaW5nQW5jaG9yID0gLTE7XHJcbiAgcHJpdmF0ZSBtX25Sb3dHcmlkRHJhZ2dpbmcgPSAtMTtcclxuXHJcbiAgcHJpdmF0ZSBtX25XaGVlbENvdW50IDogbnVtYmVyID0gMDtcclxuXHJcblxyXG4gIHByaXZhdGUgbV9hclhZTW91c2VPbkNlbGxEb3duID0gWy0yLCAtMl07XHJcbiAgcHJpdmF0ZSBtX2FyWFlNb3VzZU9uQ2VsbFVwID0gWy0xLCAtMV07XHJcbiAgcHJpdmF0ZSBtX2JTb3J0ZWRBc2NlbmRpbmcgOiBib29sZWFuIHwgbnVsbCA9IG51bGw7XHJcblxyXG4gIHByaXZhdGUgbV9jZWxsQ3VycmVudCA6IERHLkdyaWRDZWxsIHwgbnVsbCA9IG51bGw7XHJcblxyXG4gIGNvbnN0cnVjdG9yKGNvbEdyaWQgOiBERy5HcmlkQ29sdW1uKSB7XHJcblxyXG4gICAgTW91c2VEaXNwYXRjaGVyLmNyZWF0ZSgpO1xyXG5cclxuICAgIGNvbnN0IGdyaWQgPSBnZXRHcmlkKGNvbEdyaWQpO1xyXG4gICAgaWYoZ3JpZCA9PT0gbnVsbCkge1xyXG4gICAgICB0aHJvdyBuZXcgRXJyb3IoXCJDb2x1bW4gJ1wiICsgY29sR3JpZC5uYW1lICsgXCInIGlzIG5vdCBhdHRhY2hlZCB0byB0aGUgZ3JpZC5cIik7XHJcbiAgICB9XHJcblxyXG4gICAgaWYoIVBpbm5lZFV0aWxzLmlzUGlubmFibGVDb2x1bW4oY29sR3JpZCkpIHtcclxuICAgICAgdGhyb3cgbmV3IEVycm9yKFwiQ29sdW1uICdcIiArIGNvbEdyaWQubmFtZSArIFwiJyBjYW5ub3QgYmUgcGlubmVkLiBJdCBlaXRoZXIgcGlubmVkIG9yIEhUTUwuXCIpO1xyXG4gICAgfVxyXG5cclxuICAgIC8vbGV0IG5Sb3dNaW4gPSBncmlkLm1pblZpc2libGVSb3c7XHJcbiAgICAvL2xldCBuUm93TWF4ID0gZ3JpZC5tYXhWaXNpYmxlUm93O1xyXG4gICAgLy9sZXQgbkNvbE1pbiA9IGdyaWQubWluVmlzaWJsZUNvbHVtbjtcclxuICAgIC8vbGV0IG5Db2xNYXggPSBncmlkLm1heFZpc2libGVDb2x1bW47XHJcbiAgICAvL2NvbnN0IGl0ID0gZ3JpZC5waW5uZWRSb3dzO1xyXG4gICAgLy9jb25zdCBhciA9IEFycmF5LmZyb20oaXQpO1xyXG5cclxuICAgIHRoaXMubV9mRGV2aWNlUGl4ZWxSYXRpbyA9IHdpbmRvdy5kZXZpY2VQaXhlbFJhdGlvO1xyXG5cclxuICAgIGNvbnN0IGRhcnQgPSBERy50b0RhcnQoZ3JpZCk7XHJcblxyXG4gICAgaWYoZGFydC5tX2FyUGlubmVkQ29scyA9PT0gdW5kZWZpbmVkKVxyXG4gICAgICBkYXJ0Lm1fYXJQaW5uZWRDb2xzID0gW107XHJcblxyXG4gICAgaWYoZGFydC5tX2FyUGlubmVkQ29scy5sZW5ndGggPT09IDAgJiYgIUdyaWRVdGlscy5pc1Jvd0hlYWRlcihjb2xHcmlkKSkge1xyXG4gICAgICBjb25zdCBjb2xHcmlkMCA9IGdyaWQuY29sdW1ucy5ieUluZGV4KDApO1xyXG4gICAgICBpZihjb2xHcmlkMCAhPT0gbnVsbCAmJiBjb2xHcmlkMCAhPT0gdW5kZWZpbmVkKVxyXG4gICAgICBuZXcgUGlubmVkQ29sdW1uKGNvbEdyaWQwKTtcclxuICAgIH1cclxuXHJcbiAgICBjb25zdCBuV1RvdGFsUGlubmVkQ29scyA9IFBpbm5lZFV0aWxzLmdldFRvdGFsUGlubmVkQ29sc1dpZHRoKGdyaWQpO1xyXG4gICAgZGFydC5tX2FyUGlubmVkQ29scy5wdXNoKHRoaXMpO1xyXG5cclxuICAgIGNvbnN0IHZpZXdUYWJsZSA9IGdyaWQudmlldztcclxuICAgIGNvbnN0IGRmcmFtZSA9IGdyaWQuZGF0YUZyYW1lO1xyXG5cclxuICAgIGNvbnN0IG5XID0gY29sR3JpZC53aWR0aDtcclxuICAgIHRoaXMubV9jb2xHcmlkID0gY29sR3JpZDtcclxuICAgIHRoaXMubV9uV2lkdGhCdWcgPSAtMTtcclxuICAgIHRyeSB7XHJcbiAgICAgIGNvbEdyaWQudmlzaWJsZSA9IGZhbHNlO1xyXG4gICAgfVxyXG4gICAgY2F0Y2goZSkge1xyXG4gICAgICAvL0RHIGJ1Z1xyXG4gICAgICBjb25zb2xlLmVycm9yKFwiRVJST1I6IENvdWxkbid0IGhpZGUgY29sdW1uICdcIiArIGNvbEdyaWQubmFtZSArIFwiJyBkdWUgdG8gYSBERyBidWcuIEF0dGVtcHQgdG8gc2V0IHRoZSB3aWR0aCB0byAwXCIpO1xyXG4gICAgICB0cnkge1xyXG4gICAgICAgIHRoaXMubV9uV2lkdGhCdWcgPSBjb2xHcmlkLndpZHRoO1xyXG4gICAgICAgIGNvbEdyaWQud2lkdGggPSAwO1xyXG4gICAgICB9IGNhdGNoIChlKSB7XHJcbiAgICAgICAgLy9ERyBidWdcclxuICAgICAgICBjb25zb2xlLmVycm9yKFwiRVJST1I6IENvdWxkbid0IHNldCB0aGUgd2lkdGggdG8gMCBmb3IgY29sdW1uICdcIiArIGNvbEdyaWQubmFtZSArIFwiJyBkdWUgdG8gYSBERyBidWcuIFRoaXMgY291bGQgYmUgaWdub3JlZCBpZiB0aGUgY29sdW1uIHZpc3VhbGx5IGxvb2tzIG9rLlwiKTtcclxuICAgICAgfVxyXG4gICAgfVxyXG5cclxuICAgIGlmKCFHcmlkVXRpbHMuaXNSb3dIZWFkZXIoY29sR3JpZCkpIHtcclxuICAgICAgaWYgKGNvbEdyaWQuc2V0dGluZ3MgPT09IG51bGwgfHwgY29sR3JpZC5zZXR0aW5ncyA9PT0gdW5kZWZpbmVkKVxyXG4gICAgICAgIGNvbEdyaWQuc2V0dGluZ3MgPSB7fTtcclxuXHJcbiAgICAgIGNvbEdyaWQuc2V0dGluZ3MuaXNQaW5uZWQgPSB0cnVlOyAvL3RoaXMgd2lsbCBiZSBzYXZlZCB3aXRoIHRoZSBsYXlvdXRcclxuICAgICAgY29sR3JpZC5zZXR0aW5ncy5pZHhQaW5uZWQgPSBkYXJ0Lm1fYXJQaW5uZWRDb2xzLmxlbmd0aCAtIDE7XHJcbiAgICB9XHJcblxyXG4gICAgZ3JpZC5jYW52YXMuc3R5bGUubGVmdCA9IChncmlkLmNhbnZhcy5vZmZzZXRMZWZ0ICsgblcpLnRvU3RyaW5nKCkgKyBcInB4XCI7XHJcbiAgICBncmlkLm92ZXJsYXkuc3R5bGUubGVmdD0gKGdyaWQub3ZlcmxheS5vZmZzZXRMZWZ0ICsgblcpLnRvU3RyaW5nKCkgKyBcInB4XCI7XHJcblxyXG4gICAgZ3JpZC5jYW52YXMuc3R5bGUud2lkdGggPSAoZ3JpZC5jYW52YXMub2Zmc2V0V2lkdGggLSBuVykudG9TdHJpbmcoKSArIFwicHhcIjtcclxuICAgIGdyaWQub3ZlcmxheS5zdHlsZS53aWR0aD0gKGdyaWQub3ZlcmxheS5vZmZzZXRXaWR0aCAtIG5XKS50b1N0cmluZygpICsgXCJweFwiO1xyXG5cclxuICAgIGNvbnN0IG5IZWlnaHQgPSBncmlkLmNhbnZhcy5oZWlnaHQ7Ly9jYW52YXMgcGl4ZWwgaGVpZ2h0XHJcbiAgICBjb25zdCBlQ2FudmFzVGhpcyA9IHVpLmNhbnZhcyhuVyp3aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbywgbkhlaWdodCk7XHJcbiAgICBjb25zdCB0YWJJbmRleCA9ICBncmlkLmNhbnZhcy5nZXRBdHRyaWJ1dGUoXCJ0YWJJbmRleFwiKTtcclxuICAgIGlmKHRhYkluZGV4ICE9PSBudWxsKVxyXG4gICAgIGVDYW52YXNUaGlzLnNldEF0dHJpYnV0ZShcInRhYkluZGV4XCIsIHRhYkluZGV4KTtcclxuXHJcbiAgICBlQ2FudmFzVGhpcy5zdHlsZS5wb3NpdGlvbiA9IFwiYWJzb2x1dGVcIjtcclxuICAgIGVDYW52YXNUaGlzLnN0eWxlLmxlZnQgPSBuV1RvdGFsUGlubmVkQ29scyArIFwicHhcIjtcclxuICAgIGVDYW52YXNUaGlzLnN0eWxlLnRvcCA9IGdyaWQuY2FudmFzLm9mZnNldFRvcCArIFwicHhcIjtcclxuICAgIGVDYW52YXNUaGlzLnN0eWxlLndpZHRoID0gblcgKyBcInB4XCI7XHJcbiAgICBlQ2FudmFzVGhpcy5zdHlsZS5oZWlnaHQgPSBNYXRoLnJvdW5kKG5IZWlnaHQvd2luZG93LmRldmljZVBpeGVsUmF0aW8pICsgXCJweFwiO1xyXG5cclxuICAgIC8vY29uc29sZS5sb2coXCJoIFwiICsgZ3JpZC5jYW52YXMuaGVpZ2h0ICsgXCIgb2Zmc2V0IFwiICsgZ3JpZC5jYW52YXMub2Zmc2V0SGVpZ2h0KTtcclxuXHJcbiAgICBpZihncmlkLmNhbnZhcy5wYXJlbnROb2RlID09PSBudWxsKVxyXG4gICAgICB0aHJvdyBuZXcgRXJyb3IoXCJQYXJlbnQgbm9kZSBmb3IgY2FudmFzIGNhbm5vdCBiZSBudWxsLlwiKTtcclxuXHJcbiAgICBncmlkLmNhbnZhcy5wYXJlbnROb2RlLmluc2VydEJlZm9yZShlQ2FudmFzVGhpcywgZ3JpZC5jYW52YXMpO1xyXG4gICAgdGhpcy5tX3Jvb3QgPSBlQ2FudmFzVGhpcztcclxuXHJcblxyXG4gICAgY29uc3QgY29sR3JpZDAgPSBncmlkLmNvbHVtbnMuYnlJbmRleCgwKTtcclxuICAgIGlmKGNvbEdyaWQwICE9PSBudWxsICYmIGNvbEdyaWQwICE9PSB1bmRlZmluZWQpIHsvL0RHIEJ1ZyBmcm9tIHJlYWRpbmcgbGF5b3V0XHJcbiAgICB0cnl7XHJcbiAgICAgICAgY29sR3JpZDAudmlzaWJsZSA9IGZhbHNlO1xyXG4gICAgICB9XHJcbiAgICAgIGNhdGNoKGUpIHtcclxuICAgICAgICBjb25zb2xlLmVycm9yKFwiRVJST1I6IENvdWxkbid0IGhpZGUgcm93IGhlYWRlci5cIik7XHJcbiAgICAgIH1cclxuICAgIH1cclxuXHJcblxyXG4gICAgLy9PblJlc2l6ZSBSb3cgaGVhZGVyXHJcbiAgICBjb25zdCBoZWFkZXJUaGlzID0gdGhpczsvKlxyXG4gICAgdGhpcy5tX29ic2VydmVyUmVzaXplID0gbmV3IFJlc2l6ZU9ic2VydmVyKGVudHJpZXMgPT4ge1xyXG4gICAgICBjb25zdCBnID0gaGVhZGVyVGhpcy5tX3Jvb3QuZ2V0Q29udGV4dCgnMmQnKTtcclxuICAgICAgZm9yIChsZXQgZW50cnkgb2YgZW50cmllcykge1xyXG4gICAgICAgIGhlYWRlclRoaXMucGFpbnQoZywgZ3JpZCk7XHJcbiAgICAgIH1cclxuICAgIH0pO1xyXG4gICAgdGhpcy5tX29ic2VydmVyUmVzaXplLm9ic2VydmUoaGVhZGVyVGhpcy5tX3Jvb3QpOyovXHJcblxyXG5cclxuXHJcbiAgICAvL09uUmVzaXplIEdyaWRcclxuICAgIHRoaXMubV9vYnNlcnZlclJlc2l6ZUdyaWQgPSBuZXcgUmVzaXplT2JzZXJ2ZXIoZnVuY3Rpb24gKGVudHJpZXMgOiBhbnkpIHtcclxuXHJcbiAgICAgIGNvbnN0IGJDdXJyZW50ID0gIERHLnRvRGFydChncm9rLnNoZWxsLnYpID09PSBERy50b0RhcnQodmlld1RhYmxlKTtcclxuICAgICAgaWYoIWJDdXJyZW50KVxyXG4gICAgICAgIHJldHVybjtcclxuXHJcbiAgICAgIGlmKGhlYWRlclRoaXMubV9mRGV2aWNlUGl4ZWxSYXRpbyAhPT0gd2luZG93LmRldmljZVBpeGVsUmF0aW8gfHwgZ3JpZC5jYW52YXMuaGVpZ2h0ICE9PSBlQ2FudmFzVGhpcy5oZWlnaHQpIHtcclxuICAgICAgICBlQ2FudmFzVGhpcy53aWR0aCA9IG5XKndpbmRvdy5kZXZpY2VQaXhlbFJhdGlvO1xyXG4gICAgICAgIGVDYW52YXNUaGlzLmhlaWdodCA9IGdyaWQuY2FudmFzLmhlaWdodDtcclxuICAgICAgICBlQ2FudmFzVGhpcy5zdHlsZS50b3AgPSBncmlkLmNhbnZhcy5vZmZzZXRUb3AgKyBcInB4XCI7XHJcbiAgICAgICAgZUNhbnZhc1RoaXMuc3R5bGUud2lkdGggPSBuVyArIFwicHhcIjtcclxuICAgICAgICBlQ2FudmFzVGhpcy5zdHlsZS5oZWlnaHQgPSBNYXRoLnJvdW5kKGdyaWQuY2FudmFzLmhlaWdodC93aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbykgKyBcInB4XCI7XHJcblxyXG4gICAgICAgIGhlYWRlclRoaXMubV9mRGV2aWNlUGl4ZWxSYXRpbyA9IHdpbmRvdy5kZXZpY2VQaXhlbFJhdGlvO1xyXG4gICAgICB9XHJcblxyXG4gICAgICAvL2NvbnNvbGUubG9nKFwiR3JpZCBSZXNpemU6IFwiICsgZ3JpZC5jYW52YXMuaGVpZ2h0ICsgXCIgXCIgKyB3aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbyk7XHJcbiAgICAgIC8vZUNhbnZhc1RoaXMuc3R5bGUuaGVpZ2h0ID0gZ3JpZC5yb290LnN0eWxlLmhlaWdodDtcclxuLypcclxuICAgICAgY29uc3QgZUNhbnZhc05ldyA9IHVpLmNhbnZhcyhuVywgZ3JpZC5yb290Lm9mZnNldEhlaWdodCk7XHJcbiAgICAgIGlmKGhlYWRlclRoaXMubV9yb290LnBhcmVudE5vZGUgIT09IG51bGwpIHtcclxuICAgICAgICBoZWFkZXJUaGlzLm1fcm9vdC5wYXJlbnROb2RlLnJlcGxhY2VDaGlsZChlQ2FudmFzTmV3LCBoZWFkZXJUaGlzLm1fcm9vdCk7XHJcbiAgICAgICAgaGVhZGVyVGhpcy5tX3Jvb3QgPSBlQ2FudmFzTmV3O1xyXG4gICAgICB9Ki9cclxuICAgICAgLy9oZWFkZXJUaGlzLm1fcm9vdC5oZWlnaHQgPSBncmlkLnJvb3Qub2Zmc2V0SGVpZ2h0O1xyXG4gICAgICBjb25zdCBnID0gZUNhbnZhc1RoaXMuZ2V0Q29udGV4dCgnMmQnKTtcclxuICAgICAgZm9yIChsZXQgZW50cnkgb2YgZW50cmllcykge1xyXG4gICAgICAgIHNldFRpbWVvdXQoKCk9PiB7aGVhZGVyVGhpcy5wYWludChnLCBncmlkKTt9LCAxMDApO1xyXG4gICAgICB9XHJcbiAgICB9KTtcclxuXHJcbiAgICB0aGlzLm1fb2JzZXJ2ZXJSZXNpemVHcmlkPy5vYnNlcnZlKGdyaWQuY2FudmFzKTtcclxuXHJcblxyXG4gICAgdGhpcy5tX2hhbmRsZXJLZXlEb3duID0gcnhqcy5mcm9tRXZlbnQ8S2V5Ym9hcmRFdmVudD4oZUNhbnZhc1RoaXMsICdrZXlkb3duJykuc3Vic2NyaWJlKChlIDogS2V5Ym9hcmRFdmVudCkgPT4ge1xyXG5cclxuICAgICAgLy9hbGVydCgndXAnKTtcclxuICAgICAgc2V0VGltZW91dCgoKSA9PntcclxuICAgICAgICBjb25zdCBlZSA9IG5ldyBLZXlib2FyZEV2ZW50KGUudHlwZSwgZSk7XHJcbiAgICAgICAgdHJ5e2dyaWQub3ZlcmxheS5kaXNwYXRjaEV2ZW50KGVlKTt9XHJcbiAgICAgICAgY2F0Y2goZXgpIHtcclxuICAgICAgICAgIC8vY29uc29sZS5lcnJvcihleC5tZXNzYWdlKTtcclxuICAgICAgICB9XHJcbiAgICAgIH0sIDEpO1xyXG5cclxuICAgIH0pO1xyXG5cclxuXHJcbiAgICBjb25zdCBzY3JvbGxWZXJ0ID0gZ3JpZC52ZXJ0U2Nyb2xsO1xyXG4gICAgdGhpcy5tX2hhbmRsZXJWU2Nyb2xsID0gc2Nyb2xsVmVydC5vblZhbHVlc0NoYW5nZWQuc3Vic2NyaWJlKCgpID0+IHtcclxuICAgICAgY29uc3QgZyA9IGVDYW52YXNUaGlzLmdldENvbnRleHQoJzJkJyk7XHJcbiAgICAgIGhlYWRlclRoaXMucGFpbnQoZywgZ3JpZCk7XHJcbiAgICB9KTtcclxuXHJcbiAgICB0aGlzLm1faGFuZGxlclJvd3NGaWx0ZXJpbmcgPSBkZnJhbWUub25Sb3dzRmlsdGVyaW5nLnN1YnNjcmliZSgoKSA9PiB7XHJcbiAgICAgIHNldFRpbWVvdXQoKCkgPT4ge1xyXG4gICAgICAgIGNvbnN0IGcgPSBlQ2FudmFzVGhpcy5nZXRDb250ZXh0KCcyZCcpO1xyXG4gICAgICAgIGhlYWRlclRoaXMucGFpbnQoZywgZ3JpZCk7XHJcbiAgICAgIH0sIDEwMCk7XHJcblxyXG4gICAgfSk7XHJcblxyXG4gICAgdGhpcy5tX2hhbmRsZXJDdXJyUm93ID0gZGZyYW1lLm9uQ3VycmVudFJvd0NoYW5nZWQuc3Vic2NyaWJlKCgpID0+IHtcclxuICAgICAgICBjb25zdCBnID0gZUNhbnZhc1RoaXMuZ2V0Q29udGV4dCgnMmQnKTtcclxuICAgICAgICBoZWFkZXJUaGlzLnBhaW50KGcsIGdyaWQpO1xyXG4gICAgICB9XHJcbiAgICApO1xyXG5cclxuICAgIHRoaXMubV9oYW5kbGVyU2VsID0gZGZyYW1lLm9uU2VsZWN0aW9uQ2hhbmdlZC5zdWJzY3JpYmUoKGUgOiBhbnkpID0+IHtcclxuICAgICAgICBjb25zdCBnID0gZUNhbnZhc1RoaXMuZ2V0Q29udGV4dCgnMmQnKTtcclxuICAgICAgICBoZWFkZXJUaGlzLnBhaW50KGcsIGdyaWQpO1xyXG4gICAgICB9XHJcbiAgICApO1xyXG5cclxuICAgIHRoaXMubV9oYW5kbGVyQ29sc1JlbW92ZWQgPSBkZnJhbWUub25Db2x1bW5zUmVtb3ZlZC5zdWJzY3JpYmUoKGUgOiBDb2x1bW5zQXJncykgPT4ge1xyXG5cclxuICAgICAgICAgIGlmKGhlYWRlclRoaXMubV9jb2xHcmlkID09PSBudWxsKVxyXG4gICAgICAgICAgICByZXR1cm47XHJcbiAgICAgICAgICBmb3IobGV0IG5DPTA7IG5DPGUuY29sdW1ucy5sZW5ndGg7ICsrbkMpIHtcclxuICAgICAgICAgICAgaWYoZS5jb2x1bW5zW25DXS5uYW1lID09PSBoZWFkZXJUaGlzLm1fY29sR3JpZC5uYW1lKVxyXG4gICAgICAgICAgICAgIGhlYWRlclRoaXMuY2xvc2UoKTtcclxuICAgICAgICAgIH1cclxuICAgICAgICB9XHJcbiAgICApO1xyXG5cclxuICAgIHRoaXMubV9oYW5kbGVyQ29sTmFtZUNoYW5nZWQgPSBkZnJhbWUub25Db2x1bW5OYW1lQ2hhbmdlZC5zdWJzY3JpYmUoKGUgOiBhbnkpID0+IHtcclxuXHJcbiAgICAgICAgICBjb25zdCBkYXJ0ID0gdG9EYXJ0KGUpO1xyXG4gICAgICAgICAgY29uc3Qgc3RyQ29sTmFtZU9sZCA9IGRhcnQubmV3TmFtZTtcclxuICAgICAgICAgIGlmKHN0ckNvbE5hbWVPbGQgPT09IGhlYWRlclRoaXMubV9jb2xHcmlkPy5uYW1lKSB7XHJcbiAgICAgICAgICAgIGNvbnN0IGcgPSBlQ2FudmFzVGhpcy5nZXRDb250ZXh0KCcyZCcpO1xyXG4gICAgICAgICAgICBoZWFkZXJUaGlzLnBhaW50KGcsIGdyaWQpO1xyXG4gICAgICAgICAgfVxyXG4gICAgICAgIH1cclxuICAgICk7XHJcblxyXG5cclxuLypcclxuICAgIHRoaXMubV9oYW5kbGVyRmlsdGVyID0gZGZyYW1lLm9uUm93c0ZpbHRlcmVkLnN1YnNjcmliZSgoZSA6IGFueSkgPT4ge1xyXG4gICAgICAgIGNvbnN0IGcgPSBlQ2FudmFzVGhpcy5nZXRDb250ZXh0KCcyZCcpO1xyXG4gICAgICAgIGhlYWRlclRoaXMucGFpbnQoZywgZ3JpZCk7XHJcbiAgICAgIH1cclxuICAgICk7XHJcbiovXHJcblxyXG4gICAgdGhpcy5tX2hhbmRsZXJSb3dzUmVzaXplZCA9IGdyaWQub25Sb3dzUmVzaXplZC5zdWJzY3JpYmUoKGUgOiBhbnkpID0+IHtcclxuICAgICAgICBjb25zdCBnID0gZUNhbnZhc1RoaXMuZ2V0Q29udGV4dCgnMmQnKTtcclxuICAgICAgICBoZWFkZXJUaGlzLnBhaW50KGcsIGdyaWQpO1xyXG4gICAgICB9XHJcbiAgICApO1xyXG5cclxuICAgIHRoaXMubV9oYW5kbGVyUm93c1NvcnRlZCA9IGdyaWQub25Sb3dzU29ydGVkLnN1YnNjcmliZSgoZSA6IGFueSkgPT4ge1xyXG4gICAgICAgIGNvbnN0IGcgPSBlQ2FudmFzVGhpcy5nZXRDb250ZXh0KCcyZCcpO1xyXG4gICAgICAgIGhlYWRlclRoaXMucGFpbnQoZywgZ3JpZCk7XHJcbiAgICAgIH1cclxuICAgICk7XHJcbiAgfVxyXG5cclxuICBpc1Bpbm5lZCgpIDogYm9vbGVhbiB7XHJcbiAgICByZXR1cm4gdGhpcy5tX2NvbEdyaWQgIT09IG51bGw7XHJcbiAgfVxyXG5cclxuICBnZXRHcmlkQ29sdW1uKCkgOiBERy5HcmlkQ29sdW1uIHwgbnVsbHtcclxuICAgIHJldHVybiB0aGlzLm1fY29sR3JpZDtcclxuICB9XHJcblxyXG4gIGdldFdpZHRoKCkgOiBudW1iZXIge1xyXG4gICAgcmV0dXJuIHRoaXMubV9yb290ID09PSBudWxsID8gLTEgOiB0aGlzLm1fcm9vdC5vZmZzZXRXaWR0aDtcclxuICB9XHJcblxyXG4gIGdldFJvb3QoKSA6IEhUTUxDYW52YXNFbGVtZW50IHwgbnVsbCB7XHJcbiAgICByZXR1cm4gdGhpcy5tX3Jvb3Q7XHJcbiAgfVxyXG5cclxuICBwdWJsaWMgY2xvc2UoKSA6IHZvaWQge1xyXG5cclxuICAgIGlmKHRoaXMubV9jb2xHcmlkID09PSBudWxsKSB7XHJcbiAgICAgIHRocm93IG5ldyBFcnJvcihcIkNvbHVtbiBoYXMgYWxyZWFkeSBiZWVuIHVucGlubmVkXCIpO1xyXG4gICAgfVxyXG5cclxuICAgIGlmKHRoaXMubV9vYnNlcnZlclJlc2l6ZUdyaWQgIT09IG51bGwpIHtcclxuICAgICAgdGhpcy5tX29ic2VydmVyUmVzaXplR3JpZC5kaXNjb25uZWN0KCk7XHJcbiAgICAgIHRoaXMubV9vYnNlcnZlclJlc2l6ZUdyaWQgPSBudWxsO1xyXG4gICAgfVxyXG4vKm15IGNoYW5nZXNcclxuICAgIGlmKHRoaXMubV9vYnNlcnZlclJlc2l6ZSAhPT0gbnVsbCkge1xyXG4gICAgICB0aGlzLm1fb2JzZXJ2ZXJSZXNpemUuZGlzY29ubmVjdCgpO1xyXG4gICAgICB0aGlzLm1fb2JzZXJ2ZXJSZXNpemUgPSBudWxsO1xyXG4gICAgfVxyXG4gICAgKi9cclxuXHJcbiAgICB0aGlzLm1faGFuZGxlcktleURvd24/LnVuc3Vic2NyaWJlKCk7XHJcbiAgICB0aGlzLm1faGFuZGxlcktleURvd24gPSBudWxsO1xyXG5cclxuICAgIHRoaXMubV9oYW5kbGVyQ29sc1JlbW92ZWQ/LnVuc3Vic2NyaWJlKCk7XHJcbiAgICB0aGlzLm1faGFuZGxlckNvbHNSZW1vdmVkID0gbnVsbDtcclxuXHJcbiAgICB0aGlzLm1faGFuZGxlckNvbE5hbWVDaGFuZ2VkPy51bnN1YnNjcmliZSgpO1xyXG4gICAgdGhpcy5tX2hhbmRsZXJDb2xOYW1lQ2hhbmdlZCA9IG51bGw7XHJcblxyXG4gICAgdGhpcy5tX2hhbmRsZXJWU2Nyb2xsPy51bnN1YnNjcmliZSgpO1xyXG4gICAgdGhpcy5tX2hhbmRsZXJWU2Nyb2xsID0gbnVsbDtcclxuXHJcbiAgICB0aGlzLm1faGFuZGxlclJvd3NSZXNpemVkPy51bnN1YnNjcmliZSgpO1xyXG4gICAgdGhpcy5tX2hhbmRsZXJSb3dzUmVzaXplZCA9IG51bGw7XHJcblxyXG4gICAgdGhpcy5tX2hhbmRsZXJSb3dzU29ydGVkPy51bnN1YnNjcmliZSgpO1xyXG4gICAgdGhpcy5tX2hhbmRsZXJSb3dzU29ydGVkID0gbnVsbDtcclxuXHJcbiAgICB0aGlzLm1faGFuZGxlclJvd3NGaWx0ZXJpbmc/LnVuc3Vic2NyaWJlKCk7XHJcbiAgICB0aGlzLm1faGFuZGxlclJvd3NGaWx0ZXJpbmcgPSBudWxsO1xyXG5cclxuICAgIHRoaXMubV9oYW5kbGVyQ3VyclJvdz8udW5zdWJzY3JpYmUoKTtcclxuICAgIHRoaXMubV9oYW5kbGVyQ3VyclJvdyA9IG51bGw7XHJcblxyXG4gICAgdGhpcy5tX2hhbmRsZXJTZWw/LnVuc3Vic2NyaWJlKCk7XHJcbiAgICB0aGlzLm1faGFuZGxlclNlbCA9IG51bGw7XHJcblxyXG4gICAgY29uc3QgZ3JpZCA9IGdldEdyaWQodGhpcy5tX2NvbEdyaWQpO1xyXG4gICAgaWYoZ3JpZCA9PT0gbnVsbCl7XHJcbiAgICAgIHRocm93IG5ldyBFcnJvcihcIkNvbHVtbiAnXCIgKyB0aGlzLm1fY29sR3JpZC5uYW1lICsgXCInIGlzIGRpc2Nvbm5lY3RlZCBmcm9tIGdyaWQuXCIpO1xyXG4gICAgfVxyXG5cclxuICAgIGNvbnN0IGRhcnQgPSBERy50b0RhcnQoZ3JpZCk7XHJcbiAgICBjb25zdCBhciA9IGRhcnQubV9hclBpbm5lZENvbHM7XHJcbiAgICBjb25zdCBuSWR4ID0gYXIuaW5kZXhPZih0aGlzKTtcclxuICAgIGFyLnNwbGljZShuSWR4LCAxKTtcclxuXHJcbiAgICBpZih0aGlzLm1fcm9vdCA9PT0gbnVsbClcclxuICAgICAgdGhyb3cgbmV3IEVycm9yKCdSb290IGNhbm5vdCBiZSBudWxsJyk7XHJcblxyXG4gICAgbGV0IG5JZHhQaW5uZWQgPSAtMTtcclxuICAgIGxldCBjb2xHcmlkVG1wPSBudWxsO1xyXG4gICAgZm9yKGxldCBuPW5JZHg7IG48YXIubGVuZ3RoOyArK24pIHtcclxuICAgICAgY29sR3JpZFRtcCA9IGFyW25dO1xyXG4gICAgICBjb2xHcmlkVG1wLm1fcm9vdC5zdHlsZS5sZWZ0ID0gKGNvbEdyaWRUbXAubV9yb290Lm9mZnNldExlZnQgLSB0aGlzLm1fcm9vdC5vZmZzZXRXaWR0aCkudG9TdHJpbmcoKSArIFwicHhcIjtcclxuXHJcbiAgICAgIG5JZHhQaW5uZWQgPSAgY29sR3JpZFRtcC5tX2NvbEdyaWQuc2V0dGluZ3MuaWR4UGlubmVkO1xyXG4gICAgICBjb2xHcmlkVG1wLm1fY29sR3JpZC5zZXR0aW5ncy5pZHhQaW5uZWQgPSBuO1xyXG4gICAgfVxyXG5cclxuICAgIGlmKCFHcmlkVXRpbHMuaXNSb3dIZWFkZXIodGhpcy5tX2NvbEdyaWQpKSB7XHJcbiAgICAgIHRoaXMubV9jb2xHcmlkLnNldHRpbmdzLmlkeFBpbm5lZCA9IC0xO1xyXG4gICAgICB0aGlzLm1fY29sR3JpZC5zZXR0aW5ncy5pc1Bpbm5lZCA9IGZhbHNlO1xyXG4gICAgfVxyXG5cclxuXHJcbiAgICBpZih0aGlzLm1fbldpZHRoQnVnID49IDApIHtcclxuICAgICAgdHJ5IHtcclxuICAgICAgICB0aGlzLm1fY29sR3JpZC53aWR0aCA9IHRoaXMubV9uV2lkdGhCdWc7XHJcbiAgICAgIH1cclxuICAgICAgY2F0Y2goZSkge1xyXG4gICAgICAgIC8vREcgYnVnXHJcbiAgICAgICAgY29uc29sZS5lcnJvcihcIkVSUk9SOiBDb3VsZG4ndCBzZXQgdGhlIHdpZHRoIHRvIFwiICsgdGhpcy5tX25XaWR0aEJ1ZyArIFwiIGZvciBjb2x1bW4gJ1wiICsgdGhpcy5tX2NvbEdyaWQubmFtZSArIFwiJyBkdWUgdG8gYSBERyBidWcuIFRoaXMgY291bGQgYmUgaWdub3JlZCBpZiB0aGUgY29sdW1uIHZpc3VhbGx5IGxvb2tzIG9rLlwiKTtcclxuICAgICAgfVxyXG4gICAgfVxyXG5cclxuICAgIHRyeSB7XHJcbiAgICAgIHRoaXMubV9jb2xHcmlkLndpZHRoID0gdGhpcy5tX3Jvb3Qub2Zmc2V0V2lkdGg7XHJcbiAgICAgIHRoaXMubV9jb2xHcmlkLnZpc2libGUgPSB0cnVlO1xyXG4gICAgfVxyXG4gICAgY2F0Y2goZSkge1xyXG4gICAgICAvL0RHIGJ1Z1xyXG4gICAgICBjb25zb2xlLmVycm9yKFwiRVJST1I6IENvdWxkbid0IHNob3cgY29sdW1uICdcIiArIHRoaXMubV9jb2xHcmlkLm5hbWUgKyBcIicgZHVlIHRvIGEgREcgYnVnLiBUaGlzIGNvdWxkIGJlIGlnbm9yZWQgaWYgdGhlIGNvbHVtbiB2aXN1YWxseSBsb29rcyBvay5cIik7XHJcbiAgICB9XHJcblxyXG4gICAgZ3JpZC5jYW52YXMuc3R5bGUubGVmdCA9IChncmlkLmNhbnZhcy5vZmZzZXRMZWZ0IC0gdGhpcy5tX3Jvb3Qub2Zmc2V0V2lkdGgpLnRvU3RyaW5nKCkgKyBcInB4XCI7XHJcbiAgICBncmlkLm92ZXJsYXkuc3R5bGUubGVmdD0gKGdyaWQub3ZlcmxheS5vZmZzZXRMZWZ0IC0gdGhpcy5tX3Jvb3Qub2Zmc2V0V2lkdGgpLnRvU3RyaW5nKCkgKyBcInB4XCI7XHJcbiAgICBncmlkLmNhbnZhcy5zdHlsZS53aWR0aCA9IChncmlkLmNhbnZhcy5vZmZzZXRXaWR0aCArIHRoaXMubV9yb290Lm9mZnNldFdpZHRoKS50b1N0cmluZygpICsgXCJweFwiO1xyXG4gICAgZ3JpZC5vdmVybGF5LnN0eWxlLndpZHRoPSAoZ3JpZC5vdmVybGF5Lm9mZnNldFdpZHRoICsgdGhpcy5tX3Jvb3Qub2Zmc2V0V2lkdGgpLnRvU3RyaW5nKCkgKyBcInB4XCI7XHJcblxyXG4gICAgaWYodGhpcy5tX3Jvb3QucGFyZW50Tm9kZSAhPT0gbnVsbClcclxuICAgICB0aGlzLm1fcm9vdC5wYXJlbnROb2RlLnJlbW92ZUNoaWxkKHRoaXMubV9yb290KTtcclxuXHJcbiAgICB0aGlzLm1fcm9vdCA9IG51bGw7XHJcblxyXG4gICAgaWYgKGRhcnQubV9hclBpbm5lZENvbHMubGVuZ3RoID09PSAxICYmIGRhcnQubV9hclBpbm5lZENvbHNbMF0ubV9jb2xHcmlkLmlkeCA9PT0gMCAmJiB0aGlzLm1fY29sR3JpZC5pZHggIT09IDApIHtcclxuXHJcbiAgICAgICAgLy8gdHJ5e2NvbEdyaWQwLnZpc2libGUgPSB0cnVlO31cclxuICAgICAgICB0cnkge1xyXG4gICAgICAgICAgZGFydC5tX2FyUGlubmVkQ29sc1swXS5jbG9zZSgpO1xyXG4gICAgICAgIH0gY2F0Y2ggKGUpIHtcclxuICAgICAgICAgIGNvbnNvbGUuZXJyb3IoXCJFUlJPUjogQ291bGRuJ3QgY2xvc2UgcGlubmVkIGNvbHVtbiAnXCIgKyBkYXJ0Lm1fYXJQaW5uZWRDb2xzWzBdLm1fY29sR3JpZC5uYW1lICsgXCInIFwiKTtcclxuICAgICAgICB9XHJcbiAgICB9XHJcbiAgICB0aGlzLm1fY29sR3JpZCA9IG51bGw7XHJcbiAgfVxyXG5cclxuXHJcbiAgcHVibGljIG9uTW91c2VFbnRlcihlIDogTW91c2VFdmVudCkgOiB2b2lkIHtcclxuICAgIGlmKERFQlVHKVxyXG4gICAgICBjb25zb2xlLmxvZygnTW91c2UgRW50ZXIgUGlubmVkIENvbHVtbjogJyArIHRoaXMuZ2V0R3JpZENvbHVtbigpPy5uYW1lKTtcclxuICB9XHJcblxyXG4gIHB1YmxpYyBvbk1vdXNlTW92ZShlIDogTW91c2VFdmVudCkgOiB2b2lkIHtcclxuICAgIGlmKERFQlVHKVxyXG4gICAgICBjb25zb2xlLmxvZygnTW91c2UgTW92ZSBQaW5uZWQgQ29sdW1uOiAnICsgdGhpcy5nZXRHcmlkQ29sdW1uKCk/Lm5hbWUpO1xyXG5cclxuICAgIGlmKHRoaXMubV9jb2xHcmlkID09PSBudWxsIHx8IHRoaXMubV9yb290ID09PSBudWxsKVxyXG4gICAgICByZXR1cm47XHJcblxyXG4gICAgY29uc3QgZ3JpZCA9IHRoaXMubV9jb2xHcmlkLmdyaWQ7XHJcbiAgICBjb25zdCB2aWV3VGFibGUgPSBncmlkLnZpZXc7XHJcblxyXG4gICAgaWYoREcudG9EYXJ0KGdyb2suc2hlbGwudikgIT09IERHLnRvRGFydCh2aWV3VGFibGUpKSB7XHJcbiAgICAgIHJldHVybjtcclxuICAgIH1cclxuXHJcblxyXG4gICAgY29uc3QgYXJYWU9uQ2VsbCA9IFstMSwtMV07XHJcblxyXG4gICAgbGV0IG5Sb3dHcmlkID0gUGlubmVkQ29sdW1uLmhpdFRlc3RSb3dzKHRoaXMubV9yb290LCBncmlkLCBlLCBmYWxzZSwgYXJYWU9uQ2VsbCk7XHJcbiAgICBpZihuUm93R3JpZCA+PSAwKSB7XHJcbiAgICAgIGNvbnN0IGNlbGwgPSBncmlkLmNlbGwodGhpcy5tX2NvbEdyaWQubmFtZSwgblJvd0dyaWQpO1xyXG4gICAgICBjb25zdCByZW5kZXJlciA9IGdldFJlbmRlcmVyKGNlbGwpO1xyXG5cclxuICAgICAgaWYgKHJlbmRlcmVyIGluc3RhbmNlb2YgR3JpZENlbGxSZW5kZXJlckV4KSB7XHJcblxyXG4gICAgICAgIGlmICh0aGlzLm1fY2VsbEN1cnJlbnQgPT09IG51bGwpIHtcclxuICAgICAgICAgIHJlbmRlcmVyLm9uTW91c2VFbnRlckV4KGNlbGwsIGUsIGFyWFlPbkNlbGxbMF0sIGFyWFlPbkNlbGxbMV0pO1xyXG4gICAgICAgIH1cclxuXHJcbiAgICAgICAgaWYgKHRoaXMubV9jZWxsQ3VycmVudCAhPT0gbnVsbCAmJiBuUm93R3JpZCAhPT0gdGhpcy5tX2NlbGxDdXJyZW50LmdyaWRSb3cpIHtcclxuICAgICAgICAgIHJlbmRlcmVyLm9uTW91c2VMZWF2ZUV4KHRoaXMubV9jZWxsQ3VycmVudCwgZSwgLTEsIC0xKTtcclxuXHJcbiAgICAgICAgICByZW5kZXJlci5vbk1vdXNlRW50ZXJFeChjZWxsLCBlLCBhclhZT25DZWxsWzBdLCBhclhZT25DZWxsWzFdKTtcclxuICAgICAgICB9XHJcblxyXG4gICAgICAgIHJlbmRlcmVyLm9uTW91c2VNb3ZlRXgoY2VsbCwgZSwgYXJYWU9uQ2VsbFswXSwgYXJYWU9uQ2VsbFsxXSk7XHJcbiAgICAgIH1cclxuXHJcbiAgICAgIHRoaXMubV9jZWxsQ3VycmVudCA9IGNlbGw7XHJcbiAgICB9XHJcbiAgICBlbHNlIGlmICh0aGlzLm1fY2VsbEN1cnJlbnQgIT09IG51bGwpIHtcclxuICAgICAgY29uc3QgcmVuZGVyZXIgPSBnZXRSZW5kZXJlcih0aGlzLm1fY2VsbEN1cnJlbnQpO1xyXG4gICAgICBpZiAocmVuZGVyZXIgaW5zdGFuY2VvZiBHcmlkQ2VsbFJlbmRlcmVyRXgpIHtcclxuICAgICAgICByZW5kZXJlci5vbk1vdXNlTGVhdmVFeCh0aGlzLm1fY2VsbEN1cnJlbnQsIGUsIC0xLCAtMSk7XHJcbiAgICAgIH1cclxuXHJcbiAgICAgIHRoaXMubV9jZWxsQ3VycmVudCA9IG51bGw7XHJcbiAgICB9XHJcblxyXG4gICAgblJvd0dyaWQgPSBQaW5uZWRDb2x1bW4uaGl0VGVzdFJvd3ModGhpcy5tX3Jvb3QsIGdyaWQsIGUsIHRydWUsIHVuZGVmaW5lZCk7XHJcbiAgICBpZiAoblJvd0dyaWQgPj0gMCkge1xyXG4gICAgICB0aGlzLm1fblJlc2l6ZVJvd0dyaWRNb3ZpbmcgPSBuUm93R3JpZDtcclxuICAgICAgZG9jdW1lbnQuYm9keS5zdHlsZS5jdXJzb3IgPSBcInJvdy1yZXNpemVcIjtcclxuICAgICAgcmV0dXJuO1xyXG4gICAgfVxyXG5cclxuICAgIGlmKHRoaXMubV9uUmVzaXplUm93R3JpZE1vdmluZyA+PSAwKSB7XHJcbiAgICAgIHRoaXMubV9uUmVzaXplUm93R3JpZE1vdmluZyA9IC0xO1xyXG4gICAgICBkb2N1bWVudC5ib2R5LnN0eWxlLmN1cnNvciA9IFwiYXV0b1wiO1xyXG4gICAgfVxyXG5cclxuXHJcbiAgICAvL0hhbWJ1cmdlciBNZW51XHJcbiAgICBjb25zdCBjb2xHcmlkID0gdGhpcy5nZXRHcmlkQ29sdW1uKCk7XHJcbiAgICBpZihjb2xHcmlkID09PSBudWxsIHx8IGNvbEdyaWQubmFtZSA9PT0gJycpXHJcbiAgICAgIHJldHVybjtcclxuXHJcbiAgICBjb25zdCBlRGl2SGFtYiA9IEdyaWRVdGlscy5nZXRUb29sSWNvbkRpdihjb2xHcmlkLmdyaWQpO1xyXG4gICAgY29uc3QgbkhDb2xIZWFkZXIgPSBHcmlkVXRpbHMuZ2V0R3JpZENvbHVtbkhlYWRlckhlaWdodChjb2xHcmlkLmdyaWQpO1xyXG4gICAgaWYoMCA8PSBlLm9mZnNldFkgJiYgZS5vZmZzZXRZIDwgbkhDb2xIZWFkZXIpIHtcclxuXHJcbiAgICAgIGVEaXZIYW1iPy5zdHlsZS5yZW1vdmVQcm9wZXJ0eSgndmlzaWJpbGl0eScpO1xyXG4gICAgICBlRGl2SGFtYj8uc2V0QXR0cmlidXRlKCdjb2x1bW5fbmFtZScsIGNvbEdyaWQubmFtZSk7XHJcbiAgICAgIC8vY29uc29sZS5sb2coJ1Rvb2xzSWNvbiBmb3IgY29sdW1uICcgKyBjb2xHcmlkLm5hbWUpO1xyXG4gICAgICAvLyBAdHMtaWdub3JlXHJcbiAgICAgIGVEaXZIYW1iPy5zdHlsZS5sZWZ0ID0gKFBpbm5lZFV0aWxzLmdldFBpbm5lZENvbHVtbkxlZnQodGhpcykgKyB0aGlzLmdldFdpZHRoKCkgLSAxOCkgKyAncHgnO1xyXG4gICAgICAvLyBAdHMtaWdub3JlXHJcbiAgICAgIGVEaXZIYW1iPy5zdHlsZS50b3AgPSAoR3JpZFV0aWxzLmdldEdyaWRDb2x1bW5IZWFkZXJIZWlnaHQoY29sR3JpZC5ncmlkKSAtIDE2KSArIFwicHhcIjtcclxuICAgIH0gZWxzZSB7XHJcbiAgICAgIGNvbnN0IGNvbEdyaWQgPSB0aGlzLmdldEdyaWRDb2x1bW4oKTtcclxuICAgICAgaWYoY29sR3JpZCAhPSBudWxsKSB7XHJcbiAgICAgICAgICBlRGl2SGFtYj8uc2V0QXR0cmlidXRlKCdjb2x1bW5fbmFtZScsICcnKTtcclxuICAgICAgICAgIC8vIEB0cy1pZ25vcmVcclxuICAgICAgICAgIGVEaXZIYW1iPy5zdHlsZS52aXNpYmlsaXR5ID0gJ2hpZGRlbic7XHJcbiAgICAgICAgfVxyXG4gICAgfVxyXG4gIH1cclxuXHJcbiAgcHVibGljIG9uTW91c2VEcmFnKGUgOiBNb3VzZUV2ZW50KSA6IHZvaWQge1xyXG4gICAgaWYoREVCVUcpXHJcbiAgICAgY29uc29sZS5sb2coJ01vdXNlIERyYWcgUGlubmVkIENvbHVtbjogJyArIHRoaXMuZ2V0R3JpZENvbHVtbigpPy5uYW1lKTtcclxuXHJcbiAgICBpZih0aGlzLm1fY29sR3JpZCA9PT0gbnVsbCB8fCB0aGlzLm1fcm9vdCA9PT0gbnVsbClcclxuICAgIHJldHVybjtcclxuXHJcbiAgICBjb25zdCBncmlkID0gdGhpcy5tX2NvbEdyaWQuZ3JpZDtcclxuICAgIGNvbnN0IHZpZXdUYWJsZSA9IGdyaWQudmlldztcclxuXHJcbiAgICBpZihERy50b0RhcnQoZ3Jvay5zaGVsbC52KSAhPT0gREcudG9EYXJ0KHZpZXdUYWJsZSkpIHtcclxuICAgICAgcmV0dXJuO1xyXG4gICAgfVxyXG5cclxuICAgIGNvbnN0IGJSZXNpemluZyA9IHRoaXMubV9uUmVzaXplUm93R3JpZERyYWdnaW5nID49IDA7XHJcbiAgICBpZiAoYlJlc2l6aW5nKSB7XHJcblxyXG4gICAgICAvL2NvbnNvbGUubG9nKFwiRHJhZ2dpbmcgOiBcIiArIGhlYWRlclRoaXMubV9zdHJDb2xOYW1lKTtcclxuICAgICAgY29uc3QgbllEaWZmID0gZS5jbGllbnRZIC0gdGhpcy5tX25ZUmVzaXplRHJhZ2dpbmdBbmNob3I7XHJcbiAgICAgIGxldCBuSFJvd0dyaWQgPSB0aGlzLm1fbkhSZXNpemVSb3dzQmVmb3JlRHJhZyArIG5ZRGlmZjtcclxuXHJcbiAgICAgIGlmIChuSFJvd0dyaWQgPCBQaW5uZWRDb2x1bW4uTUlOX1JPV19IRUlHSFQpXHJcbiAgICAgICAgbkhSb3dHcmlkID0gUGlubmVkQ29sdW1uLk1JTl9ST1dfSEVJR0hUO1xyXG4gICAgICBlbHNlIGlmIChuSFJvd0dyaWQgPiBQaW5uZWRDb2x1bW4uTUFYX1JPV19IRUlHSFQpXHJcbiAgICAgICAgbkhSb3dHcmlkID0gUGlubmVkQ29sdW1uLk1BWF9ST1dfSEVJR0hUO1xyXG5cclxuICAgICAgY29uc3QgZUNhbnZhc1RoaXMgPSB0aGlzLm1fcm9vdDtcclxuXHJcbiAgICAgIGxldCBnID0gZUNhbnZhc1RoaXMuZ2V0Q29udGV4dCgnMmQnKTtcclxuICAgICAgaWYoZyA9PT0gbnVsbClcclxuICAgICAgICByZXR1cm47XHJcblxyXG4gICAgICBnLmZpbGxTdHlsZSA9IFwid2hpdGVcIjtcclxuICAgICAgY29uc3QgbkhIZWFkZXJDb2xzID0gR3JpZFV0aWxzLmdldEdyaWRDb2x1bW5IZWFkZXJIZWlnaHQoZ3JpZCk7XHJcbiAgICAgIGcuZmlsbFJlY3QoMCxuSEhlYWRlckNvbHMsIGVDYW52YXNUaGlzLm9mZnNldFdpZHRoLCBlQ2FudmFzVGhpcy5vZmZzZXRIZWlnaHQpO1xyXG5cclxuICAgICAgZ3JpZC5zZXRPcHRpb25zKHtcclxuICAgICAgICByb3dIZWlnaHQ6IG5IUm93R3JpZCAvL3RoaXMgd29uJ3QgdHJpZ2dlciBvblJvd3NSZXppemVkIGV2ZW50LCB3aGljaCBpcyBhIERHIGJ1Z1xyXG4gICAgICB9KTtcclxuXHJcbiAgICAgIG5vdGlmeUFsbFBpbm5lZENvbHNSb3dzUmVzaXplZCh0aGlzLCBuSFJvd0dyaWQsIHRydWUpO1xyXG4gICAgICBub3RpZnlBbGxDb2xzUm93c1Jlc2l6ZWQoZ3JpZCwgbkhSb3dHcmlkLCB0cnVlKTtcclxuXHJcbiAgICAgIGxldCBoZWFkZXIgPSBudWxsO1xyXG4gICAgICBjb25zdCBhciA9IGdyaWQuZGFydC5tX2FyUGlubmVkQ29scztcclxuICAgICAgZm9yKGxldCBuPTA7IG48YXIubGVuZ3RoOyArK24pIHtcclxuICAgICAgICBoZWFkZXIgPSBhcltuXTtcclxuICAgICAgICBnID0gaGVhZGVyLm1fcm9vdC5nZXRDb250ZXh0KCcyZCcpO1xyXG4gICAgICAgIGhlYWRlci5wYWludChnLCBncmlkKTtcclxuICAgICAgfVxyXG5cclxuICAgICAgdHJ5IHtcclxuICAgICAgICBjb25zdCBjb2xHcmlkMCA9IGdyaWQuY29sdW1ucy5ieUluZGV4KDApO1xyXG4gICAgICAgIGlmIChjb2xHcmlkMCAhPT0gbnVsbClcclxuICAgICAgICAgIGNvbEdyaWQwLnZpc2libGUgPSBmYWxzZTsvL3RlbXBvcmFyeSBhZGRyZXNzZWQgdGhlIERHIGJ1Z1xyXG4gICAgICB9XHJcbiAgICAgIGNhdGNoKGUpIHtcclxuICAgICAgICAvL0RHIGJ1Z1xyXG4gICAgICB9XHJcbiAgICAgIHJldHVybjtcclxuICAgIH1cclxuXHJcblxyXG4gIH1cclxuXHJcbiAgcHVibGljIG9uTW91c2VMZWF2ZShlIDogTW91c2VFdmVudCwgYk92ZXJsYXAgOiBib29sZWFuKSA6IHZvaWQge1xyXG4gICAgaWYoREVCVUcpXHJcbiAgICAgY29uc29sZS5sb2coJ01vdXNlIExlZnQgUGlubmVkIENvbHVtbjogJyArIHRoaXMuZ2V0R3JpZENvbHVtbigpPy5uYW1lICsgJyAgb3ZlcmxhcDogJyArIGJPdmVybGFwKTtcclxuXHJcbiAgICBpZih0aGlzLm1fblJlc2l6ZVJvd0dyaWRNb3ZpbmcgPj0gMCkge1xyXG4gICAgICB0aGlzLm1fblJlc2l6ZVJvd0dyaWRNb3ZpbmcgPSAtMTtcclxuICAgICAgZG9jdW1lbnQuYm9keS5zdHlsZS5jdXJzb3IgPSBcImF1dG9cIjtcclxuICAgIH1cclxuXHJcbiAgICBpZih0aGlzLm1fY2VsbEN1cnJlbnQgIT09IG51bGwpIHtcclxuICAgICAgY29uc3QgcmVuZGVyZXIgPSBnZXRSZW5kZXJlcih0aGlzLm1fY2VsbEN1cnJlbnQpO1xyXG4gICAgICBpZiAocmVuZGVyZXIgaW5zdGFuY2VvZiBHcmlkQ2VsbFJlbmRlcmVyRXgpIHtcclxuICAgICAgICBjb25zdCBlTW91c2UgPSBlIGFzIE1vdXNlRXZlbnQ7XHJcbiAgICAgICAgcmVuZGVyZXIub25Nb3VzZUxlYXZlRXgodGhpcy5tX2NlbGxDdXJyZW50LCBlTW91c2UsIC0xLCAtMSk7XHJcbiAgICAgIH1cclxuICAgICAgdGhpcy5tX2NlbGxDdXJyZW50ID0gbnVsbDtcclxuICAgIH1cclxuXHJcbiAgICBjb25zdCBjb2xHcmlkID0gdGhpcy5nZXRHcmlkQ29sdW1uKCk7XHJcbiAgICBpZihjb2xHcmlkICE9IG51bGwgJiYgIWJPdmVybGFwKSB7XHJcbiAgICAgIGNvbnN0IGVEaXZIYW1iID0gR3JpZFV0aWxzLmdldFRvb2xJY29uRGl2KGNvbEdyaWQuZ3JpZCk7XHJcbiAgICAgIGVEaXZIYW1iPy5zZXRBdHRyaWJ1dGUoJ2NvbHVtbl9uYW1lJywgJycpO1xyXG4gICAgICAvLyBAdHMtaWdub3JlXHJcbiAgICAgIGVEaXZIYW1iPy5zdHlsZS52aXNpYmlsaXR5ID0gJ2hpZGRlbic7XHJcbiAgICB9XHJcblxyXG5cclxuICB9XHJcblxyXG4gIHB1YmxpYyBvbk1vdXNlRGJsQ2xpY2soZSA6IE1vdXNlRXZlbnQpIDogdm9pZCB7XHJcbiAgICBpZihERUJVRylcclxuICAgICBjb25zb2xlLmxvZygnTW91c2UgRGJsIENsaWNrZWQgUGlubmVkIENvbHVtbjogJyArIHRoaXMuZ2V0R3JpZENvbHVtbigpPy5uYW1lKTtcclxuXHJcbiAgICBpZih0aGlzLm1fY29sR3JpZCA9PT0gbnVsbCB8fCB0aGlzLm1fcm9vdCA9PT0gbnVsbClcclxuICAgICAgcmV0dXJuO1xyXG5cclxuICAgIGNvbnN0IGdyaWQgPSB0aGlzLm1fY29sR3JpZC5ncmlkO1xyXG4gICAgY29uc3Qgdmlld1RhYmxlID0gZ3JpZD8udmlldztcclxuXHJcbiAgICBpZiAoREcudG9EYXJ0KGdyb2suc2hlbGwudikgIT09IERHLnRvRGFydCh2aWV3VGFibGUpKSB7XHJcbiAgICAgIHJldHVybjtcclxuICAgIH1cclxuXHJcbiAgICBpZih0aGlzLm1fY29sR3JpZD8ubmFtZSA9PT0gJycpXHJcbiAgICAgIHJldHVybjtcclxuXHJcbiAgICBpZih0aGlzLm1fYlNvcnRlZEFzY2VuZGluZyA9PSBudWxsKVxyXG4gICAgICB0aGlzLm1fYlNvcnRlZEFzY2VuZGluZyA9IHRydWU7XHJcbiAgICBlbHNlIGlmKHRoaXMubV9iU29ydGVkQXNjZW5kaW5nKVxyXG4gICAgICB0aGlzLm1fYlNvcnRlZEFzY2VuZGluZyA9IGZhbHNlO1xyXG4gICAgZWxzZSB0aGlzLm1fYlNvcnRlZEFzY2VuZGluZyA9IHRydWU7XHJcblxyXG4gICAgY29uc3QgbkhIZWFkZXJDb2xzID0gR3JpZFV0aWxzLmdldEdyaWRDb2x1bW5IZWFkZXJIZWlnaHQoZ3JpZCk7XHJcblxyXG4gICAgaWYoMCA8PSBlLm9mZnNldFggJiYgZS5vZmZzZXRYIDw9IHRoaXMubV9yb290Lm9mZnNldFdpZHRoICYmXHJcbiAgICAgICAgMCA8PSBlLm9mZnNldFkgJiYgZS5vZmZzZXRZIDw9IG5ISGVhZGVyQ29scykgICAvL29uIHRoZSByb3dzIGhlYWRlclxyXG4gICAge1xyXG4gICAgICBncmlkPy5zb3J0KFt0aGlzLm1fY29sR3JpZD8ubmFtZV0sIFt0aGlzLm1fYlNvcnRlZEFzY2VuZGluZ10pO1xyXG4gICAgfVxyXG4gIH1cclxuXHJcbiAgcHVibGljIG9uTW91c2VEb3duKGUgOiBNb3VzZUV2ZW50KSA6IHZvaWQge1xyXG4gICAgaWYoREVCVUcpXHJcbiAgICAgY29uc29sZS5sb2coJ01vdXNlIERvd24gUGlubmVkIENvbHVtbjogJyArIHRoaXMuZ2V0R3JpZENvbHVtbigpPy5uYW1lKTtcclxuLypcclxuICAgIGlmKGUudmlldyAhPSBudWxsKSB7XHJcbiAgICAgIGNvbnN0IGVlID0gZG9jdW1lbnQuY3JlYXRlRXZlbnQoIFwiTW91c2VFdmVudFwiICk7XHJcbiAgICAgIGVlLmluaXRNb3VzZUV2ZW50KGUudHlwZSwgZS5idWJibGVzLCBlLmNhbmNlbGFibGUsIGUudmlldywgZS5kZXRhaWwsIGUuc2NyZWVuWCArIDEwMCwgZS5zY3JlZW5ZLCBlLmNsaWVudFggKyAxMDAsIGUuY2xpZW50WSwgZS5jdHJsS2V5LCBlLmFsdEtleSwgZS5zaGlmdEtleSwgZS5tZXRhS2V5LCBlLmJ1dHRvbiwgZS5yZWxhdGVkVGFyZ2V0KTtcclxuICAgICAgdGhpcy5tX2NvbEdyaWQ/LmdyaWQucm9vdC5kaXNwYXRjaEV2ZW50KGVlKTtcclxuICAgICAgcmV0dXJuO1xyXG4gICAgfVxyXG4qL1xyXG5cclxuICAgIGlmKHRoaXMubV9jb2xHcmlkID09PSBudWxsKVxyXG4gICAgICByZXR1cm47XHJcblxyXG4gICAgY29uc3QgZ3JpZCA9IHRoaXMubV9jb2xHcmlkPy5ncmlkO1xyXG4gICAgY29uc3Qgdmlld1RhYmxlID0gZ3JpZD8udmlldztcclxuICAgIGlmKERHLnRvRGFydChncm9rLnNoZWxsLnYpICE9PSBERy50b0RhcnQodmlld1RhYmxlKSlcclxuICAgICAgcmV0dXJuO1xyXG5cclxuICAgIGlmKGUuYnV0dG9ucyAhPT0gMSlcclxuICAgICAgcmV0dXJuO1xyXG5cclxuICAgIGxldCBlQ2FudmFzVGhpcyA9IHRoaXMubV9yb290O1xyXG4gICAgaWYoZUNhbnZhc1RoaXMgPT09IG51bGwpXHJcbiAgICAgIHJldHVybjtcclxuXHJcbiAgICB0aGlzLm1fblJlc2l6ZVJvd0dyaWRNb3ZpbmcgPSAtMTtcclxuICAgIGNvbnN0IGJBZGRUb1NlbCA6IGJvb2xlYW4gPSBlLmN0cmxLZXkgfHwgZS5zaGlmdEtleTtcclxuXHJcbiAgICBsZXQgblJvd0dyaWQgPSBiQWRkVG9TZWwgPyAtMSA6IFBpbm5lZENvbHVtbi5oaXRUZXN0Um93cyhlQ2FudmFzVGhpcywgZ3JpZCwgZSwgdHJ1ZSwgdW5kZWZpbmVkKTtcclxuICAgIGlmIChuUm93R3JpZCA+PSAwKSB7XHJcbiAgICAgIGNvbnN0IG5IUm93cyA9IEdyaWRVdGlscy5nZXRHcmlkUm93SGVpZ2h0KGdyaWQpO1xyXG4gICAgICB0aGlzLm1fblJlc2l6ZVJvd0dyaWREcmFnZ2luZyA9IG5Sb3dHcmlkO1xyXG4gICAgICB0aGlzLm1fbllSZXNpemVEcmFnZ2luZ0FuY2hvciA9IGUuY2xpZW50WTtcclxuICAgICAgdGhpcy5tX25IUmVzaXplUm93c0JlZm9yZURyYWcgPSBuSFJvd3M7XHJcbiAgICB9XHJcbiAgICBlbHNlXHJcbiAgICB7XHJcblxyXG4gICAgICBuUm93R3JpZCA9IFBpbm5lZENvbHVtbi5oaXRUZXN0Um93cyhlQ2FudmFzVGhpcywgZ3JpZCwgZSwgZmFsc2UsIHRoaXMubV9hclhZTW91c2VPbkNlbGxEb3duKTtcclxuXHJcbiAgICAgIHRoaXMubV9uUm93R3JpZERyYWdnaW5nID0gblJvd0dyaWQ7XHJcbiAgICAgIHRoaXMubV9uWURyYWdnaW5nQW5jaG9yID0gZS5jbGllbnRZO1xyXG5cclxuICAgICAgY29uc3QgY2VsbCA9IGdyaWQuY2VsbCh0aGlzLm1fY29sR3JpZC5uYW1lLCBuUm93R3JpZCk7XHJcbiAgICAgIGNvbnN0IHJlbmRlcmVyID0gZ2V0UmVuZGVyZXIoY2VsbCk7XHJcbiAgICAgIGlmKHJlbmRlcmVyIGluc3RhbmNlb2YgR3JpZENlbGxSZW5kZXJlckV4KSB7XHJcbiAgICAgICAgcmVuZGVyZXIub25Nb3VzZURvd25FeChjZWxsLCBlLCB0aGlzLm1fYXJYWU1vdXNlT25DZWxsRG93blswXSwgdGhpcy5tX2FyWFlNb3VzZU9uQ2VsbERvd25bMV0pO1xyXG4gICAgICB9XHJcbiAgICB9XHJcbiAgfVxyXG5cclxuICBwdWJsaWMgb25Nb3VzZVVwKGUgOiBNb3VzZUV2ZW50KSA6IHZvaWQge1xyXG4gICAgaWYoREVCVUcpXHJcbiAgICAgY29uc29sZS5sb2coJ01vdXNlIFVwIFBpbm5lZCBDb2x1bW46ICcgKyB0aGlzLmdldEdyaWRDb2x1bW4oKT8ubmFtZSk7XHJcbi8qXHJcbiAgICBpZihlLnZpZXcgIT0gbnVsbCkge1xyXG4gICAgICBjb25zdCBlZSA9IGRvY3VtZW50LmNyZWF0ZUV2ZW50KCBcIk1vdXNlRXZlbnRcIiApO1xyXG4gICAgICBlZS5pbml0TW91c2VFdmVudChlLnR5cGUsIGUuYnViYmxlcywgZS5jYW5jZWxhYmxlLCBlLnZpZXcsIGUuZGV0YWlsLCBlLnNjcmVlblggKyAxMDAsIGUuc2NyZWVuWSwgZS5jbGllbnRYICsgMTAwLCBlLmNsaWVudFksIGUuY3RybEtleSwgZS5hbHRLZXksIGUuc2hpZnRLZXksIGUubWV0YUtleSwgZS5idXR0b24sIGUucmVsYXRlZFRhcmdldCk7XHJcbiAgICAgIHRoaXMubV9jb2xHcmlkPy5ncmlkLnJvb3QuZGlzcGF0Y2hFdmVudChlZSk7XHJcbiAgICAgIHJldHVybjtcclxuICAgIH1cclxuKi9cclxuXHJcbiAgICBpZih0aGlzLm1fY29sR3JpZCA9PT0gbnVsbCB8fCB0aGlzLm1fcm9vdCA9PSBudWxsKVxyXG4gICAgICByZXR1cm47XHJcblxyXG4gICAgY29uc3QgZ3JpZCA9IHRoaXMubV9jb2xHcmlkPy5ncmlkO1xyXG4gICAgY29uc3Qgdmlld1RhYmxlID0gZ3JpZD8udmlldztcclxuXHJcbiAgICBpZihERy50b0RhcnQoZ3Jvay5zaGVsbC52KSAhPT0gREcudG9EYXJ0KHZpZXdUYWJsZSkpIHtcclxuICAgICAgcmV0dXJuO1xyXG4gICAgfVxyXG5cclxuICAgIGlmKHRoaXMubV9uUmVzaXplUm93R3JpZERyYWdnaW5nID49IDApIHtcclxuICAgICAgY29uc3QgbkhSb3cgPSBHcmlkVXRpbHMuZ2V0R3JpZFJvd0hlaWdodChncmlkKTtcclxuICAgICAgbm90aWZ5QWxsUGlubmVkQ29sc1Jvd3NSZXNpemVkKHRoaXMsIG5IUm93LCBmYWxzZSk7XHJcbiAgICAgIG5vdGlmeUFsbENvbHNSb3dzUmVzaXplZChncmlkLCBuSFJvdywgZmFsc2UpO1xyXG4gICAgfVxyXG5cclxuICAgIHRoaXMubV9uSFJlc2l6ZVJvd3NCZWZvcmVEcmFnID0gLTE7XHJcbiAgICB0aGlzLm1fblJlc2l6ZVJvd0dyaWREcmFnZ2luZyA9IC0xO1xyXG4gICAgdGhpcy5tX25ZUmVzaXplRHJhZ2dpbmdBbmNob3IgPSAtMTtcclxuICAgIHRoaXMubV9uUmVzaXplUm93R3JpZE1vdmluZyA9IC0xO1xyXG5cclxuICAgIGRvY3VtZW50LmJvZHkuc3R5bGUuY3Vyc29yID0gXCJhdXRvXCI7XHJcblxyXG4gICAgaWYodGhpcy5tX25Sb3dHcmlkRHJhZ2dpbmcgPj0gMCkge1xyXG4gICAgICBjb25zdCBkZnJhbWUgPSBncmlkLmRhdGFGcmFtZTtcclxuICAgICAgY29uc3QgYkN0cmwgPSBlLmN0cmxLZXk7XHJcbiAgICAgIGNvbnN0IGJSYW5nZVNlbCA9IGUuc2hpZnRLZXk7XHJcblxyXG4gICAgICBsZXQgYlNlbCA9IHRydWU7XHJcblxyXG4gICAgICBjb25zdCBuUm93R3JpZCA9IFBpbm5lZENvbHVtbi5oaXRUZXN0Um93cyh0aGlzLm1fcm9vdCwgZ3JpZCwgZSwgZmFsc2UsIHRoaXMubV9hclhZTW91c2VPbkNlbGxVcCk7XHJcbiAgICAgIGlmKCFiQ3RybCAmJiAhYlJhbmdlU2VsICYmIG5Sb3dHcmlkID09PSB0aGlzLm1fblJvd0dyaWREcmFnZ2luZykgeyAvL2NsaWNrIG9uIHRoZSBzYW1lIHJvdyB3aGljaCB3aWxsIGJlY29tZSBhY3RpdmVcclxuXHJcbiAgICAgICAgbGV0IGNlbGxSSCA9IG51bGw7XHJcbiAgICAgICAgdHJ5IHtcclxuICAgICAgICAgIGNlbGxSSCA9IGdyaWQuY2VsbChcIlwiLCBuUm93R3JpZCk7XHJcbiAgICAgICAgfVxyXG4gICAgICAgIGNhdGNoKGUpIHtcclxuICAgICAgICAgIGxldCBjb2xHID0gbnVsbDtcclxuICAgICAgICAgIGNvbnN0IGxzdENvbHMgPSBncmlkLmNvbHVtbnM7XHJcbiAgICAgICAgICBmb3IobGV0IG5DPTE7IG5DPGxzdENvbHMubGVuZ3RoOyArK25DKSB7XHJcbiAgICAgICAgICAgIGNvbEcgPSBsc3RDb2xzLmJ5SW5kZXgobkMpO1xyXG4gICAgICAgICAgICBjZWxsUkggPSBjb2xHID09PSBudWxsID8gbnVsbCA6IGdyaWQuY2VsbChjb2xHLm5hbWUsIG5Sb3dHcmlkKTtcclxuICAgICAgICAgICAgaWYoY2VsbFJIICE9PSBudWxsKVxyXG4gICAgICAgICAgICAgIGJyZWFrO1xyXG4gICAgICAgICAgfVxyXG4gICAgICAgIH1cclxuICAgICAgICBpZihjZWxsUkggIT09IG51bGwpIHtcclxuICAgICAgICAgIGNvbnN0IG5Sb3dUYWJsZSA6IGFueSA9IGNlbGxSSC50YWJsZVJvd0luZGV4O1xyXG4gICAgICAgICAgaWYoblJvd1RhYmxlICE9PSBudWxsKVxyXG4gICAgICAgICAgICBkZnJhbWUuY3VycmVudFJvdyA9IG5Sb3dUYWJsZTtcclxuICAgICAgICB9XHJcbiAgICAgIH1cclxuICAgICAgZWxzZVxyXG4gICAgICB7XHJcbiAgICAgICAgY29uc3QgYml0c2V0U2VsID0gZGZyYW1lLnNlbGVjdGlvbjtcclxuXHJcblxyXG5cclxuICAgICAgICBsZXQgblJvd01pbiA9IHRoaXMubV9uUm93R3JpZERyYWdnaW5nIDwgblJvd0dyaWQgPyB0aGlzLm1fblJvd0dyaWREcmFnZ2luZyA6IG5Sb3dHcmlkO1xyXG4gICAgICAgIGxldCBuUm93TWF4ID0gdGhpcy5tX25Sb3dHcmlkRHJhZ2dpbmcgPiBuUm93R3JpZCA/IHRoaXMubV9uUm93R3JpZERyYWdnaW5nIDogblJvd0dyaWQ7XHJcblxyXG4gICAgICAgIGlmKGJDdHJsKSB7XHJcbiAgICAgICAgICBuUm93TWluID0gblJvd0dyaWQ7XHJcbiAgICAgICAgICBuUm93TWF4ID0gblJvd0dyaWQ7XHJcbiAgICAgICAgICBsZXQgYkN1clNlbCA9IGJpdHNldFNlbC5nZXQoblJvd0dyaWQpO1xyXG4gICAgICAgICAgYlNlbCA9ICFiQ3VyU2VsO1xyXG4gICAgICAgIH1cclxuICAgICAgICBlbHNlIGlmKGJSYW5nZVNlbCkge1xyXG4gICAgICAgICAgbGV0IG5Sb3dHcmlkQWN0aXZlID0gR3JpZFV0aWxzLmdldEFjdGl2ZUdyaWRSb3coZ3JpZCk7XHJcbiAgICAgICAgICBpZihuUm93R3JpZEFjdGl2ZSA9PT0gbnVsbClcclxuICAgICAgICAgICAgblJvd0dyaWRBY3RpdmUgPSAwO1xyXG5cclxuICAgICAgICAgIGlmKG5Sb3dNaW4gPT0gblJvd01heCkge1xyXG4gICAgICAgICAgICBiaXRzZXRTZWwuc2V0QWxsKGZhbHNlLCB0cnVlKTtcclxuXHJcbiAgICAgICAgICAgIG5Sb3dNaW4gPSBuUm93R3JpZEFjdGl2ZSA8IG5Sb3dHcmlkID8gblJvd0dyaWRBY3RpdmUgOiBuUm93R3JpZDtcclxuICAgICAgICAgICAgblJvd01heCA9IG5Sb3dHcmlkQWN0aXZlID4gblJvd0dyaWQgPyBuUm93R3JpZEFjdGl2ZSA6IG5Sb3dHcmlkO1xyXG4gICAgICAgICAgfVxyXG4gICAgICAgIH1cclxuXHJcblxyXG4gICAgICAgIC8vaWYoIWJDdHJsIHx8IGJSYW5nZVNlbClcclxuICAgICAgICAgIC8vYml0c2V0U2VsLnNldEFsbChmYWxzZSwgdHJ1ZSk7XHJcblxyXG4gICAgICAgIGxldCBjZWxsUkggPSBudWxsO1xyXG4gICAgICAgIGxldCBuUm93VGFibGUgPSAtMTtcclxuICAgICAgICBmb3IobGV0IG5Sb3c9blJvd01pbjsgblJvdzw9blJvd01heDsgKytuUm93KSB7XHJcblxyXG4gICAgICAgICAgdHJ5IHtcclxuICAgICAgICAgICAgY2VsbFJIID0gZ3JpZC5jZWxsKFwiXCIsIG5Sb3cpO1xyXG4gICAgICAgICAgfVxyXG4gICAgICAgICAgY2F0Y2goZSkge1xyXG4gICAgICAgICAgICBsZXQgY29sRyA9IG51bGw7XHJcbiAgICAgICAgICAgIGNvbnN0IGxzdENvbHMgPSBncmlkLmNvbHVtbnM7XHJcbiAgICAgICAgICAgIGZvcihsZXQgbkM9MTsgbkM8bHN0Q29scy5sZW5ndGg7ICsrbkMpIHtcclxuICAgICAgICAgICAgICBjb2xHID0gbHN0Q29scy5ieUluZGV4KG5DKTtcclxuICAgICAgICAgICAgICBjZWxsUkggPSBjb2xHID09PSBudWxsID8gbnVsbCA6IGdyaWQuY2VsbChjb2xHLm5hbWUsIG5Sb3dHcmlkKTtcclxuICAgICAgICAgICAgICBpZihjZWxsUkggIT09IG51bGwpXHJcbiAgICAgICAgICAgICAgICBicmVhaztcclxuICAgICAgICAgICAgfVxyXG4gICAgICAgICAgfVxyXG5cclxuICAgICAgICAgIGlmKGNlbGxSSCAhPT0gbnVsbCAmJiBjZWxsUkgudGFibGVSb3dJbmRleCAhPT0gbnVsbCkge1xyXG4gICAgICAgICAgICBuUm93VGFibGUgPSBjZWxsUkgudGFibGVSb3dJbmRleDtcclxuICAgICAgICAgICAgYml0c2V0U2VsLnNldChuUm93VGFibGUsIGJTZWwsIHRydWUpO1xyXG4gICAgICAgICAgfVxyXG4gICAgICAgIH1cclxuICAgICAgfVxyXG5cclxuICAgICAgY29uc3QgY2VsbCA9IGdyaWQuY2VsbCh0aGlzLm1fY29sR3JpZC5uYW1lLCBuUm93R3JpZCk7XHJcbiAgICAgIGNvbnN0IHJlbmRlcmVyID0gZ2V0UmVuZGVyZXIoY2VsbCk7XHJcbiAgICAgIGlmKHJlbmRlcmVyIGluc3RhbmNlb2YgR3JpZENlbGxSZW5kZXJlckV4KSB7XHJcbiAgICAgICAgcmVuZGVyZXIub25Nb3VzZVVwRXgoY2VsbCwgZSwgdGhpcy5tX2FyWFlNb3VzZU9uQ2VsbFVwWzBdLCB0aGlzLm1fYXJYWU1vdXNlT25DZWxsVXBbMV0pO1xyXG4gICAgICB9XHJcblxyXG4gICAgICBpZih0aGlzLm1fYXJYWU1vdXNlT25DZWxsVXBbMF0gPT09IHRoaXMubV9hclhZTW91c2VPbkNlbGxEb3duWzBdICYmIHRoaXMubV9hclhZTW91c2VPbkNlbGxEb3duWzFdID09PSB0aGlzLm1fYXJYWU1vdXNlT25DZWxsVXBbMV0pIHtcclxuICAgICAgICBpZihyZW5kZXJlciBpbnN0YW5jZW9mIEdyaWRDZWxsUmVuZGVyZXJFeCkge1xyXG4gICAgICAgICAgcmVuZGVyZXIub25DbGlja0V4KGNlbGwsIGUsIHRoaXMubV9hclhZTW91c2VPbkNlbGxVcFswXSwgdGhpcy5tX2FyWFlNb3VzZU9uQ2VsbFVwWzFdKTtcclxuICAgICAgICB9XHJcbiAgICAgIH1cclxuXHJcbiAgICAgIHRoaXMubV9uUm93R3JpZERyYWdnaW5nID0gLTE7XHJcbiAgICAgIHRoaXMubV9uWURyYWdnaW5nQW5jaG9yID0gLTE7XHJcbiAgICAgIHRoaXMubV9hclhZTW91c2VPbkNlbGxEb3duWzBdID0gLTI7XHJcbiAgICAgIHRoaXMubV9hclhZTW91c2VPbkNlbGxEb3duWzFdID0gLTI7XHJcbiAgICAgIHRoaXMubV9hclhZTW91c2VPbkNlbGxVcFswXSA9IC0xO1xyXG4gICAgICB0aGlzLm1fYXJYWU1vdXNlT25DZWxsVXBbMV0gPSAtMTtcclxuICAgIH1cclxuICB9XHJcblxyXG4gIHB1YmxpYyBvbkNvbnRleHRNZW51KGUgOiBNb3VzZUV2ZW50KSA6IHZvaWQge1xyXG4gICBpZihERUJVRylcclxuICAgIGNvbnNvbGUubG9nKCdDb250ZXh0IG1lbnUgUGlubmVkIENvbHVtbjogJyArIHRoaXMuZ2V0R3JpZENvbHVtbigpPy5uYW1lKTtcclxuICB9XHJcblxyXG4gIHB1YmxpYyBvbk1vdXNlV2hlZWwoZSA6IFdoZWVsRXZlbnQpIDogdm9pZCB7XHJcblxyXG4gICAgaWYodGhpcy5tX2NvbEdyaWQgPT09IG51bGwgfHwgdGhpcy5tX3Jvb3QgPT0gbnVsbClcclxuICAgICAgcmV0dXJuO1xyXG5cclxuICAgIGNvbnN0IGdyaWQgPSB0aGlzLm1fY29sR3JpZD8uZ3JpZDtcclxuICAgIGNvbnN0IHZpZXdUYWJsZSA9IGdyaWQ/LnZpZXc7XHJcblxyXG4gICAgaWYgKERHLnRvRGFydChncm9rLnNoZWxsLnYpICE9PSBERy50b0RhcnQodmlld1RhYmxlKSkge1xyXG4gICAgICByZXR1cm47XHJcbiAgICB9XHJcblxyXG4gICAgaWYoZS5kZWx0YVggIT09IDAgfHwgZS5kZWx0YVogIT09IDApIHtcclxuICAgICAgcmV0dXJuO1xyXG4gICAgfVxyXG5cclxuICAgIHNldFRpbWVvdXQoKCkgPT57XHJcbiAgICAgIGNvbnN0IGVlID0gbmV3IFdoZWVsRXZlbnQoZS50eXBlLCBlKTtcclxuICAgICAgdHJ5e2dyaWQub3ZlcmxheS5kaXNwYXRjaEV2ZW50KGVlKTt9XHJcbiAgICAgIGNhdGNoKGV4KSB7XHJcbiAgICAgICAgLy9jb25zb2xlLmVycm9yKGV4Lm1lc3NhZ2UpO1xyXG4gICAgICB9XHJcbiAgICB9LCAxKTtcclxuXHJcblxyXG4gICAgaWYodHJ1ZSlcclxuICAgICAgcmV0dXJuO1xyXG4gICAgLy9lLmNsaWVudFggPSA1O1xyXG5cclxuXHJcbiAgICBpZih0aGlzLm1fbldoZWVsQ291bnQgPT09IDEpIHtcclxuICAgICAgLy9zY3JvbGwgK1xyXG4gICAgICBjb25zdCBuUm93Q291bnQgPSBHcmlkVXRpbHMuZ2V0R3JpZFZpc2libGVSb3dDb3VudChncmlkKTtcclxuICAgICAgY29uc3Qgc2Nyb2xsWSA9IGdyaWQudmVydFNjcm9sbDtcclxuICAgICAgaWYoblJvd0NvdW50IC0xID4gc2Nyb2xsWS5tYXgpIHtcclxuICAgICAgICBzY3JvbGxZLnNldFZhbHVlcyhzY3JvbGxZLm1pblJhbmdlLCBzY3JvbGxZLm1heFJhbmdlLCBzY3JvbGxZLm1pbiArIDEsIHNjcm9sbFkubWF4ICsgMSk7XHJcbiAgICAgIH1cclxuICAgICAgdGhpcy5tX25XaGVlbENvdW50ID0gMDtcclxuICAgIH1cclxuICAgIGVsc2UgaWYodGhpcy5tX25XaGVlbENvdW50ID09PSAtMSlcclxuICAgIHtcclxuICAgICAgLy9zY3JvbGwgLVxyXG4gICAgICBjb25zdCBzY3JvbGxZID0gZ3JpZC52ZXJ0U2Nyb2xsO1xyXG4gICAgICBpZihzY3JvbGxZLm1pbiA+PTEpIHtcclxuICAgICAgICBzY3JvbGxZLnNldFZhbHVlcyhzY3JvbGxZLm1pblJhbmdlLCBzY3JvbGxZLm1heFJhbmdlLCBzY3JvbGxZLm1pbiAtIDEsIHNjcm9sbFkubWF4IC0gMSk7XHJcbiAgICAgIH1cclxuICAgICAgdGhpcy5tX25XaGVlbENvdW50ID0gMDtcclxuICAgIH1cclxuICAgIGVsc2Uge1xyXG4gICAgICB0aGlzLm1fbldoZWVsQ291bnQgPSBlLmRlbHRhWSA+IDAgPyAxIDogLTE7XHJcbiAgICB9XHJcbiAgfVxyXG5cclxuXHJcbiAgcHJpdmF0ZSBwYWludChnIDogQ2FudmFzUmVuZGVyaW5nQ29udGV4dDJEIHwgbnVsbCwgZ3JpZCA6IERHLkdyaWQpIDogdm9pZCB7XHJcbiAgICAvL2NvbnN0IG5XRGl2ID0gZW50cnkuY29udGVudEJveFNpemUgPyBlbnRyeS5jb250ZW50Qm94U2l6ZVswXS5pbmxpbmVTaXplIDogZW50cnkuY29udGVudFJlY3Qud2lkdGg7XHJcblxyXG4gICAgaWYoZyA9PT0gbnVsbCkge1xyXG4gICAgICByZXR1cm47XHJcbiAgICB9XHJcblxyXG4gICAgaWYodGhpcy5tX3Jvb3QgPT09IG51bGwpIHtcclxuICAgICAgdGhyb3cgbmV3IEVycm9yKCdSb290IGNhbm5vdCBiZSBudWxsLicpO1xyXG4gICAgfVxyXG5cclxuICAgIGlmKHRoaXMubV9jb2xHcmlkID09PSBudWxsKSB7XHJcbiAgICAgIHRocm93IG5ldyBFcnJvcignQ29sdW1uIGdyaWQgY2Fubm90IGJlIG51bGwuJyk7XHJcbiAgICB9XHJcbiAgICBjb25zdCBkZnJhbWUgPSBncmlkLmRhdGFGcmFtZTtcclxuICAgIGNvbnN0IG5XID0gdGhpcy5tX3Jvb3Qub2Zmc2V0V2lkdGg7XHJcbiAgICBjb25zdCBuSCA9IHRoaXMubV9yb290Lm9mZnNldEhlaWdodDtcclxuXHJcbiAgICBnLmZpbGxTdHlsZSA9IFwid2hpdGVcIjtcclxuICAgIGcuZmlsbFJlY3QoMCwwLCBuVyp3aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbywgbkgqd2luZG93LmRldmljZVBpeGVsUmF0aW8pO1xyXG5cclxuICAgIGlmKHRoaXMubV9jb2xHcmlkLm5hbWUgPT09IG51bGwpXHJcbiAgICAgIHJldHVybjtcclxuXHJcbiAgICBjb25zdCBiaXRzZXRGaWx0ZXIgPSBkZnJhbWUuZmlsdGVyO1xyXG4gICAgaWYoYml0c2V0RmlsdGVyLmZhbHNlQ291bnQgPT09IGRmcmFtZS5yb3dDb3VudClcclxuICAgICAgcmV0dXJuOyAvL2V2ZXJ5dGhpbmcgaXMgZmlsdGVyZWRcclxuXHJcbiAgICAvL2NvbHVtbiBIZWFkZXJcclxuICAgIGNvbnN0IG9wdGlvbnMgOiBhbnkgPSBncmlkLmdldE9wdGlvbnModHJ1ZSk7XHJcblxyXG4gICAgY29uc3QgZm9udENlbGxEZWZhdWx0ID0gb3B0aW9ucy5sb29rLmRlZmF1bHRDZWxsRm9udDtcclxuXHJcbiAgICBsZXQgZm9udCA9IG9wdGlvbnMubG9vay5jb2xIZWFkZXJGb250ID09IG51bGwgfHwgb3B0aW9ucy5sb29rLmNvbEhlYWRlckZvbnQgPT09IHVuZGVmaW5lZCA/IFwiYm9sZCAxNHB4IFZvbHRhIFRleHQsIEFyaWFsXCIgOiBvcHRpb25zLmxvb2suY29sSGVhZGVyRm9udDtcclxuICAgIGxldCBmb250U2NhbGVkID0gR3JpZFV0aWxzLnNjYWxlRm9udChmb250LCB3aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbyk7XHJcbiAgICBnLmZvbnQgPSBmb250U2NhbGVkO1xyXG5cclxuICAgIGxldCBzdHIgPSBUZXh0VXRpbHMudHJpbVRleHQodGhpcy5tX2NvbEdyaWQubmFtZSwgZywgblcpO1xyXG5cclxuICAgIGNvbnN0IHRtID0gZy5tZWFzdXJlVGV4dChzdHIpO1xyXG4gICAgY29uc3QgbldMYWJlbCA9IHRtLndpZHRoO1xyXG5cclxuICAgIGNvbnN0IG5Bc2NlbnQgPSBNYXRoLmFicyh0bS5hY3R1YWxCb3VuZGluZ0JveEFzY2VudCk7XHJcbiAgICBjb25zdCBuRGVzY2VudCA9IHRtLmFjdHVhbEJvdW5kaW5nQm94RGVzY2VudDtcclxuICAgIGNvbnN0IG5IRm9udCA9ICBuQXNjZW50ICsgbkRlc2NlbnQ7Ly8gKyAyKm5ZSW5zZXQ7XHJcblxyXG4gICAgLy9sZXQgY2VsbENIID0gZ3JpZC5jZWxsKHRoaXMubV9jb2xHcmlkLm5hbWUsIC0xKTtcclxuICAgIC8vbGV0IHJlbmRlcmVyID0gY2VsbENILnJlbmRlcmVyO1xyXG5cclxuICAgIGxldCBuWCA9IDA7XHJcbiAgICBsZXQgblkgPSAwO1xyXG4gICAgY29uc3QgbkhDSCA9IEdyaWRVdGlscy5nZXRHcmlkQ29sdW1uSGVhZGVySGVpZ2h0KGdyaWQpKndpbmRvdy5kZXZpY2VQaXhlbFJhdGlvO1xyXG4gICAgZy50ZXh0QWxpZ24gPSAnc3RhcnQnO1xyXG4gICAgZy5maWxsU3R5bGUgPSBcIkJsYWNrXCI7XHJcbiAgICBsZXQgbllPZmZzZXQgPSBNYXRoLmZsb29yKChuSENIIC0gbkhGb250KS8yKTtcclxuICAgIGNvbnN0IG5YWCA9IG5YICsgKChuVyp3aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbyAtIG5XTGFiZWwpID4+IDEpO1xyXG4gICAgbGV0IG5ZWSA9IChuWSArIG5IQ0ggLSBNYXRoLmNlaWwoMyp3aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbykpOy8vLTIqd2luZG93LmRldmljZVBpeGVsUmF0aW8pO1xyXG4gICAgLy9vbnNvbGUubG9nKFwiblhYIFwiICsgblhYICsgXCIgbllZID0gXCIgKyBuWVkgKyBcIiBDSEggXCIgKyBuSENIKTtcclxuICAgIGcuZmlsbFRleHQoc3RyLCBuWFgsIG5ZWSk7XHJcblxyXG5cclxuICAgIC8vUmVndWxhciBjZWxsc1xyXG4gICAgY29uc3QgblJvd0N1cnJlbnQgPSAgZGZyYW1lLmN1cnJlbnRSb3cuaWR4O1xyXG4gICAgY29uc3QgYml0c2V0U2VsID0gZGZyYW1lLnNlbGVjdGlvbjtcclxuXHJcbiAgICBjb25zdCBhclJvd3NNaW5NYXggPSBbLTEsLTFdO1xyXG4gICAgR3JpZFV0aWxzLmZpbGxWaXNpYmxlVmlld3BvcnRSb3dzKGFyUm93c01pbk1heCwgZ3JpZCk7XHJcbiAgICBjb25zdCBuUm93TWluID0gYXJSb3dzTWluTWF4WzBdO1xyXG4gICAgY29uc3QgblJvd01heCA9IGFyUm93c01pbk1heFsxXTtcclxuXHJcbiAgICAvL2NvbnNvbGUubG9nKG5Sb3dNaW4gKyBcIiBcIiArIG5Sb3dNYXgpO1xyXG4gICAgY29uc3QgbkhSb3cgPSBHcmlkVXRpbHMuZ2V0R3JpZFJvd0hlaWdodChncmlkKTtcclxuICAgIG5ZT2Zmc2V0ID0gbkhDSDtcclxuICAgIGNvbnN0IG5IUm93R3JpZCA9IG5IUm93KndpbmRvdy5kZXZpY2VQaXhlbFJhdGlvO1xyXG4gICAgbGV0IGNlbGxSSCA9IG51bGw7XHJcblxyXG4gICAgbGV0IG5XVyA9IG5XKndpbmRvdy5kZXZpY2VQaXhlbFJhdGlvO1xyXG4gICAgLy9jb25zdCBuSEggPSBuSFJvd0dyaWQ7XHJcbiAgICAvL215IGNoYW5nZXMgY29uc3QgblBpbm5lZFJvd0NvdW50ID0gQXJyYXkuZnJvbShncmlkLnBpbm5lZFJvd3MpLmxlbmd0aDtcclxuXHJcbiAgICBjb25zdCBhclRhYmxlUm93cyA9IG5ldyBBcnJheShuUm93TWF4IC0gblJvd01pbiArMSk7XHJcbiAgICBsZXQgblJvd1RhYmxlID0gLTE7XHJcbiAgICBsZXQgYlNlbCA9IGZhbHNlO1xyXG4gICAgZm9yKGxldCBuUkc9blJvd01pbjsgblJHPD1uUm93TWF4OyArK25SRykge1xyXG4gICAgICB0cnkge1xyXG4gICAgICAgIGNlbGxSSCA9IGdyaWQuY2VsbCh0aGlzLm1fY29sR3JpZC5uYW1lLCBuUkcpO1xyXG4gICAgICB9IGNhdGNoIChlKSAvL3RvIGFkZHJlc3MgREcgYnVnIHdoZW4gZXZlcnl0aGluZyBpcyBmaWx0ZXJlZFxyXG4gICAgICB7XHJcbiAgICAgICAgY29udGludWU7XHJcbiAgICAgIH1cclxuXHJcbiAgICAgIGlmIChjZWxsUkgudGFibGVSb3dJbmRleCA9PT0gdW5kZWZpbmVkKS8vREcgYnVnXHJcbiAgICAgICAgY29udGludWU7XHJcblxyXG4gICAgICAvL215IGNoYW5nZXMgaWYodGhpcy5tX2NvbEdyaWQubmFtZSA9PSAnJylcclxuICAgICAgIC8vbXkgY2hhbmdlc2NlbGxSSC5jdXN0b21UZXh0ID0gblJHIC0gblJvd01pbiA8IG5QaW5uZWRSb3dDb3VudCA/ICcnIDogKG5SRyAtIG5QaW5uZWRSb3dDb3VudCArMSkudG9TdHJpbmcoKTtcclxuXHJcbiAgICAgIG5Sb3dUYWJsZSA9IGNlbGxSSC50YWJsZVJvd0luZGV4ID09PSBudWxsID8gLTEgOiBjZWxsUkgudGFibGVSb3dJbmRleDtcclxuICAgICAgYXJUYWJsZVJvd3NbblJHIC0gblJvd01pbl0gPSBuUm93VGFibGU7XHJcblxyXG4gICAgICBuWVkgPSBuWU9mZnNldCArIChuUkcgLSBuUm93TWluKSAqIG5IUm93R3JpZDtcclxuXHJcbiAgICAgIGxldCByZW5kZXJlcjogYW55ID0gR3JpZFV0aWxzLmdldEdyaWRDb2x1bW5SZW5kZXJlcihjZWxsUkguZ3JpZENvbHVtbik7XHJcbiAgICAgIGlmIChyZW5kZXJlciA9PT0gbnVsbCkge1xyXG4gICAgICAgIHRyeSB7XHJcbiAgICAgICAgICByZW5kZXJlciA9IGNlbGxSSC5yZW5kZXJlcjtcclxuICAgICAgICB9IGNhdGNoIChlKSB7XHJcbiAgICAgICAgICBjb25zb2xlLmVycm9yKFwiQ291bGQgbm90IG9idGFpbiByZW5kZXJlciBmb3IgREcgY2VsbC4gREcgYnVnIFwiICsgdGhpcy5tX2NvbEdyaWQubmFtZSArIFwiIHJvdyBcIiArIG5SRyk7XHJcbiAgICAgICAgICBjb250aW51ZTtcclxuICAgICAgICB9XHJcbiAgICAgIH1cclxuXHJcbiAgICAgIGlmIChyZW5kZXJlciA9PT0gbnVsbCB8fCByZW5kZXJlciA9PT0gdW5kZWZpbmVkKSB7XHJcbiAgICAgICAgY29uc29sZS5lcnJvcihcIkNvdWxkbid0IGZpbmQgcmVuZGVyZXIgZm9yIHBpbm5lZCBjb2x1bW4gXCIgKyB0aGlzLm1fY29sR3JpZC5uYW1lICsgXCIgcm93IFwiICsgblJHKTtcclxuICAgICAgICBjb250aW51ZTtcclxuICAgICAgfVxyXG5cclxuICAgICAgLy9sZXQgbllZID0gblk7Ly8qd2luZG93LmRldmljZVBpeGVsUmF0aW87XHJcblxyXG5cclxuICAgICAgZm9udCA9IGNlbGxSSC5zdHlsZS5mb250O1xyXG4gICAgICBmb250U2NhbGVkID0gR3JpZFV0aWxzLnNjYWxlRm9udChmb250LCB3aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbyk7XHJcbiAgICAgIGlmIChmb250U2NhbGVkICE9PSBudWxsKSB7XHJcbiAgICAgICAgY2VsbFJILnN0eWxlLmZvbnQgPSBmb250U2NhbGVkO1xyXG4gICAgICB9XHJcblxyXG4gICAgICBpZiAoblcgPiAwICYmIG5IUm93R3JpZCA+IDApIHsgLy90byBhZGRyZXNzIGEgYnVnIGNhdXNlZCBlaXRoZXIgREcgb3IgY2xpZW50IGFwcFxyXG4gICAgICAgIHRyeSB7XHJcbiAgICAgICAgICBpZiAocmVuZGVyZXIubmFtZSA9PT0gJ01vbGVjdWxlJykge1xyXG4gICAgICAgICAgICByZW5kZXJlci5yZW5kZXIoZywgMCwgbllZL3dpbmRvdy5kZXZpY2VQaXhlbFJhdGlvLCBuV1cvd2luZG93LmRldmljZVBpeGVsUmF0aW8sIG5IUm93R3JpZC93aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbywgY2VsbFJILCBjZWxsUkguc3R5bGUpO1xyXG4gICAgICAgICAgfVxyXG4gICAgICAgICAgZWxzZSByZW5kZXJlci5yZW5kZXIoZywgMCwgbllZLCBuV1csIG5IUm93R3JpZCwgY2VsbFJILCBjZWxsUkguc3R5bGUpO1xyXG5cclxuICAgICAgICB9IGNhdGNoIChlKSB7XHJcbiAgICAgICAgICBjb25zb2xlLmVycm9yKFwiQ291bGQgbm90IHBhaW50IGNlbGwgZm9yIHBpbm5lZCBjb2x1bW4gXCIgKyB0aGlzLm1fY29sR3JpZC5uYW1lICsgXCIgcm93IFwiICsgblJHKTtcclxuICAgICAgICAgIGNvbnRpbnVlO1xyXG4gICAgICAgICAgLy90aHJvdyBlO1xyXG4gICAgICAgIH1cclxuICAgICAgfVxyXG4gICAgfVxyXG5cclxuXHJcbiAgICAvL1BhaW50IEdyaWRcclxuICAgIGcuc3Ryb2tlU3R5bGUgPSBcIkdhaW5zYm9yb1wiO1xyXG4gICAgZy5iZWdpblBhdGgoKTtcclxuICAgIGcubW92ZVRvKDAsIG5ZKndpbmRvdy5kZXZpY2VQaXhlbFJhdGlvKTtcclxuICAgIGcubGluZVRvKDAsIChuWSArIG5IQ0gtMSp3aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbykpO1xyXG4gICAgZy5zdHJva2UoKTtcclxuXHJcbiAgICBnLmJlZ2luUGF0aCgpO1xyXG4gICAgZy5tb3ZlVG8oMCwgbllPZmZzZXQgKyAxKTtcclxuICAgIGcubGluZVRvKG5XVywgbllPZmZzZXQgKyAxKTtcclxuICAgIGcuc3Ryb2tlKCk7XHJcblxyXG4gICAgY29uc3QgblBpbm5lZENvbENvdW50ID0gUGlubmVkVXRpbHMuZ2V0UGlubmVkQ29sdW1uQ291bnQoZ3JpZCk7XHJcbiAgICBjb25zdCBjb2xQaW5uZWQgPSBQaW5uZWRVdGlscy5nZXRQaW5uZWRDb2x1bW4oblBpbm5lZENvbENvdW50IC0xLCBncmlkKTtcclxuICAgIGNvbnN0IGJMYXN0ID0gdGhpcyA9PT0gY29sUGlubmVkO1xyXG5cclxuICAgIGZvcihsZXQgblJHPW5Sb3dNaW47IG5SRzw9blJvd01heDsgKytuUkcpXHJcbiAgICB7XHJcbiAgICAgIG5ZWSA9IG5ZT2Zmc2V0ICsgKG5SRyAtIG5Sb3dNaW4pICogbkhSb3dHcmlkO1xyXG4gICAgICAvL2lmKG9wdGlvbnMubG9vay5zaG93Um93R3JpZGxpbmVzKSB7XHJcbiAgICAgICAgZy5zdHJva2VTdHlsZSA9IFwiR2FpbnNib3JvXCI7XHJcbiAgICAgICAgZy5iZWdpblBhdGgoKTtcclxuICAgICAgICBnLm1vdmVUbygwLCBuWVkgKyBuSFJvd0dyaWQrMSk7XHJcbiAgICAgICAgZy5saW5lVG8obldXLCBuWVkgKyBuSFJvd0dyaWQrMSk7XHJcbiAgICAgICAgZy5zdHJva2UoKTtcclxuXHJcbiAgICAgICAgZy5iZWdpblBhdGgoKTtcclxuICAgICAgICBnLm1vdmVUbygwLCBuWVkpO1xyXG4gICAgICAgIGcubGluZVRvKDAsIG5ZWSArIG5IUm93R3JpZCsxKTtcclxuICAgICAgICBnLnN0cm9rZSgpO1xyXG5cclxuICAgICAgICBpZihiTGFzdCl7Ly9teSBjaGFuZ2VzICAmJiAoblJHIC0gblJvd01pbikgPj0gblBpbm5lZFJvd0NvdW50KSB7XHJcbiAgICAgICAgICBnLnN0cm9rZVN0eWxlID0gXCJibGFja1wiO1xyXG4gICAgICAgICAgZy5iZWdpblBhdGgoKTtcclxuICAgICAgICAgIGcubW92ZVRvKG5XVyAtIDEsIG5ZWSk7XHJcbiAgICAgICAgICBnLmxpbmVUbyhuV1cgLSAxLCBuWVkgKyBuSFJvd0dyaWQgKyAxKTtcclxuICAgICAgICAgIGcuc3Ryb2tlKCk7XHJcbiAgICAgICAgfVxyXG5cclxuICAgICAgLy99XHJcbiAgICAgIG5Sb3dUYWJsZSA9IGFyVGFibGVSb3dzW25SRyAtIG5Sb3dNaW5dO1xyXG4gICAgICB0cnl7YlNlbCA9IG5Sb3dUYWJsZSA9PT0gdW5kZWZpbmVkIHx8IG5Sb3dUYWJsZSA8IDAgPyBmYWxzZSA6IGJpdHNldFNlbC5nZXQoblJvd1RhYmxlKTt9XHJcbiAgICAgIGNhdGNoIChlKXtcclxuICAgICAgICBjb25zb2xlLmVycm9yKCdQYWludEVycm9yOiByb3dfbWluOiAnICsgblJvd01pbiArICcgcm93X21heDogJyArIG5Sb3dNYXggKyAnIG5SICcgKyBuUkcgKyAnICcgKyBuUm93VGFibGUpO1xyXG4gICAgICAgIHRocm93IGU7XHJcbiAgICAgIH1cclxuICAgICAgaWYoYlNlbClcclxuICAgICAge1xyXG4gICAgICAgIGcuZ2xvYmFsQWxwaGEgPSAwLjI7XHJcbiAgICAgICAgZy5maWxsU3R5bGUgPSBQaW5uZWRDb2x1bW4uU0VMRUNUSU9OX0NPTE9SO1xyXG4gICAgICAgIGcuZmlsbFJlY3QoMCwgbllZLCBuV1csIG5IUm93R3JpZCk7XHJcbiAgICAgICAgZy5nbG9iYWxBbHBoYSA9IDE7XHJcbiAgICAgIH1cclxuXHJcbiAgICAgIGlmKG5Sb3dDdXJyZW50ID09PSBuUm93VGFibGUpXHJcbiAgICAgIHtcclxuICAgICAgICBnLmdsb2JhbEFscGhhID0gMC4yO1xyXG4gICAgICAgIGcuZmlsbFN0eWxlID0gUGlubmVkQ29sdW1uLkFDVElWRV9DRUxMX0NPTE9SO1xyXG4gICAgICAgIGcuZmlsbFJlY3QoMCwgbllZLCBuV1csIG5IUm93R3JpZCk7XHJcbiAgICAgICAgZy5nbG9iYWxBbHBoYSA9IDE7XHJcbiAgICAgIH1cclxuICAgIH0vL2ZvclxyXG4gIH1cclxuXHJcblxyXG4gIHByaXZhdGUgc3RhdGljIGhpdFRlc3RSb3dzKGVDYW52YXNQaW5uZWQgOiBIVE1MQ2FudmFzRWxlbWVudCwgZ3JpZCA6IERHLkdyaWQsIGUgOiBNb3VzZUV2ZW50LCBiQm9yZGVyIDogYm9vbGVhbiwgYXJYWU9uQ2VsbCA6IEFycmF5PG51bWJlcj4gfCB1bmRlZmluZWQpXHJcbiAge1xyXG4gICAgY29uc3QgcmVjdCA9IGVDYW52YXNQaW5uZWQuZ2V0Qm91bmRpbmdDbGllbnRSZWN0KCk7XHJcbiAgICBjb25zdCBzY3JvbGxMZWZ0PSB3aW5kb3cucGFnZVhPZmZzZXQgfHwgZG9jdW1lbnQuZG9jdW1lbnRFbGVtZW50LnNjcm9sbExlZnQ7XHJcbiAgICBjb25zdCBzY3JvbGxUb3AgPSB3aW5kb3cucGFnZVlPZmZzZXQgfHwgZG9jdW1lbnQuZG9jdW1lbnRFbGVtZW50LnNjcm9sbFRvcDtcclxuICAgIGNvbnN0IG5ZID0gcmVjdC50b3AgICsgc2Nyb2xsVG9wO1xyXG4gICAgY29uc3QgblggPSByZWN0LmxlZnQgKyBzY3JvbGxMZWZ0O1xyXG5cclxuICAgIGlmKG5YIDw9IGUuY2xpZW50WCAmJiBlLmNsaWVudFggPD0gblggKyBlQ2FudmFzUGlubmVkLm9mZnNldFdpZHRoKSAgIC8vb24gdGhlIHJvd3MgaGVhZGVyXHJcbiAgICB7XHJcbiAgICAgIGNvbnN0IG5ISGVhZGVyQ29scyA9IEdyaWRVdGlscy5nZXRHcmlkQ29sdW1uSGVhZGVySGVpZ2h0KGdyaWQpO1xyXG4gICAgICBjb25zdCBuSFJvd0dyaWQgPSBHcmlkVXRpbHMuZ2V0R3JpZFJvd0hlaWdodChncmlkKTtcclxuXHJcbiAgICAgIGNvbnN0IGFyTWluTWF4Um93cyA9IFstMSwtMV07XHJcbiAgICAgIEdyaWRVdGlscy5maWxsVmlzaWJsZVZpZXdwb3J0Um93cyhhck1pbk1heFJvd3MsIGdyaWQpO1xyXG4gICAgICBjb25zdCBuUm93TWluID0gYXJNaW5NYXhSb3dzWzBdO1xyXG4gICAgICBjb25zdCBuUm93TWF4ID0gYXJNaW5NYXhSb3dzWzFdO1xyXG5cclxuICAgICAgY29uc3QgbllNb3VzZU9uSGVhZGVyID0gZS5jbGllbnRZIC0gblk7XHJcblxyXG4gICAgICBsZXQgbllCb3JkZXIgPSAtMTtcclxuICAgICAgbGV0IG5ZRGlmZiA9IC0xO1xyXG5cclxuICAgICAgZm9yKGxldCBuUm93PW5Sb3dNaW47IG5Sb3c8PSBuUm93TWF4OyArK25Sb3cpXHJcbiAgICAgIHtcclxuICAgICAgICBuWUJvcmRlciA9IG5ISGVhZGVyQ29scyArIChuUm93IC0gblJvd01pbisxKSpuSFJvd0dyaWQ7XHJcbiAgICAgICAgbllEaWZmID0gbllNb3VzZU9uSGVhZGVyIC0gbllCb3JkZXI7XHJcblxyXG4gICAgICAgIGlmKGJCb3JkZXIgJiYgTWF0aC5hYnMobllEaWZmKSA8PSBQaW5uZWRDb2x1bW4uWV9SRVNJWkVfU0VOU0lUSVZJVFkpXHJcbiAgICAgICAge1xyXG4gICAgICAgICAgcmV0dXJuIG5Sb3c7XHJcbiAgICAgICAgfVxyXG5cclxuICAgICAgICBpZighYkJvcmRlciAmJiBuWUJvcmRlciAtIG5IUm93R3JpZCA8PSBuWU1vdXNlT25IZWFkZXIgJiYgbllNb3VzZU9uSGVhZGVyIDw9IG5ZQm9yZGVyKSB7XHJcblxyXG4gICAgICAgICAgaWYoYXJYWU9uQ2VsbCAhPT0gdW5kZWZpbmVkKSB7XHJcbiAgICAgICAgICAgIGFyWFlPbkNlbGxbMF0gPSBlLmNsaWVudFggLSBuWDtcclxuICAgICAgICAgICAgYXJYWU9uQ2VsbFsxXSA9IG5ZTW91c2VPbkhlYWRlciAtIG5ZQm9yZGVyICsgbkhSb3dHcmlkO1xyXG4gICAgICAgICAgfVxyXG5cclxuICAgICAgICAgIHJldHVybiBuUm93O1xyXG4gICAgICAgIH1cclxuICAgICAgfVxyXG4gICAgfVxyXG5cclxuICAgIHJldHVybiAtMTtcclxuICB9XHJcbn1cclxuIl19