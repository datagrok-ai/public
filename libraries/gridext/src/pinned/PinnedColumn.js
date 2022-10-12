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
        this.m_handlerKeyDown = rxjs.fromEvent(eCanvasThis, 'mousemove').subscribe((e) => {
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
        /*
           if(e.button === 2) {
        
             const eDivPOpup : HTMLElement | null = GridUtils.getGridDartPopupMenu();
             eDivPOpup?.setAttribute('column_name', this.m_colGrid.name);
             let d = 0;
             return;
           }*/
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
            const bAddToSel = e.ctrlKey;
            const bRangeSel = e.shiftKey;
            const nRowGrid = PinnedColumn.hitTestRows(this.m_root, grid, e, false, this.m_arXYMouseOnCellUp);
            if (!bAddToSel && !bRangeSel && nRowGrid === this.m_nRowGridDragging) { //click on the same row which will become active
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
                if (!bAddToSel || bRangeSel)
                    bitsetSel.setAll(false, true);
                let nRowMin = this.m_nRowGridDragging < nRowGrid ? this.m_nRowGridDragging : nRowGrid;
                let nRowMax = this.m_nRowGridDragging > nRowGrid ? this.m_nRowGridDragging : nRowGrid;
                if (bRangeSel) {
                    let nRowGridActive = GridUtils.getActiveGridRow(grid);
                    if (nRowGridActive === null)
                        nRowGridActive = 0;
                    nRowMin = nRowGridActive < nRowGrid ? nRowGridActive : nRowGrid;
                    nRowMax = nRowGridActive > nRowGrid ? nRowGridActive : nRowGrid;
                }
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
                        bitsetSel.set(nRowTable, true, true);
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
        //if(options.look.showRowGridlines) {
        //}
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
        for (let nRG = nRowMin; nRG <= nRowMax; ++nRG) {
            nYY = nYOffset + (nRG - nRowMin) * nHRowGrid;
            //if(options.look.showRowGridlines) {
            g.beginPath();
            g.moveTo(0, nYY + nHRowGrid + 1);
            g.lineTo(nWW, nYY + nHRowGrid + 1);
            g.stroke();
            g.beginPath();
            g.moveTo(0, nYY);
            g.lineTo(0, nYY + nHRowGrid + 1);
            g.stroke();
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
        }
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
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJmaWxlIjoiUGlubmVkQ29sdW1uLmpzIiwic291cmNlUm9vdCI6IiIsInNvdXJjZXMiOlsiUGlubmVkQ29sdW1uLnRzIl0sIm5hbWVzIjpbXSwibWFwcGluZ3MiOiJBQUFBLE9BQU8sS0FBSyxJQUFJLE1BQU0sbUJBQW1CLENBQUM7QUFDMUMsT0FBTyxLQUFLLEVBQUUsTUFBTSxpQkFBaUIsQ0FBQztBQUN0QyxPQUFPLEtBQUssRUFBRSxNQUFNLGlCQUFpQixDQUFDO0FBQ3RDLE9BQU8sS0FBSyxTQUFTLE1BQU0sb0JBQW9CLENBQUM7QUFDaEQsT0FBTyxLQUFLLFNBQVMsTUFBTSxvQkFBb0IsQ0FBQztBQUNoRCxPQUFPLEVBQUMsVUFBVSxFQUFDLE1BQU0scUJBQXFCLENBQUM7QUFDL0MsT0FBTyxLQUFLLElBQUksTUFBTSxNQUFNLENBQUM7QUFDN0IsT0FBTyxFQUFFLGtCQUFrQixFQUFDLE1BQU0sZ0NBQWdDLENBQUM7QUFDbkUsT0FBTyxLQUFLLFdBQVcsTUFBTSxlQUFlLENBQUM7QUFFN0MsT0FBTyxFQUFDLGVBQWUsRUFBQyxNQUFNLHVCQUF1QixDQUFDO0FBQ3RELE9BQU8sRUFBYyxNQUFNLEVBQUMsTUFBTSxpQkFBaUIsQ0FBQztBQUNwRCw0Q0FBNEM7QUFHNUM7Ozs7Ozs7Ozs7Ozs7Ozs7RUFnQkU7QUFFRixTQUFTLFdBQVcsQ0FBQyxJQUFrQjtJQUNyQyxNQUFNLE9BQU8sR0FBRyxJQUFJLENBQUMsVUFBVSxDQUFDO0lBQ2hDLElBQUksT0FBTyxLQUFLLElBQUksSUFBSSxPQUFPLEtBQUssU0FBUyxFQUFFO1FBQzdDLE1BQU0sSUFBSSxLQUFLLENBQUMsNENBQTRDLENBQUMsQ0FBQztLQUMvRDtJQUVELElBQUksUUFBUSxHQUFHLFNBQVMsQ0FBQyxxQkFBcUIsQ0FBQyxPQUFPLENBQUMsQ0FBQztJQUN4RCxJQUFHLFFBQVEsWUFBWSxrQkFBa0IsRUFBRTtRQUN6QyxPQUFPLFFBQVEsQ0FBQztLQUNqQjtJQUVELE9BQU8sSUFBSSxDQUFDLFFBQVEsQ0FBQztBQUN2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QU0sR0FBRyxlQUFlLEdBQUcsUUFBUSxDQUFDO2dCQUVwQyxJQUFHLE9BQU8sSUFBSSxJQUFJLENBQUMsR0FBRyxDQUFDLE1BQU0sQ0FBQyxJQUFJLFlBQVksQ0FBQyxvQkFBb0IsRUFDbkU7b0JBQ0UsT0FBTyxJQUFJLENBQUM7aUJBQ2I7Z0JBRUQsSUFBRyxDQUFDLE9BQU8sSUFBSSxRQUFRLEdBQUcsU0FBUyxJQUFJLGVBQWUsSUFBSSxlQUFlLElBQUksUUFBUSxFQUFFO29CQUVyRixJQUFHLFVBQVUsS0FBSyxTQUFTLEVBQUU7d0JBQzNCLFVBQVUsQ0FBQyxDQUFDLENBQUMsR0FBRyxDQUFDLENBQUMsT0FBTyxHQUFHLEVBQUUsQ0FBQzt3QkFDL0IsVUFBVSxDQUFDLENBQUMsQ0FBQyxHQUFHLGVBQWUsR0FBRyxRQUFRLEdBQUcsU0FBUyxDQUFDO3FCQUN4RDtvQkFFRCxPQUFPLElBQUksQ0FBQztpQkFDYjthQUNGO1NBQ0Y7UUFFRCxPQUFPLENBQUMsQ0FBQyxDQUFDO0lBQ1osQ0FBQzs7QUEzbENjLDJCQUFjLEdBQUcsRUFBRSxDQUFDO0FBQ3BCLDJCQUFjLEdBQUcsR0FBRyxDQUFDO0FBQ3JCLDRCQUFlLEdBQUcsVUFBVSxDQUFDLEtBQUssQ0FBQyxVQUFVLENBQUMsWUFBWSxDQUFDLENBQUMsQ0FBQyw2QkFBNkI7QUFDMUYsOEJBQWlCLEdBQUcsVUFBVSxDQUFDLEtBQUssQ0FBQyxVQUFVLENBQUMsVUFBVSxDQUFDLENBQUMsQ0FBQyw2QkFBNkI7QUFDMUYsaUNBQW9CLEdBQUcsQ0FBQyxDQUFDIiwic291cmNlc0NvbnRlbnQiOlsiaW1wb3J0ICogYXMgZ3JvayBmcm9tICdkYXRhZ3Jvay1hcGkvZ3Jvayc7XHJcbmltcG9ydCAqIGFzIERHIGZyb20gJ2RhdGFncm9rLWFwaS9kZyc7XHJcbmltcG9ydCAqIGFzIHVpIGZyb20gJ2RhdGFncm9rLWFwaS91aSc7XHJcbmltcG9ydCAqIGFzIEdyaWRVdGlscyBmcm9tICcuLi91dGlscy9HcmlkVXRpbHMnO1xyXG5pbXBvcnQgKiBhcyBUZXh0VXRpbHMgZnJvbSAnLi4vdXRpbHMvVGV4dFV0aWxzJztcclxuaW1wb3J0IHtDb2xvclV0aWxzfSBmcm9tICcuLi91dGlscy9Db2xvclV0aWxzJztcclxuaW1wb3J0ICogYXMgcnhqcyBmcm9tICdyeGpzJztcclxuaW1wb3J0IHsgR3JpZENlbGxSZW5kZXJlckV4fSBmcm9tIFwiLi4vcmVuZGVyZXIvR3JpZENlbGxSZW5kZXJlckV4XCI7XHJcbmltcG9ydCAqIGFzIFBpbm5lZFV0aWxzIGZyb20gXCIuL1Bpbm5lZFV0aWxzXCI7XHJcbmltcG9ydCB7Z2V0R3JpZERhcnRQb3B1cE1lbnUsIGlzSGl0VGVzdE9uRWxlbWVudH0gZnJvbSBcIi4uL3V0aWxzL0dyaWRVdGlsc1wiO1xyXG5pbXBvcnQge01vdXNlRGlzcGF0Y2hlcn0gZnJvbSBcIi4uL3VpL01vdXNlRGlzcGF0Y2hlclwiO1xyXG5pbXBvcnQge0NvbHVtbnNBcmdzLCB0b0RhcnR9IGZyb20gXCJkYXRhZ3Jvay1hcGkvZGdcIjtcclxuLy9pbXBvcnQge1RhYmxlVmlld30gZnJvbSBcImRhdGFncm9rLWFwaS9kZ1wiO1xyXG5cclxuXHJcbi8qXHJcbmNvbnN0IGhTdWJzY3JpYmVyICA9IGdyb2suZXZlbnRzLm9uVmlld0xheW91dEFwcGxpZWQuc3Vic2NyaWJlKChsYXlvdXQgOiBERy5WaWV3TGF5b3V0KSA9PiB7XHJcbiAgY29uc3QgdmlldyA6IERHLlRhYmxlVmlldyA9IGxheW91dC52aWV3IGFzIFRhYmxlVmlldztcclxuICBjb25zdCBpdFZpZXdlcnMgPSB2aWV3LnZpZXdlcnM7XHJcbiAgY29uc3QgYXJWaWV3ZXJzID0gQXJyYXkuZnJvbShpdFZpZXdlcnMpO1xyXG5cclxuICBsZXQgdmlld2VyID0gbnVsbDtcclxuICBjb25zdCBuVmlld2VyQ291bnQgPSBhclZpZXdlcnMubGVuZ3RoO1xyXG4gIGZvciAobGV0IG4gPSAwOyBuIDwgblZpZXdlckNvdW50OyArK24pIHtcclxuICAgIHZpZXdlciA9IGFyVmlld2Vyc1tuXTtcclxuICAgIGlmICh2aWV3ZXIudHlwZSAhPT0gXCJHcmlkXCIpXHJcbiAgICAgIGNvbnRpbnVlO1xyXG5cclxuICAgIFBpbm5lZFV0aWxzLmluc3RhbGxQaW5uZWRDb2x1bW5zKHZpZXdlciBhcyBERy5HcmlkKTtcclxuICB9XHJcbn0pO1xyXG4qL1xyXG5cclxuZnVuY3Rpb24gZ2V0UmVuZGVyZXIoY2VsbCA6IERHLkdyaWRDZWxsKSA6IEdyaWRDZWxsUmVuZGVyZXJFeCB8IERHLkdyaWRDZWxsUmVuZGVyZXIge1xyXG4gIGNvbnN0IGNvbEdyaWQgPSBjZWxsLmdyaWRDb2x1bW47XHJcbiAgaWYgKGNvbEdyaWQgPT09IG51bGwgfHwgY29sR3JpZCA9PT0gdW5kZWZpbmVkKSB7XHJcbiAgICB0aHJvdyBuZXcgRXJyb3IoJ0dyaWQgY2VsbCBpcyBkZXRhY2hlZCBmcm9tIHRoZSBHcmlkIGNvbHVtbicpO1xyXG4gIH1cclxuXHJcbiAgbGV0IHJlbmRlcmVyID0gR3JpZFV0aWxzLmdldEdyaWRDb2x1bW5SZW5kZXJlcihjb2xHcmlkKTtcclxuICBpZihyZW5kZXJlciBpbnN0YW5jZW9mIEdyaWRDZWxsUmVuZGVyZXJFeCkge1xyXG4gICAgcmV0dXJuIHJlbmRlcmVyO1xyXG4gIH1cclxuXHJcbiAgcmV0dXJuIGNlbGwucmVuZGVyZXI7XHJcbn1cclxuXHJcblxyXG5mdW5jdGlvbiBnZXRHcmlkKGNvbEdyaWQgOiBERy5HcmlkQ29sdW1uKSA6IERHLkdyaWQgfCBudWxsIHtcclxuICBsZXQgZ3JpZCA6IERHLkdyaWQgfCBudWxsID0gY29sR3JpZC5ncmlkO1xyXG4gIGlmKCBncmlkID09PSBudWxsKSB7XHJcbiAgICBncmlkID0gR3JpZFV0aWxzLmdldEluc3RhbGxlZEdyaWRGb3JDb2x1bW4oY29sR3JpZCk7XHJcbiAgICBpZihncmlkIGluc3RhbmNlb2YgREcuR3JpZClcclxuICAgICAgcmV0dXJuIGdyaWQ7XHJcbiAgfVxyXG5cclxuICByZXR1cm4gZ3JpZDtcclxufVxyXG5cclxuXHJcbmZ1bmN0aW9uIG5vdGlmeUFsbENvbHNSb3dzUmVzaXplZChncmlkIDogREcuR3JpZCwgbkhSb3dzIDogbnVtYmVyLCBiQWRqdXN0aW5nIDogYm9vbGVhbikgOiB2b2lkIHtcclxuXHJcbiAgbGV0IHJlbmRlcmVyIDogR3JpZENlbGxSZW5kZXJlckV4IHwgbnVsbCA9IG51bGxcclxuICBsZXQgY29sR3JpZCA9IG51bGw7XHJcbiAgY29uc3QgbHN0Q29sc0dyaWQgPSBncmlkLmNvbHVtbnM7XHJcbiAgY29uc3QgbkNvbENvdW50ID0gbHN0Q29sc0dyaWQubGVuZ3RoO1xyXG4gIGZvcihsZXQgbkNvbD0wOyBuQ29sPG5Db2xDb3VudDsgKytuQ29sKSB7XHJcbiAgICBjb2xHcmlkID0gbHN0Q29sc0dyaWQuYnlJbmRleChuQ29sKTtcclxuICAgIGlmKGNvbEdyaWQgPT09IG51bGwgfHwgIWNvbEdyaWQudmlzaWJsZSl7XHJcbiAgICAgIGNvbnRpbnVlXHJcbiAgICB9XHJcblxyXG4gICAgcmVuZGVyZXIgPSBHcmlkVXRpbHMuZ2V0R3JpZENvbHVtblJlbmRlcmVyKGNvbEdyaWQpO1xyXG4gICAgaWYgKHJlbmRlcmVyIGluc3RhbmNlb2YgR3JpZENlbGxSZW5kZXJlckV4KSB7XHJcbiAgICAgIHJlbmRlcmVyLm9uUmVzaXplSGVpZ2h0KGNvbEdyaWQsIGdyaWQsIG5IUm93cywgYkFkanVzdGluZyk7XHJcbiAgICB9XHJcbiAgfVxyXG59XHJcblxyXG5cclxuZnVuY3Rpb24gbm90aWZ5QWxsUGlubmVkQ29sc1Jvd3NSZXNpemVkKGNvbFBpbm5lZFNvdXJjZSA6IFBpbm5lZENvbHVtbiwgbkhSb3dzIDogbnVtYmVyLCBiQWRqdXN0aW5nIDogYm9vbGVhbikgOiB2b2lkIHtcclxuXHJcbiAgY29uc3QgY29sR3JpZFNvdXJjZSAgPSBjb2xQaW5uZWRTb3VyY2UuZ2V0R3JpZENvbHVtbigpO1xyXG4gIGlmKGNvbEdyaWRTb3VyY2UgPT09IG51bGwpe1xyXG4gICAgcmV0dXJuO1xyXG4gIH1cclxuXHJcbiAgY29uc3QgZ3JpZCA9IGdldEdyaWQoY29sR3JpZFNvdXJjZSk7XHJcbiAgY29uc3QgZGFydCA9IERHLnRvRGFydChncmlkKTtcclxuICBpZihkYXJ0Lm1fYXJQaW5uZWRDb2xzID09PSB1bmRlZmluZWQpIHtcclxuICAgIHRocm93IG5ldyBFcnJvcignUGlubmVkIENvbHVtbnMgYXJlIG5vdCBpbnN0YWxsZWQuJyk7XHJcbiAgfVxyXG5cclxuICBsZXQgcmVuZGVyZXIgOiBHcmlkQ2VsbFJlbmRlcmVyRXggfCBudWxsID0gbnVsbFxyXG4gIGxldCBjb2xQaW5uZWQgPSBudWxsO1xyXG4gIGxldCBjb2xHcmlkID0gbnVsbDtcclxuICBjb25zdCBuUGlubmVkQ29sQ291bnQgPSBkYXJ0Lm1fYXJQaW5uZWRDb2xzLmxlbmd0aDtcclxuICBmb3IobGV0IG5Db2xQaW49MDsgbkNvbFBpbjxuUGlubmVkQ29sQ291bnQ7ICsrbkNvbFBpbikge1xyXG4gICAgY29sUGlubmVkID0gZGFydC5tX2FyUGlubmVkQ29sc1tuQ29sUGluXTtcclxuICAgIGNvbEdyaWQgPSBjb2xQaW5uZWQubV9jb2xHcmlkO1xyXG4gICAgaWYoY29sR3JpZCA9PT0gbnVsbCkge1xyXG4gICAgICB0aHJvdyBuZXcgRXJyb3IoJ1Bpbm5lZCBDb2x1bW4gaXMgZGV0YWNoZWQuJyk7XHJcbiAgICB9XHJcblxyXG4gICAgcmVuZGVyZXIgPSBHcmlkVXRpbHMuZ2V0R3JpZENvbHVtblJlbmRlcmVyKGNvbEdyaWQpO1xyXG4gICAgaWYgKHJlbmRlcmVyIGluc3RhbmNlb2YgR3JpZENlbGxSZW5kZXJlckV4ICAmJiBjb2xQaW5uZWQubV9yb290ICE9PSBudWxsICYmIGdyaWQgIT09IG51bGwpIHtcclxuICAgICAgcmVuZGVyZXIub25SZXNpemVIZWlnaHQoY29sUGlubmVkLCBncmlkLCBuSFJvd3MsIGJBZGp1c3RpbmcpO1xyXG4gICAgfVxyXG4gIH1cclxufVxyXG5cclxuXHJcbmNvbnN0IERFQlVHIDogYm9vbGVhbiA9IGZhbHNlO1xyXG5cclxuXHJcbmV4cG9ydCBjbGFzcyBQaW5uZWRDb2x1bW4ge1xyXG5cclxuICBwcml2YXRlIHN0YXRpYyBNSU5fUk9XX0hFSUdIVCA9IDIwO1xyXG4gIHByaXZhdGUgc3RhdGljIE1BWF9ST1dfSEVJR0hUID0gNTAwO1xyXG4gIHByaXZhdGUgc3RhdGljIFNFTEVDVElPTl9DT0xPUiA9IENvbG9yVXRpbHMudG9SZ2IoQ29sb3JVdGlscy5jb2xTZWxlY3Rpb24pOyAvL1wicmdiYSgyMzcsIDIyMCwgODgsIDAuMTUpXCI7XHJcbiAgcHJpdmF0ZSBzdGF0aWMgQUNUSVZFX0NFTExfQ09MT1IgPSBDb2xvclV0aWxzLnRvUmdiKENvbG9yVXRpbHMuY3VycmVudFJvdyk7IC8vXCJyZ2JhKDE1MywgMjM3LCA4MiwgMC4yNSlcIjtcclxuICBwcml2YXRlIHN0YXRpYyBZX1JFU0laRV9TRU5TSVRJVklUWSA9IDI7XHJcblxyXG4gIHByaXZhdGUgbV9mRGV2aWNlUGl4ZWxSYXRpbyA6IG51bWJlcjtcclxuICBwcml2YXRlIG1fY29sR3JpZCA6IERHLkdyaWRDb2x1bW4gfCBudWxsO1xyXG4gIHByaXZhdGUgbV9yb290IDogSFRNTENhbnZhc0VsZW1lbnQgfCBudWxsO1xyXG4gIHByaXZhdGUgbV9uV2lkdGhCdWcgOiBudW1iZXI7XHJcbiAgLy9wcml2YXRlIG1fb2JzZXJ2ZXJSZXNpemUgOiBSZXNpemVPYnNlcnZlciB8IG51bGw7XHJcbiAgcHJpdmF0ZSBtX29ic2VydmVyUmVzaXplR3JpZCA6IFJlc2l6ZU9ic2VydmVyIHwgbnVsbDtcclxuICBwcml2YXRlIG1faGFuZGxlcktleURvd24gOiByeGpzLlN1YnNjcmlwdGlvbiB8IG51bGw7XHJcbiAgcHJpdmF0ZSBtX2hhbmRsZXJDb2xzUmVtb3ZlZCA6IHJ4anMuU3Vic2NyaXB0aW9uIHwgbnVsbDtcclxuICBwcml2YXRlIG1faGFuZGxlckNvbE5hbWVDaGFuZ2VkIDogcnhqcy5TdWJzY3JpcHRpb24gfCBudWxsO1xyXG4gIHByaXZhdGUgbV9oYW5kbGVyVlNjcm9sbCA6IHJ4anMuU3Vic2NyaXB0aW9uIHwgbnVsbDtcclxuICBwcml2YXRlIG1faGFuZGxlclJvd3NGaWx0ZXJpbmcgOiByeGpzLlN1YnNjcmlwdGlvbiB8IG51bGw7XHJcbiAgcHJpdmF0ZSBtX2hhbmRsZXJDdXJyUm93IDogcnhqcy5TdWJzY3JpcHRpb24gfCBudWxsO1xyXG4gIHByaXZhdGUgbV9oYW5kbGVyU2VsIDogcnhqcy5TdWJzY3JpcHRpb24gfCBudWxsO1xyXG4gIC8vcHJpdmF0ZSBtX2hhbmRsZXJGaWx0ZXIgOiBhbnk7XHJcbiAgcHJpdmF0ZSBtX2hhbmRsZXJSb3dzUmVzaXplZCA6IHJ4anMuU3Vic2NyaXB0aW9uIHwgbnVsbDtcclxuICBwcml2YXRlIG1faGFuZGxlclJvd3NTb3J0ZWQgOiByeGpzLlN1YnNjcmlwdGlvbiB8IG51bGw7XHJcblxyXG4gIHByaXZhdGUgbV9uSFJlc2l6ZVJvd3NCZWZvcmVEcmFnID0gLTE7XHJcbiAgcHJpdmF0ZSBtX25SZXNpemVSb3dHcmlkRHJhZ2dpbmcgPSAtMTtcclxuICBwcml2YXRlIG1fbllSZXNpemVEcmFnZ2luZ0FuY2hvciA9IC0xO1xyXG4gIHByaXZhdGUgbV9uUmVzaXplUm93R3JpZE1vdmluZyA9IC0xO1xyXG5cclxuICBwcml2YXRlIG1fbllEcmFnZ2luZ0FuY2hvciA9IC0xO1xyXG4gIHByaXZhdGUgbV9uUm93R3JpZERyYWdnaW5nID0gLTE7XHJcblxyXG4gIHByaXZhdGUgbV9uV2hlZWxDb3VudCA6IG51bWJlciA9IDA7XHJcblxyXG5cclxuICBwcml2YXRlIG1fYXJYWU1vdXNlT25DZWxsRG93biA9IFstMiwgLTJdO1xyXG4gIHByaXZhdGUgbV9hclhZTW91c2VPbkNlbGxVcCA9IFstMSwgLTFdO1xyXG4gIHByaXZhdGUgbV9iU29ydGVkQXNjZW5kaW5nIDogYm9vbGVhbiB8IG51bGwgPSBudWxsO1xyXG5cclxuICBwcml2YXRlIG1fY2VsbEN1cnJlbnQgOiBERy5HcmlkQ2VsbCB8IG51bGwgPSBudWxsO1xyXG5cclxuICBjb25zdHJ1Y3Rvcihjb2xHcmlkIDogREcuR3JpZENvbHVtbikge1xyXG5cclxuICAgIE1vdXNlRGlzcGF0Y2hlci5jcmVhdGUoKTtcclxuXHJcbiAgICBjb25zdCBncmlkID0gZ2V0R3JpZChjb2xHcmlkKTtcclxuICAgIGlmKGdyaWQgPT09IG51bGwpIHtcclxuICAgICAgdGhyb3cgbmV3IEVycm9yKFwiQ29sdW1uICdcIiArIGNvbEdyaWQubmFtZSArIFwiJyBpcyBub3QgYXR0YWNoZWQgdG8gdGhlIGdyaWQuXCIpO1xyXG4gICAgfVxyXG5cclxuICAgIGlmKCFQaW5uZWRVdGlscy5pc1Bpbm5hYmxlQ29sdW1uKGNvbEdyaWQpKSB7XHJcbiAgICAgIHRocm93IG5ldyBFcnJvcihcIkNvbHVtbiAnXCIgKyBjb2xHcmlkLm5hbWUgKyBcIicgY2Fubm90IGJlIHBpbm5lZC4gSXQgZWl0aGVyIHBpbm5lZCBvciBIVE1MLlwiKTtcclxuICAgIH1cclxuXHJcbiAgICAvL2xldCBuUm93TWluID0gZ3JpZC5taW5WaXNpYmxlUm93O1xyXG4gICAgLy9sZXQgblJvd01heCA9IGdyaWQubWF4VmlzaWJsZVJvdztcclxuICAgIC8vbGV0IG5Db2xNaW4gPSBncmlkLm1pblZpc2libGVDb2x1bW47XHJcbiAgICAvL2xldCBuQ29sTWF4ID0gZ3JpZC5tYXhWaXNpYmxlQ29sdW1uO1xyXG5cclxuXHJcbiAgICB0aGlzLm1fZkRldmljZVBpeGVsUmF0aW8gPSB3aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbztcclxuXHJcbiAgICBjb25zdCBkYXJ0ID0gREcudG9EYXJ0KGdyaWQpO1xyXG5cclxuICAgIGlmKGRhcnQubV9hclBpbm5lZENvbHMgPT09IHVuZGVmaW5lZClcclxuICAgICAgZGFydC5tX2FyUGlubmVkQ29scyA9IFtdO1xyXG5cclxuICAgIGlmKGRhcnQubV9hclBpbm5lZENvbHMubGVuZ3RoID09PSAwICYmICFHcmlkVXRpbHMuaXNSb3dIZWFkZXIoY29sR3JpZCkpIHtcclxuICAgICAgY29uc3QgY29sR3JpZDAgPSBncmlkLmNvbHVtbnMuYnlJbmRleCgwKTtcclxuICAgICAgaWYoY29sR3JpZDAgIT09IG51bGwgJiYgY29sR3JpZDAgIT09IHVuZGVmaW5lZClcclxuICAgICAgbmV3IFBpbm5lZENvbHVtbihjb2xHcmlkMCk7XHJcbiAgICB9XHJcblxyXG4gICAgY29uc3QgbldUb3RhbFBpbm5lZENvbHMgPSBQaW5uZWRVdGlscy5nZXRUb3RhbFBpbm5lZENvbHNXaWR0aChncmlkKTtcclxuICAgIGRhcnQubV9hclBpbm5lZENvbHMucHVzaCh0aGlzKTtcclxuXHJcbiAgICBjb25zdCB2aWV3VGFibGUgPSBncmlkLnZpZXc7XHJcbiAgICBjb25zdCBkZnJhbWUgPSBncmlkLmRhdGFGcmFtZTtcclxuXHJcbiAgICBjb25zdCBuVyA9IGNvbEdyaWQud2lkdGg7XHJcbiAgICB0aGlzLm1fY29sR3JpZCA9IGNvbEdyaWQ7XHJcbiAgICB0aGlzLm1fbldpZHRoQnVnID0gLTE7XHJcbiAgICB0cnkge1xyXG4gICAgICBjb2xHcmlkLnZpc2libGUgPSBmYWxzZTtcclxuICAgIH1cclxuICAgIGNhdGNoKGUpIHtcclxuICAgICAgLy9ERyBidWdcclxuICAgICAgY29uc29sZS5lcnJvcihcIkVSUk9SOiBDb3VsZG4ndCBoaWRlIGNvbHVtbiAnXCIgKyBjb2xHcmlkLm5hbWUgKyBcIicgZHVlIHRvIGEgREcgYnVnLiBBdHRlbXB0IHRvIHNldCB0aGUgd2lkdGggdG8gMFwiKTtcclxuICAgICAgdHJ5IHtcclxuICAgICAgICB0aGlzLm1fbldpZHRoQnVnID0gY29sR3JpZC53aWR0aDtcclxuICAgICAgICBjb2xHcmlkLndpZHRoID0gMDtcclxuICAgICAgfSBjYXRjaCAoZSkge1xyXG4gICAgICAgIC8vREcgYnVnXHJcbiAgICAgICAgY29uc29sZS5lcnJvcihcIkVSUk9SOiBDb3VsZG4ndCBzZXQgdGhlIHdpZHRoIHRvIDAgZm9yIGNvbHVtbiAnXCIgKyBjb2xHcmlkLm5hbWUgKyBcIicgZHVlIHRvIGEgREcgYnVnLiBUaGlzIGNvdWxkIGJlIGlnbm9yZWQgaWYgdGhlIGNvbHVtbiB2aXN1YWxseSBsb29rcyBvay5cIik7XHJcbiAgICAgIH1cclxuICAgIH1cclxuXHJcbiAgICBpZighR3JpZFV0aWxzLmlzUm93SGVhZGVyKGNvbEdyaWQpKSB7XHJcbiAgICAgIGlmIChjb2xHcmlkLnNldHRpbmdzID09PSBudWxsIHx8IGNvbEdyaWQuc2V0dGluZ3MgPT09IHVuZGVmaW5lZClcclxuICAgICAgICBjb2xHcmlkLnNldHRpbmdzID0ge307XHJcblxyXG4gICAgICBjb2xHcmlkLnNldHRpbmdzLmlzUGlubmVkID0gdHJ1ZTsgLy90aGlzIHdpbGwgYmUgc2F2ZWQgd2l0aCB0aGUgbGF5b3V0XHJcbiAgICAgIGNvbEdyaWQuc2V0dGluZ3MuaWR4UGlubmVkID0gZGFydC5tX2FyUGlubmVkQ29scy5sZW5ndGggLSAxO1xyXG4gICAgfVxyXG5cclxuICAgIGdyaWQuY2FudmFzLnN0eWxlLmxlZnQgPSAoZ3JpZC5jYW52YXMub2Zmc2V0TGVmdCArIG5XKS50b1N0cmluZygpICsgXCJweFwiO1xyXG4gICAgZ3JpZC5vdmVybGF5LnN0eWxlLmxlZnQ9IChncmlkLm92ZXJsYXkub2Zmc2V0TGVmdCArIG5XKS50b1N0cmluZygpICsgXCJweFwiO1xyXG5cclxuICAgIGdyaWQuY2FudmFzLnN0eWxlLndpZHRoID0gKGdyaWQuY2FudmFzLm9mZnNldFdpZHRoIC0gblcpLnRvU3RyaW5nKCkgKyBcInB4XCI7XHJcbiAgICBncmlkLm92ZXJsYXkuc3R5bGUud2lkdGg9IChncmlkLm92ZXJsYXkub2Zmc2V0V2lkdGggLSBuVykudG9TdHJpbmcoKSArIFwicHhcIjtcclxuXHJcbiAgICBjb25zdCBuSGVpZ2h0ID0gZ3JpZC5jYW52YXMuaGVpZ2h0Oy8vY2FudmFzIHBpeGVsIGhlaWdodFxyXG4gICAgY29uc3QgZUNhbnZhc1RoaXMgPSB1aS5jYW52YXMoblcqd2luZG93LmRldmljZVBpeGVsUmF0aW8sIG5IZWlnaHQpO1xyXG4gICAgY29uc3QgdGFiSW5kZXggPSAgZ3JpZC5jYW52YXMuZ2V0QXR0cmlidXRlKFwidGFiSW5kZXhcIik7XHJcbiAgICBpZih0YWJJbmRleCAhPT0gbnVsbClcclxuICAgICBlQ2FudmFzVGhpcy5zZXRBdHRyaWJ1dGUoXCJ0YWJJbmRleFwiLCB0YWJJbmRleCk7XHJcblxyXG4gICAgZUNhbnZhc1RoaXMuc3R5bGUucG9zaXRpb24gPSBcImFic29sdXRlXCI7XHJcbiAgICBlQ2FudmFzVGhpcy5zdHlsZS5sZWZ0ID0gbldUb3RhbFBpbm5lZENvbHMgKyBcInB4XCI7XHJcbiAgICBlQ2FudmFzVGhpcy5zdHlsZS50b3AgPSBncmlkLmNhbnZhcy5vZmZzZXRUb3AgKyBcInB4XCI7XHJcbiAgICBlQ2FudmFzVGhpcy5zdHlsZS53aWR0aCA9IG5XICsgXCJweFwiO1xyXG4gICAgZUNhbnZhc1RoaXMuc3R5bGUuaGVpZ2h0ID0gTWF0aC5yb3VuZChuSGVpZ2h0L3dpbmRvdy5kZXZpY2VQaXhlbFJhdGlvKSArIFwicHhcIjtcclxuXHJcbiAgICAvL2NvbnNvbGUubG9nKFwiaCBcIiArIGdyaWQuY2FudmFzLmhlaWdodCArIFwiIG9mZnNldCBcIiArIGdyaWQuY2FudmFzLm9mZnNldEhlaWdodCk7XHJcblxyXG4gICAgaWYoZ3JpZC5jYW52YXMucGFyZW50Tm9kZSA9PT0gbnVsbClcclxuICAgICAgdGhyb3cgbmV3IEVycm9yKFwiUGFyZW50IG5vZGUgZm9yIGNhbnZhcyBjYW5ub3QgYmUgbnVsbC5cIik7XHJcblxyXG4gICAgZ3JpZC5jYW52YXMucGFyZW50Tm9kZS5pbnNlcnRCZWZvcmUoZUNhbnZhc1RoaXMsIGdyaWQuY2FudmFzKTtcclxuICAgIHRoaXMubV9yb290ID0gZUNhbnZhc1RoaXM7XHJcblxyXG5cclxuICAgIGNvbnN0IGNvbEdyaWQwID0gZ3JpZC5jb2x1bW5zLmJ5SW5kZXgoMCk7XHJcbiAgICBpZihjb2xHcmlkMCAhPT0gbnVsbCAmJiBjb2xHcmlkMCAhPT0gdW5kZWZpbmVkKSB7Ly9ERyBCdWcgZnJvbSByZWFkaW5nIGxheW91dFxyXG4gICAgdHJ5e1xyXG4gICAgICAgIGNvbEdyaWQwLnZpc2libGUgPSBmYWxzZTtcclxuICAgICAgfVxyXG4gICAgICBjYXRjaChlKSB7XHJcbiAgICAgICAgY29uc29sZS5lcnJvcihcIkVSUk9SOiBDb3VsZG4ndCBoaWRlIHJvdyBoZWFkZXIuXCIpO1xyXG4gICAgICB9XHJcbiAgICB9XHJcblxyXG5cclxuICAgIC8vT25SZXNpemUgUm93IGhlYWRlclxyXG4gICAgY29uc3QgaGVhZGVyVGhpcyA9IHRoaXM7LypcclxuICAgIHRoaXMubV9vYnNlcnZlclJlc2l6ZSA9IG5ldyBSZXNpemVPYnNlcnZlcihlbnRyaWVzID0+IHtcclxuICAgICAgY29uc3QgZyA9IGhlYWRlclRoaXMubV9yb290LmdldENvbnRleHQoJzJkJyk7XHJcbiAgICAgIGZvciAobGV0IGVudHJ5IG9mIGVudHJpZXMpIHtcclxuICAgICAgICBoZWFkZXJUaGlzLnBhaW50KGcsIGdyaWQpO1xyXG4gICAgICB9XHJcbiAgICB9KTtcclxuICAgIHRoaXMubV9vYnNlcnZlclJlc2l6ZS5vYnNlcnZlKGhlYWRlclRoaXMubV9yb290KTsqL1xyXG5cclxuXHJcblxyXG4gICAgLy9PblJlc2l6ZSBHcmlkXHJcbiAgICB0aGlzLm1fb2JzZXJ2ZXJSZXNpemVHcmlkID0gbmV3IFJlc2l6ZU9ic2VydmVyKGZ1bmN0aW9uIChlbnRyaWVzIDogYW55KSB7XHJcblxyXG4gICAgICBjb25zdCBiQ3VycmVudCA9ICBERy50b0RhcnQoZ3Jvay5zaGVsbC52KSA9PT0gREcudG9EYXJ0KHZpZXdUYWJsZSk7XHJcbiAgICAgIGlmKCFiQ3VycmVudClcclxuICAgICAgICByZXR1cm47XHJcblxyXG4gICAgICBpZihoZWFkZXJUaGlzLm1fZkRldmljZVBpeGVsUmF0aW8gIT09IHdpbmRvdy5kZXZpY2VQaXhlbFJhdGlvIHx8IGdyaWQuY2FudmFzLmhlaWdodCAhPT0gZUNhbnZhc1RoaXMuaGVpZ2h0KSB7XHJcbiAgICAgICAgZUNhbnZhc1RoaXMud2lkdGggPSBuVyp3aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbztcclxuICAgICAgICBlQ2FudmFzVGhpcy5oZWlnaHQgPSBncmlkLmNhbnZhcy5oZWlnaHQ7XHJcbiAgICAgICAgZUNhbnZhc1RoaXMuc3R5bGUudG9wID0gZ3JpZC5jYW52YXMub2Zmc2V0VG9wICsgXCJweFwiO1xyXG4gICAgICAgIGVDYW52YXNUaGlzLnN0eWxlLndpZHRoID0gblcgKyBcInB4XCI7XHJcbiAgICAgICAgZUNhbnZhc1RoaXMuc3R5bGUuaGVpZ2h0ID0gTWF0aC5yb3VuZChncmlkLmNhbnZhcy5oZWlnaHQvd2luZG93LmRldmljZVBpeGVsUmF0aW8pICsgXCJweFwiO1xyXG5cclxuICAgICAgICBoZWFkZXJUaGlzLm1fZkRldmljZVBpeGVsUmF0aW8gPSB3aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbztcclxuICAgICAgfVxyXG5cclxuICAgICAgLy9jb25zb2xlLmxvZyhcIkdyaWQgUmVzaXplOiBcIiArIGdyaWQuY2FudmFzLmhlaWdodCArIFwiIFwiICsgd2luZG93LmRldmljZVBpeGVsUmF0aW8pO1xyXG4gICAgICAvL2VDYW52YXNUaGlzLnN0eWxlLmhlaWdodCA9IGdyaWQucm9vdC5zdHlsZS5oZWlnaHQ7XHJcbi8qXHJcbiAgICAgIGNvbnN0IGVDYW52YXNOZXcgPSB1aS5jYW52YXMoblcsIGdyaWQucm9vdC5vZmZzZXRIZWlnaHQpO1xyXG4gICAgICBpZihoZWFkZXJUaGlzLm1fcm9vdC5wYXJlbnROb2RlICE9PSBudWxsKSB7XHJcbiAgICAgICAgaGVhZGVyVGhpcy5tX3Jvb3QucGFyZW50Tm9kZS5yZXBsYWNlQ2hpbGQoZUNhbnZhc05ldywgaGVhZGVyVGhpcy5tX3Jvb3QpO1xyXG4gICAgICAgIGhlYWRlclRoaXMubV9yb290ID0gZUNhbnZhc05ldztcclxuICAgICAgfSovXHJcbiAgICAgIC8vaGVhZGVyVGhpcy5tX3Jvb3QuaGVpZ2h0ID0gZ3JpZC5yb290Lm9mZnNldEhlaWdodDtcclxuICAgICAgY29uc3QgZyA9IGVDYW52YXNUaGlzLmdldENvbnRleHQoJzJkJyk7XHJcbiAgICAgIGZvciAobGV0IGVudHJ5IG9mIGVudHJpZXMpIHtcclxuICAgICAgICBzZXRUaW1lb3V0KCgpPT4ge2hlYWRlclRoaXMucGFpbnQoZywgZ3JpZCk7fSwgMTAwKTtcclxuICAgICAgfVxyXG4gICAgfSk7XHJcblxyXG4gICAgdGhpcy5tX29ic2VydmVyUmVzaXplR3JpZD8ub2JzZXJ2ZShncmlkLmNhbnZhcyk7XHJcblxyXG5cclxuICAgIHRoaXMubV9oYW5kbGVyS2V5RG93biA9IHJ4anMuZnJvbUV2ZW50PEtleWJvYXJkRXZlbnQ+KGVDYW52YXNUaGlzLCAnbW91c2Vtb3ZlJykuc3Vic2NyaWJlKChlIDogS2V5Ym9hcmRFdmVudCkgPT4ge1xyXG5cclxuICAgICAgLy9hbGVydCgndXAnKTtcclxuICAgICAgc2V0VGltZW91dCgoKSA9PntcclxuICAgICAgICBjb25zdCBlZSA9IG5ldyBLZXlib2FyZEV2ZW50KGUudHlwZSwgZSk7XHJcbiAgICAgICAgdHJ5e2dyaWQub3ZlcmxheS5kaXNwYXRjaEV2ZW50KGVlKTt9XHJcbiAgICAgICAgY2F0Y2goZXgpIHtcclxuICAgICAgICAgIC8vY29uc29sZS5lcnJvcihleC5tZXNzYWdlKTtcclxuICAgICAgICB9XHJcbiAgICAgIH0sIDEpO1xyXG5cclxuICAgIH0pO1xyXG5cclxuXHJcbiAgICBjb25zdCBzY3JvbGxWZXJ0ID0gZ3JpZC52ZXJ0U2Nyb2xsO1xyXG4gICAgdGhpcy5tX2hhbmRsZXJWU2Nyb2xsID0gc2Nyb2xsVmVydC5vblZhbHVlc0NoYW5nZWQuc3Vic2NyaWJlKCgpID0+IHtcclxuICAgICAgY29uc3QgZyA9IGVDYW52YXNUaGlzLmdldENvbnRleHQoJzJkJyk7XHJcbiAgICAgIGhlYWRlclRoaXMucGFpbnQoZywgZ3JpZCk7XHJcbiAgICB9KTtcclxuXHJcbiAgICB0aGlzLm1faGFuZGxlclJvd3NGaWx0ZXJpbmcgPSBkZnJhbWUub25Sb3dzRmlsdGVyaW5nLnN1YnNjcmliZSgoKSA9PiB7XHJcbiAgICAgIHNldFRpbWVvdXQoKCkgPT4ge1xyXG4gICAgICAgIGNvbnN0IGcgPSBlQ2FudmFzVGhpcy5nZXRDb250ZXh0KCcyZCcpO1xyXG4gICAgICAgIGhlYWRlclRoaXMucGFpbnQoZywgZ3JpZCk7XHJcbiAgICAgIH0sIDEwMCk7XHJcblxyXG4gICAgfSk7XHJcblxyXG4gICAgdGhpcy5tX2hhbmRsZXJDdXJyUm93ID0gZGZyYW1lLm9uQ3VycmVudFJvd0NoYW5nZWQuc3Vic2NyaWJlKCgpID0+IHtcclxuICAgICAgICBjb25zdCBnID0gZUNhbnZhc1RoaXMuZ2V0Q29udGV4dCgnMmQnKTtcclxuICAgICAgICBoZWFkZXJUaGlzLnBhaW50KGcsIGdyaWQpO1xyXG4gICAgICB9XHJcbiAgICApO1xyXG5cclxuICAgIHRoaXMubV9oYW5kbGVyU2VsID0gZGZyYW1lLm9uU2VsZWN0aW9uQ2hhbmdlZC5zdWJzY3JpYmUoKGUgOiBhbnkpID0+IHtcclxuICAgICAgICBjb25zdCBnID0gZUNhbnZhc1RoaXMuZ2V0Q29udGV4dCgnMmQnKTtcclxuICAgICAgICBoZWFkZXJUaGlzLnBhaW50KGcsIGdyaWQpO1xyXG4gICAgICB9XHJcbiAgICApO1xyXG5cclxuICAgIHRoaXMubV9oYW5kbGVyQ29sc1JlbW92ZWQgPSBkZnJhbWUub25Db2x1bW5zUmVtb3ZlZC5zdWJzY3JpYmUoKGUgOiBDb2x1bW5zQXJncykgPT4ge1xyXG5cclxuICAgICAgICAgIGlmKGhlYWRlclRoaXMubV9jb2xHcmlkID09PSBudWxsKVxyXG4gICAgICAgICAgICByZXR1cm47XHJcbiAgICAgICAgICBmb3IobGV0IG5DPTA7IG5DPGUuY29sdW1ucy5sZW5ndGg7ICsrbkMpIHtcclxuICAgICAgICAgICAgaWYoZS5jb2x1bW5zW25DXS5uYW1lID09PSBoZWFkZXJUaGlzLm1fY29sR3JpZC5uYW1lKVxyXG4gICAgICAgICAgICAgIGhlYWRlclRoaXMuY2xvc2UoKTtcclxuICAgICAgICAgIH1cclxuICAgICAgICB9XHJcbiAgICApO1xyXG5cclxuICAgIHRoaXMubV9oYW5kbGVyQ29sTmFtZUNoYW5nZWQgPSBkZnJhbWUub25Db2x1bW5OYW1lQ2hhbmdlZC5zdWJzY3JpYmUoKGUgOiBhbnkpID0+IHtcclxuXHJcbiAgICAgICAgICBjb25zdCBkYXJ0ID0gdG9EYXJ0KGUpO1xyXG4gICAgICAgICAgY29uc3Qgc3RyQ29sTmFtZU9sZCA9IGRhcnQubmV3TmFtZTtcclxuICAgICAgICAgIGlmKHN0ckNvbE5hbWVPbGQgPT09IGhlYWRlclRoaXMubV9jb2xHcmlkPy5uYW1lKSB7XHJcbiAgICAgICAgICAgIGNvbnN0IGcgPSBlQ2FudmFzVGhpcy5nZXRDb250ZXh0KCcyZCcpO1xyXG4gICAgICAgICAgICBoZWFkZXJUaGlzLnBhaW50KGcsIGdyaWQpO1xyXG4gICAgICAgICAgfVxyXG4gICAgICAgIH1cclxuICAgICk7XHJcblxyXG5cclxuLypcclxuICAgIHRoaXMubV9oYW5kbGVyRmlsdGVyID0gZGZyYW1lLm9uUm93c0ZpbHRlcmVkLnN1YnNjcmliZSgoZSA6IGFueSkgPT4ge1xyXG4gICAgICAgIGNvbnN0IGcgPSBlQ2FudmFzVGhpcy5nZXRDb250ZXh0KCcyZCcpO1xyXG4gICAgICAgIGhlYWRlclRoaXMucGFpbnQoZywgZ3JpZCk7XHJcbiAgICAgIH1cclxuICAgICk7XHJcbiovXHJcblxyXG4gICAgdGhpcy5tX2hhbmRsZXJSb3dzUmVzaXplZCA9IGdyaWQub25Sb3dzUmVzaXplZC5zdWJzY3JpYmUoKGUgOiBhbnkpID0+IHtcclxuICAgICAgICBjb25zdCBnID0gZUNhbnZhc1RoaXMuZ2V0Q29udGV4dCgnMmQnKTtcclxuICAgICAgICBoZWFkZXJUaGlzLnBhaW50KGcsIGdyaWQpO1xyXG4gICAgICB9XHJcbiAgICApO1xyXG5cclxuICAgIHRoaXMubV9oYW5kbGVyUm93c1NvcnRlZCA9IGdyaWQub25Sb3dzU29ydGVkLnN1YnNjcmliZSgoZSA6IGFueSkgPT4ge1xyXG4gICAgICAgIGNvbnN0IGcgPSBlQ2FudmFzVGhpcy5nZXRDb250ZXh0KCcyZCcpO1xyXG4gICAgICAgIGhlYWRlclRoaXMucGFpbnQoZywgZ3JpZCk7XHJcbiAgICAgIH1cclxuICAgICk7XHJcbiAgfVxyXG5cclxuICBpc1Bpbm5lZCgpIDogYm9vbGVhbiB7XHJcbiAgICByZXR1cm4gdGhpcy5tX2NvbEdyaWQgIT09IG51bGw7XHJcbiAgfVxyXG5cclxuICBnZXRHcmlkQ29sdW1uKCkgOiBERy5HcmlkQ29sdW1uIHwgbnVsbHtcclxuICAgIHJldHVybiB0aGlzLm1fY29sR3JpZDtcclxuICB9XHJcblxyXG4gIGdldFdpZHRoKCkgOiBudW1iZXIge1xyXG4gICAgcmV0dXJuIHRoaXMubV9yb290ID09PSBudWxsID8gLTEgOiB0aGlzLm1fcm9vdC5vZmZzZXRXaWR0aDtcclxuICB9XHJcblxyXG4gIGdldFJvb3QoKSA6IEhUTUxDYW52YXNFbGVtZW50IHwgbnVsbCB7XHJcbiAgICByZXR1cm4gdGhpcy5tX3Jvb3Q7XHJcbiAgfVxyXG5cclxuICBwdWJsaWMgY2xvc2UoKSA6IHZvaWQge1xyXG5cclxuICAgIGlmKHRoaXMubV9jb2xHcmlkID09PSBudWxsKSB7XHJcbiAgICAgIHRocm93IG5ldyBFcnJvcihcIkNvbHVtbiBoYXMgYWxyZWFkeSBiZWVuIHVucGlubmVkXCIpO1xyXG4gICAgfVxyXG5cclxuICAgIGlmKHRoaXMubV9vYnNlcnZlclJlc2l6ZUdyaWQgIT09IG51bGwpIHtcclxuICAgICAgdGhpcy5tX29ic2VydmVyUmVzaXplR3JpZC5kaXNjb25uZWN0KCk7XHJcbiAgICAgIHRoaXMubV9vYnNlcnZlclJlc2l6ZUdyaWQgPSBudWxsO1xyXG4gICAgfVxyXG4vKm15IGNoYW5nZXNcclxuICAgIGlmKHRoaXMubV9vYnNlcnZlclJlc2l6ZSAhPT0gbnVsbCkge1xyXG4gICAgICB0aGlzLm1fb2JzZXJ2ZXJSZXNpemUuZGlzY29ubmVjdCgpO1xyXG4gICAgICB0aGlzLm1fb2JzZXJ2ZXJSZXNpemUgPSBudWxsO1xyXG4gICAgfVxyXG4gICAgKi9cclxuXHJcbiAgICB0aGlzLm1faGFuZGxlcktleURvd24/LnVuc3Vic2NyaWJlKCk7XHJcbiAgICB0aGlzLm1faGFuZGxlcktleURvd24gPSBudWxsO1xyXG5cclxuICAgIHRoaXMubV9oYW5kbGVyQ29sc1JlbW92ZWQ/LnVuc3Vic2NyaWJlKCk7XHJcbiAgICB0aGlzLm1faGFuZGxlckNvbHNSZW1vdmVkID0gbnVsbDtcclxuXHJcbiAgICB0aGlzLm1faGFuZGxlckNvbE5hbWVDaGFuZ2VkPy51bnN1YnNjcmliZSgpO1xyXG4gICAgdGhpcy5tX2hhbmRsZXJDb2xOYW1lQ2hhbmdlZCA9IG51bGw7XHJcblxyXG4gICAgdGhpcy5tX2hhbmRsZXJWU2Nyb2xsPy51bnN1YnNjcmliZSgpO1xyXG4gICAgdGhpcy5tX2hhbmRsZXJWU2Nyb2xsID0gbnVsbDtcclxuXHJcbiAgICB0aGlzLm1faGFuZGxlclJvd3NSZXNpemVkPy51bnN1YnNjcmliZSgpO1xyXG4gICAgdGhpcy5tX2hhbmRsZXJSb3dzUmVzaXplZCA9IG51bGw7XHJcblxyXG4gICAgdGhpcy5tX2hhbmRsZXJSb3dzU29ydGVkPy51bnN1YnNjcmliZSgpO1xyXG4gICAgdGhpcy5tX2hhbmRsZXJSb3dzU29ydGVkID0gbnVsbDtcclxuXHJcbiAgICB0aGlzLm1faGFuZGxlclJvd3NGaWx0ZXJpbmc/LnVuc3Vic2NyaWJlKCk7XHJcbiAgICB0aGlzLm1faGFuZGxlclJvd3NGaWx0ZXJpbmcgPSBudWxsO1xyXG5cclxuICAgIHRoaXMubV9oYW5kbGVyQ3VyclJvdz8udW5zdWJzY3JpYmUoKTtcclxuICAgIHRoaXMubV9oYW5kbGVyQ3VyclJvdyA9IG51bGw7XHJcblxyXG4gICAgdGhpcy5tX2hhbmRsZXJTZWw/LnVuc3Vic2NyaWJlKCk7XHJcbiAgICB0aGlzLm1faGFuZGxlclNlbCA9IG51bGw7XHJcblxyXG4gICAgY29uc3QgZ3JpZCA9IGdldEdyaWQodGhpcy5tX2NvbEdyaWQpO1xyXG4gICAgaWYoZ3JpZCA9PT0gbnVsbCl7XHJcbiAgICAgIHRocm93IG5ldyBFcnJvcihcIkNvbHVtbiAnXCIgKyB0aGlzLm1fY29sR3JpZC5uYW1lICsgXCInIGlzIGRpc2Nvbm5lY3RlZCBmcm9tIGdyaWQuXCIpO1xyXG4gICAgfVxyXG5cclxuICAgIGNvbnN0IGRhcnQgPSBERy50b0RhcnQoZ3JpZCk7XHJcbiAgICBjb25zdCBhciA9IGRhcnQubV9hclBpbm5lZENvbHM7XHJcbiAgICBjb25zdCBuSWR4ID0gYXIuaW5kZXhPZih0aGlzKTtcclxuICAgIGFyLnNwbGljZShuSWR4LCAxKTtcclxuXHJcbiAgICBpZih0aGlzLm1fcm9vdCA9PT0gbnVsbClcclxuICAgICAgdGhyb3cgbmV3IEVycm9yKCdSb290IGNhbm5vdCBiZSBudWxsJyk7XHJcblxyXG4gICAgbGV0IG5JZHhQaW5uZWQgPSAtMTtcclxuICAgIGxldCBjb2xHcmlkVG1wPSBudWxsO1xyXG4gICAgZm9yKGxldCBuPW5JZHg7IG48YXIubGVuZ3RoOyArK24pIHtcclxuICAgICAgY29sR3JpZFRtcCA9IGFyW25dO1xyXG4gICAgICBjb2xHcmlkVG1wLm1fcm9vdC5zdHlsZS5sZWZ0ID0gKGNvbEdyaWRUbXAubV9yb290Lm9mZnNldExlZnQgLSB0aGlzLm1fcm9vdC5vZmZzZXRXaWR0aCkudG9TdHJpbmcoKSArIFwicHhcIjtcclxuXHJcbiAgICAgIG5JZHhQaW5uZWQgPSAgY29sR3JpZFRtcC5tX2NvbEdyaWQuc2V0dGluZ3MuaWR4UGlubmVkO1xyXG4gICAgICBjb2xHcmlkVG1wLm1fY29sR3JpZC5zZXR0aW5ncy5pZHhQaW5uZWQgPSBuO1xyXG4gICAgfVxyXG5cclxuICAgIGlmKCFHcmlkVXRpbHMuaXNSb3dIZWFkZXIodGhpcy5tX2NvbEdyaWQpKSB7XHJcbiAgICAgIHRoaXMubV9jb2xHcmlkLnNldHRpbmdzLmlkeFBpbm5lZCA9IC0xO1xyXG4gICAgICB0aGlzLm1fY29sR3JpZC5zZXR0aW5ncy5pc1Bpbm5lZCA9IGZhbHNlO1xyXG4gICAgfVxyXG5cclxuXHJcbiAgICBpZih0aGlzLm1fbldpZHRoQnVnID49IDApIHtcclxuICAgICAgdHJ5IHtcclxuICAgICAgICB0aGlzLm1fY29sR3JpZC53aWR0aCA9IHRoaXMubV9uV2lkdGhCdWc7XHJcbiAgICAgIH1cclxuICAgICAgY2F0Y2goZSkge1xyXG4gICAgICAgIC8vREcgYnVnXHJcbiAgICAgICAgY29uc29sZS5lcnJvcihcIkVSUk9SOiBDb3VsZG4ndCBzZXQgdGhlIHdpZHRoIHRvIFwiICsgdGhpcy5tX25XaWR0aEJ1ZyArIFwiIGZvciBjb2x1bW4gJ1wiICsgdGhpcy5tX2NvbEdyaWQubmFtZSArIFwiJyBkdWUgdG8gYSBERyBidWcuIFRoaXMgY291bGQgYmUgaWdub3JlZCBpZiB0aGUgY29sdW1uIHZpc3VhbGx5IGxvb2tzIG9rLlwiKTtcclxuICAgICAgfVxyXG4gICAgfVxyXG5cclxuICAgIHRyeSB7XHJcbiAgICAgIHRoaXMubV9jb2xHcmlkLndpZHRoID0gdGhpcy5tX3Jvb3Qub2Zmc2V0V2lkdGg7XHJcbiAgICAgIHRoaXMubV9jb2xHcmlkLnZpc2libGUgPSB0cnVlO1xyXG4gICAgfVxyXG4gICAgY2F0Y2goZSkge1xyXG4gICAgICAvL0RHIGJ1Z1xyXG4gICAgICBjb25zb2xlLmVycm9yKFwiRVJST1I6IENvdWxkbid0IHNob3cgY29sdW1uICdcIiArIHRoaXMubV9jb2xHcmlkLm5hbWUgKyBcIicgZHVlIHRvIGEgREcgYnVnLiBUaGlzIGNvdWxkIGJlIGlnbm9yZWQgaWYgdGhlIGNvbHVtbiB2aXN1YWxseSBsb29rcyBvay5cIik7XHJcbiAgICB9XHJcblxyXG4gICAgZ3JpZC5jYW52YXMuc3R5bGUubGVmdCA9IChncmlkLmNhbnZhcy5vZmZzZXRMZWZ0IC0gdGhpcy5tX3Jvb3Qub2Zmc2V0V2lkdGgpLnRvU3RyaW5nKCkgKyBcInB4XCI7XHJcbiAgICBncmlkLm92ZXJsYXkuc3R5bGUubGVmdD0gKGdyaWQub3ZlcmxheS5vZmZzZXRMZWZ0IC0gdGhpcy5tX3Jvb3Qub2Zmc2V0V2lkdGgpLnRvU3RyaW5nKCkgKyBcInB4XCI7XHJcbiAgICBncmlkLmNhbnZhcy5zdHlsZS53aWR0aCA9IChncmlkLmNhbnZhcy5vZmZzZXRXaWR0aCArIHRoaXMubV9yb290Lm9mZnNldFdpZHRoKS50b1N0cmluZygpICsgXCJweFwiO1xyXG4gICAgZ3JpZC5vdmVybGF5LnN0eWxlLndpZHRoPSAoZ3JpZC5vdmVybGF5Lm9mZnNldFdpZHRoICsgdGhpcy5tX3Jvb3Qub2Zmc2V0V2lkdGgpLnRvU3RyaW5nKCkgKyBcInB4XCI7XHJcblxyXG4gICAgaWYodGhpcy5tX3Jvb3QucGFyZW50Tm9kZSAhPT0gbnVsbClcclxuICAgICB0aGlzLm1fcm9vdC5wYXJlbnROb2RlLnJlbW92ZUNoaWxkKHRoaXMubV9yb290KTtcclxuXHJcbiAgICB0aGlzLm1fcm9vdCA9IG51bGw7XHJcblxyXG4gICAgaWYgKGRhcnQubV9hclBpbm5lZENvbHMubGVuZ3RoID09PSAxICYmIGRhcnQubV9hclBpbm5lZENvbHNbMF0ubV9jb2xHcmlkLmlkeCA9PT0gMCAmJiB0aGlzLm1fY29sR3JpZC5pZHggIT09IDApIHtcclxuXHJcbiAgICAgICAgLy8gdHJ5e2NvbEdyaWQwLnZpc2libGUgPSB0cnVlO31cclxuICAgICAgICB0cnkge1xyXG4gICAgICAgICAgZGFydC5tX2FyUGlubmVkQ29sc1swXS5jbG9zZSgpO1xyXG4gICAgICAgIH0gY2F0Y2ggKGUpIHtcclxuICAgICAgICAgIGNvbnNvbGUuZXJyb3IoXCJFUlJPUjogQ291bGRuJ3QgY2xvc2UgcGlubmVkIGNvbHVtbiAnXCIgKyBkYXJ0Lm1fYXJQaW5uZWRDb2xzWzBdLm1fY29sR3JpZC5uYW1lICsgXCInIFwiKTtcclxuICAgICAgICB9XHJcbiAgICB9XHJcbiAgICB0aGlzLm1fY29sR3JpZCA9IG51bGw7XHJcbiAgfVxyXG5cclxuXHJcbiAgcHVibGljIG9uTW91c2VFbnRlcihlIDogTW91c2VFdmVudCkgOiB2b2lkIHtcclxuICAgIGlmKERFQlVHKVxyXG4gICAgICBjb25zb2xlLmxvZygnTW91c2UgRW50ZXIgUGlubmVkIENvbHVtbjogJyArIHRoaXMuZ2V0R3JpZENvbHVtbigpPy5uYW1lKTtcclxuICB9XHJcblxyXG4gIHB1YmxpYyBvbk1vdXNlTW92ZShlIDogTW91c2VFdmVudCkgOiB2b2lkIHtcclxuICAgIGlmKERFQlVHKVxyXG4gICAgICBjb25zb2xlLmxvZygnTW91c2UgTW92ZSBQaW5uZWQgQ29sdW1uOiAnICsgdGhpcy5nZXRHcmlkQ29sdW1uKCk/Lm5hbWUpO1xyXG5cclxuICAgIGlmKHRoaXMubV9jb2xHcmlkID09PSBudWxsIHx8IHRoaXMubV9yb290ID09PSBudWxsKVxyXG4gICAgICByZXR1cm47XHJcblxyXG4gICAgY29uc3QgZ3JpZCA9IHRoaXMubV9jb2xHcmlkLmdyaWQ7XHJcbiAgICBjb25zdCB2aWV3VGFibGUgPSBncmlkLnZpZXc7XHJcblxyXG4gICAgaWYoREcudG9EYXJ0KGdyb2suc2hlbGwudikgIT09IERHLnRvRGFydCh2aWV3VGFibGUpKSB7XHJcbiAgICAgIHJldHVybjtcclxuICAgIH1cclxuXHJcblxyXG4gICAgY29uc3QgYXJYWU9uQ2VsbCA9IFstMSwtMV07XHJcblxyXG4gICAgbGV0IG5Sb3dHcmlkID0gUGlubmVkQ29sdW1uLmhpdFRlc3RSb3dzKHRoaXMubV9yb290LCBncmlkLCBlLCBmYWxzZSwgYXJYWU9uQ2VsbCk7XHJcbiAgICBpZihuUm93R3JpZCA+PSAwKSB7XHJcbiAgICAgIGNvbnN0IGNlbGwgPSBncmlkLmNlbGwodGhpcy5tX2NvbEdyaWQubmFtZSwgblJvd0dyaWQpO1xyXG4gICAgICBjb25zdCByZW5kZXJlciA9IGdldFJlbmRlcmVyKGNlbGwpO1xyXG5cclxuICAgICAgaWYgKHJlbmRlcmVyIGluc3RhbmNlb2YgR3JpZENlbGxSZW5kZXJlckV4KSB7XHJcblxyXG4gICAgICAgIGlmICh0aGlzLm1fY2VsbEN1cnJlbnQgPT09IG51bGwpIHtcclxuICAgICAgICAgIHJlbmRlcmVyLm9uTW91c2VFbnRlckV4KGNlbGwsIGUsIGFyWFlPbkNlbGxbMF0sIGFyWFlPbkNlbGxbMV0pO1xyXG4gICAgICAgIH1cclxuXHJcbiAgICAgICAgaWYgKHRoaXMubV9jZWxsQ3VycmVudCAhPT0gbnVsbCAmJiBuUm93R3JpZCAhPT0gdGhpcy5tX2NlbGxDdXJyZW50LmdyaWRSb3cpIHtcclxuICAgICAgICAgIHJlbmRlcmVyLm9uTW91c2VMZWF2ZUV4KHRoaXMubV9jZWxsQ3VycmVudCwgZSwgLTEsIC0xKTtcclxuXHJcbiAgICAgICAgICByZW5kZXJlci5vbk1vdXNlRW50ZXJFeChjZWxsLCBlLCBhclhZT25DZWxsWzBdLCBhclhZT25DZWxsWzFdKTtcclxuICAgICAgICB9XHJcblxyXG4gICAgICAgIHJlbmRlcmVyLm9uTW91c2VNb3ZlRXgoY2VsbCwgZSwgYXJYWU9uQ2VsbFswXSwgYXJYWU9uQ2VsbFsxXSk7XHJcbiAgICAgIH1cclxuXHJcbiAgICAgIHRoaXMubV9jZWxsQ3VycmVudCA9IGNlbGw7XHJcbiAgICB9XHJcbiAgICBlbHNlIGlmICh0aGlzLm1fY2VsbEN1cnJlbnQgIT09IG51bGwpIHtcclxuICAgICAgY29uc3QgcmVuZGVyZXIgPSBnZXRSZW5kZXJlcih0aGlzLm1fY2VsbEN1cnJlbnQpO1xyXG4gICAgICBpZiAocmVuZGVyZXIgaW5zdGFuY2VvZiBHcmlkQ2VsbFJlbmRlcmVyRXgpIHtcclxuICAgICAgICByZW5kZXJlci5vbk1vdXNlTGVhdmVFeCh0aGlzLm1fY2VsbEN1cnJlbnQsIGUsIC0xLCAtMSk7XHJcbiAgICAgIH1cclxuXHJcbiAgICAgIHRoaXMubV9jZWxsQ3VycmVudCA9IG51bGw7XHJcbiAgICB9XHJcblxyXG4gICAgblJvd0dyaWQgPSBQaW5uZWRDb2x1bW4uaGl0VGVzdFJvd3ModGhpcy5tX3Jvb3QsIGdyaWQsIGUsIHRydWUsIHVuZGVmaW5lZCk7XHJcbiAgICBpZiAoblJvd0dyaWQgPj0gMCkge1xyXG4gICAgICB0aGlzLm1fblJlc2l6ZVJvd0dyaWRNb3ZpbmcgPSBuUm93R3JpZDtcclxuICAgICAgZG9jdW1lbnQuYm9keS5zdHlsZS5jdXJzb3IgPSBcInJvdy1yZXNpemVcIjtcclxuICAgICAgcmV0dXJuO1xyXG4gICAgfVxyXG5cclxuICAgIGlmKHRoaXMubV9uUmVzaXplUm93R3JpZE1vdmluZyA+PSAwKSB7XHJcbiAgICAgIHRoaXMubV9uUmVzaXplUm93R3JpZE1vdmluZyA9IC0xO1xyXG4gICAgICBkb2N1bWVudC5ib2R5LnN0eWxlLmN1cnNvciA9IFwiYXV0b1wiO1xyXG4gICAgfVxyXG5cclxuXHJcbiAgICAvL0hhbWJ1cmdlciBNZW51XHJcbiAgICBjb25zdCBjb2xHcmlkID0gdGhpcy5nZXRHcmlkQ29sdW1uKCk7XHJcbiAgICBpZihjb2xHcmlkID09PSBudWxsIHx8IGNvbEdyaWQubmFtZSA9PT0gJycpXHJcbiAgICAgIHJldHVybjtcclxuXHJcbiAgICBjb25zdCBlRGl2SGFtYiA9IEdyaWRVdGlscy5nZXRUb29sSWNvbkRpdihjb2xHcmlkLmdyaWQpO1xyXG4gICAgY29uc3QgbkhDb2xIZWFkZXIgPSBHcmlkVXRpbHMuZ2V0R3JpZENvbHVtbkhlYWRlckhlaWdodChjb2xHcmlkLmdyaWQpO1xyXG4gICAgaWYoMCA8PSBlLm9mZnNldFkgJiYgZS5vZmZzZXRZIDwgbkhDb2xIZWFkZXIpIHtcclxuXHJcbiAgICAgIGVEaXZIYW1iPy5zdHlsZS5yZW1vdmVQcm9wZXJ0eSgndmlzaWJpbGl0eScpO1xyXG4gICAgICBlRGl2SGFtYj8uc2V0QXR0cmlidXRlKCdjb2x1bW5fbmFtZScsIGNvbEdyaWQubmFtZSk7XHJcbiAgICAgIC8vY29uc29sZS5sb2coJ1Rvb2xzSWNvbiBmb3IgY29sdW1uICcgKyBjb2xHcmlkLm5hbWUpO1xyXG4gICAgICAvLyBAdHMtaWdub3JlXHJcbiAgICAgIGVEaXZIYW1iPy5zdHlsZS5sZWZ0ID0gKFBpbm5lZFV0aWxzLmdldFBpbm5lZENvbHVtbkxlZnQodGhpcykgKyB0aGlzLmdldFdpZHRoKCkgLSAxOCkgKyAncHgnO1xyXG4gICAgICAvLyBAdHMtaWdub3JlXHJcbiAgICAgIGVEaXZIYW1iPy5zdHlsZS50b3AgPSAoR3JpZFV0aWxzLmdldEdyaWRDb2x1bW5IZWFkZXJIZWlnaHQoY29sR3JpZC5ncmlkKSAtIDE2KSArIFwicHhcIjtcclxuICAgIH0gZWxzZSB7XHJcbiAgICAgIGNvbnN0IGNvbEdyaWQgPSB0aGlzLmdldEdyaWRDb2x1bW4oKTtcclxuICAgICAgaWYoY29sR3JpZCAhPSBudWxsKSB7XHJcbiAgICAgICAgICBlRGl2SGFtYj8uc2V0QXR0cmlidXRlKCdjb2x1bW5fbmFtZScsICcnKTtcclxuICAgICAgICAgIC8vIEB0cy1pZ25vcmVcclxuICAgICAgICAgIGVEaXZIYW1iPy5zdHlsZS52aXNpYmlsaXR5ID0gJ2hpZGRlbic7XHJcbiAgICAgICAgfVxyXG4gICAgfVxyXG4gIH1cclxuXHJcbiAgcHVibGljIG9uTW91c2VEcmFnKGUgOiBNb3VzZUV2ZW50KSA6IHZvaWQge1xyXG4gICAgaWYoREVCVUcpXHJcbiAgICAgY29uc29sZS5sb2coJ01vdXNlIERyYWcgUGlubmVkIENvbHVtbjogJyArIHRoaXMuZ2V0R3JpZENvbHVtbigpPy5uYW1lKTtcclxuXHJcbiAgICBpZih0aGlzLm1fY29sR3JpZCA9PT0gbnVsbCB8fCB0aGlzLm1fcm9vdCA9PT0gbnVsbClcclxuICAgIHJldHVybjtcclxuXHJcbiAgICBjb25zdCBncmlkID0gdGhpcy5tX2NvbEdyaWQuZ3JpZDtcclxuICAgIGNvbnN0IHZpZXdUYWJsZSA9IGdyaWQudmlldztcclxuXHJcbiAgICBpZihERy50b0RhcnQoZ3Jvay5zaGVsbC52KSAhPT0gREcudG9EYXJ0KHZpZXdUYWJsZSkpIHtcclxuICAgICAgcmV0dXJuO1xyXG4gICAgfVxyXG5cclxuICAgIGNvbnN0IGJSZXNpemluZyA9IHRoaXMubV9uUmVzaXplUm93R3JpZERyYWdnaW5nID49IDA7XHJcbiAgICBpZiAoYlJlc2l6aW5nKSB7XHJcblxyXG4gICAgICAvL2NvbnNvbGUubG9nKFwiRHJhZ2dpbmcgOiBcIiArIGhlYWRlclRoaXMubV9zdHJDb2xOYW1lKTtcclxuICAgICAgY29uc3QgbllEaWZmID0gZS5jbGllbnRZIC0gdGhpcy5tX25ZUmVzaXplRHJhZ2dpbmdBbmNob3I7XHJcbiAgICAgIGxldCBuSFJvd0dyaWQgPSB0aGlzLm1fbkhSZXNpemVSb3dzQmVmb3JlRHJhZyArIG5ZRGlmZjtcclxuXHJcbiAgICAgIGlmIChuSFJvd0dyaWQgPCBQaW5uZWRDb2x1bW4uTUlOX1JPV19IRUlHSFQpXHJcbiAgICAgICAgbkhSb3dHcmlkID0gUGlubmVkQ29sdW1uLk1JTl9ST1dfSEVJR0hUO1xyXG4gICAgICBlbHNlIGlmIChuSFJvd0dyaWQgPiBQaW5uZWRDb2x1bW4uTUFYX1JPV19IRUlHSFQpXHJcbiAgICAgICAgbkhSb3dHcmlkID0gUGlubmVkQ29sdW1uLk1BWF9ST1dfSEVJR0hUO1xyXG5cclxuICAgICAgY29uc3QgZUNhbnZhc1RoaXMgPSB0aGlzLm1fcm9vdDtcclxuXHJcbiAgICAgIGxldCBnID0gZUNhbnZhc1RoaXMuZ2V0Q29udGV4dCgnMmQnKTtcclxuICAgICAgaWYoZyA9PT0gbnVsbClcclxuICAgICAgICByZXR1cm47XHJcblxyXG4gICAgICBnLmZpbGxTdHlsZSA9IFwid2hpdGVcIjtcclxuICAgICAgY29uc3QgbkhIZWFkZXJDb2xzID0gR3JpZFV0aWxzLmdldEdyaWRDb2x1bW5IZWFkZXJIZWlnaHQoZ3JpZCk7XHJcbiAgICAgIGcuZmlsbFJlY3QoMCxuSEhlYWRlckNvbHMsIGVDYW52YXNUaGlzLm9mZnNldFdpZHRoLCBlQ2FudmFzVGhpcy5vZmZzZXRIZWlnaHQpO1xyXG5cclxuICAgICAgZ3JpZC5zZXRPcHRpb25zKHtcclxuICAgICAgICByb3dIZWlnaHQ6IG5IUm93R3JpZCAvL3RoaXMgd29uJ3QgdHJpZ2dlciBvblJvd3NSZXppemVkIGV2ZW50LCB3aGljaCBpcyBhIERHIGJ1Z1xyXG4gICAgICB9KTtcclxuXHJcbiAgICAgIG5vdGlmeUFsbFBpbm5lZENvbHNSb3dzUmVzaXplZCh0aGlzLCBuSFJvd0dyaWQsIHRydWUpO1xyXG4gICAgICBub3RpZnlBbGxDb2xzUm93c1Jlc2l6ZWQoZ3JpZCwgbkhSb3dHcmlkLCB0cnVlKTtcclxuXHJcbiAgICAgIGxldCBoZWFkZXIgPSBudWxsO1xyXG4gICAgICBjb25zdCBhciA9IGdyaWQuZGFydC5tX2FyUGlubmVkQ29scztcclxuICAgICAgZm9yKGxldCBuPTA7IG48YXIubGVuZ3RoOyArK24pIHtcclxuICAgICAgICBoZWFkZXIgPSBhcltuXTtcclxuICAgICAgICBnID0gaGVhZGVyLm1fcm9vdC5nZXRDb250ZXh0KCcyZCcpO1xyXG4gICAgICAgIGhlYWRlci5wYWludChnLCBncmlkKTtcclxuICAgICAgfVxyXG5cclxuICAgICAgdHJ5IHtcclxuICAgICAgICBjb25zdCBjb2xHcmlkMCA9IGdyaWQuY29sdW1ucy5ieUluZGV4KDApO1xyXG4gICAgICAgIGlmIChjb2xHcmlkMCAhPT0gbnVsbClcclxuICAgICAgICAgIGNvbEdyaWQwLnZpc2libGUgPSBmYWxzZTsvL3RlbXBvcmFyeSBhZGRyZXNzZWQgdGhlIERHIGJ1Z1xyXG4gICAgICB9XHJcbiAgICAgIGNhdGNoKGUpIHtcclxuICAgICAgICAvL0RHIGJ1Z1xyXG4gICAgICB9XHJcbiAgICAgIHJldHVybjtcclxuICAgIH1cclxuXHJcblxyXG4gIH1cclxuXHJcbiAgcHVibGljIG9uTW91c2VMZWF2ZShlIDogTW91c2VFdmVudCwgYk92ZXJsYXAgOiBib29sZWFuKSA6IHZvaWQge1xyXG4gICAgaWYoREVCVUcpXHJcbiAgICAgY29uc29sZS5sb2coJ01vdXNlIExlZnQgUGlubmVkIENvbHVtbjogJyArIHRoaXMuZ2V0R3JpZENvbHVtbigpPy5uYW1lICsgJyAgb3ZlcmxhcDogJyArIGJPdmVybGFwKTtcclxuXHJcbiAgICBpZih0aGlzLm1fblJlc2l6ZVJvd0dyaWRNb3ZpbmcgPj0gMCkge1xyXG4gICAgICB0aGlzLm1fblJlc2l6ZVJvd0dyaWRNb3ZpbmcgPSAtMTtcclxuICAgICAgZG9jdW1lbnQuYm9keS5zdHlsZS5jdXJzb3IgPSBcImF1dG9cIjtcclxuICAgIH1cclxuXHJcbiAgICBpZih0aGlzLm1fY2VsbEN1cnJlbnQgIT09IG51bGwpIHtcclxuICAgICAgY29uc3QgcmVuZGVyZXIgPSBnZXRSZW5kZXJlcih0aGlzLm1fY2VsbEN1cnJlbnQpO1xyXG4gICAgICBpZiAocmVuZGVyZXIgaW5zdGFuY2VvZiBHcmlkQ2VsbFJlbmRlcmVyRXgpIHtcclxuICAgICAgICBjb25zdCBlTW91c2UgPSBlIGFzIE1vdXNlRXZlbnQ7XHJcbiAgICAgICAgcmVuZGVyZXIub25Nb3VzZUxlYXZlRXgodGhpcy5tX2NlbGxDdXJyZW50LCBlTW91c2UsIC0xLCAtMSk7XHJcbiAgICAgIH1cclxuICAgICAgdGhpcy5tX2NlbGxDdXJyZW50ID0gbnVsbDtcclxuICAgIH1cclxuXHJcbiAgICBjb25zdCBjb2xHcmlkID0gdGhpcy5nZXRHcmlkQ29sdW1uKCk7XHJcbiAgICBpZihjb2xHcmlkICE9IG51bGwgJiYgIWJPdmVybGFwKSB7XHJcbiAgICAgIGNvbnN0IGVEaXZIYW1iID0gR3JpZFV0aWxzLmdldFRvb2xJY29uRGl2KGNvbEdyaWQuZ3JpZCk7XHJcbiAgICAgIGVEaXZIYW1iPy5zZXRBdHRyaWJ1dGUoJ2NvbHVtbl9uYW1lJywgJycpO1xyXG4gICAgICAvLyBAdHMtaWdub3JlXHJcbiAgICAgIGVEaXZIYW1iPy5zdHlsZS52aXNpYmlsaXR5ID0gJ2hpZGRlbic7XHJcbiAgICB9XHJcblxyXG5cclxuICB9XHJcblxyXG4gIHB1YmxpYyBvbk1vdXNlRGJsQ2xpY2soZSA6IE1vdXNlRXZlbnQpIDogdm9pZCB7XHJcbiAgICBpZihERUJVRylcclxuICAgICBjb25zb2xlLmxvZygnTW91c2UgRGJsIENsaWNrZWQgUGlubmVkIENvbHVtbjogJyArIHRoaXMuZ2V0R3JpZENvbHVtbigpPy5uYW1lKTtcclxuXHJcbiAgICBpZih0aGlzLm1fY29sR3JpZCA9PT0gbnVsbCB8fCB0aGlzLm1fcm9vdCA9PT0gbnVsbClcclxuICAgICAgcmV0dXJuO1xyXG5cclxuICAgIGNvbnN0IGdyaWQgPSB0aGlzLm1fY29sR3JpZC5ncmlkO1xyXG4gICAgY29uc3Qgdmlld1RhYmxlID0gZ3JpZD8udmlldztcclxuXHJcbiAgICBpZiAoREcudG9EYXJ0KGdyb2suc2hlbGwudikgIT09IERHLnRvRGFydCh2aWV3VGFibGUpKSB7XHJcbiAgICAgIHJldHVybjtcclxuICAgIH1cclxuXHJcbiAgICBpZih0aGlzLm1fY29sR3JpZD8ubmFtZSA9PT0gJycpXHJcbiAgICAgIHJldHVybjtcclxuXHJcbiAgICBpZih0aGlzLm1fYlNvcnRlZEFzY2VuZGluZyA9PSBudWxsKVxyXG4gICAgICB0aGlzLm1fYlNvcnRlZEFzY2VuZGluZyA9IHRydWU7XHJcbiAgICBlbHNlIGlmKHRoaXMubV9iU29ydGVkQXNjZW5kaW5nKVxyXG4gICAgICB0aGlzLm1fYlNvcnRlZEFzY2VuZGluZyA9IGZhbHNlO1xyXG4gICAgZWxzZSB0aGlzLm1fYlNvcnRlZEFzY2VuZGluZyA9IHRydWU7XHJcblxyXG4gICAgY29uc3QgbkhIZWFkZXJDb2xzID0gR3JpZFV0aWxzLmdldEdyaWRDb2x1bW5IZWFkZXJIZWlnaHQoZ3JpZCk7XHJcblxyXG4gICAgaWYoMCA8PSBlLm9mZnNldFggJiYgZS5vZmZzZXRYIDw9IHRoaXMubV9yb290Lm9mZnNldFdpZHRoICYmXHJcbiAgICAgICAgMCA8PSBlLm9mZnNldFkgJiYgZS5vZmZzZXRZIDw9IG5ISGVhZGVyQ29scykgICAvL29uIHRoZSByb3dzIGhlYWRlclxyXG4gICAge1xyXG4gICAgICBncmlkPy5zb3J0KFt0aGlzLm1fY29sR3JpZD8ubmFtZV0sIFt0aGlzLm1fYlNvcnRlZEFzY2VuZGluZ10pO1xyXG4gICAgfVxyXG4gIH1cclxuXHJcbiAgcHVibGljIG9uTW91c2VEb3duKGUgOiBNb3VzZUV2ZW50KSA6IHZvaWQge1xyXG4gICAgaWYoREVCVUcpXHJcbiAgICAgY29uc29sZS5sb2coJ01vdXNlIERvd24gUGlubmVkIENvbHVtbjogJyArIHRoaXMuZ2V0R3JpZENvbHVtbigpPy5uYW1lKTtcclxuLypcclxuICAgIGlmKGUudmlldyAhPSBudWxsKSB7XHJcbiAgICAgIGNvbnN0IGVlID0gZG9jdW1lbnQuY3JlYXRlRXZlbnQoIFwiTW91c2VFdmVudFwiICk7XHJcbiAgICAgIGVlLmluaXRNb3VzZUV2ZW50KGUudHlwZSwgZS5idWJibGVzLCBlLmNhbmNlbGFibGUsIGUudmlldywgZS5kZXRhaWwsIGUuc2NyZWVuWCArIDEwMCwgZS5zY3JlZW5ZLCBlLmNsaWVudFggKyAxMDAsIGUuY2xpZW50WSwgZS5jdHJsS2V5LCBlLmFsdEtleSwgZS5zaGlmdEtleSwgZS5tZXRhS2V5LCBlLmJ1dHRvbiwgZS5yZWxhdGVkVGFyZ2V0KTtcclxuICAgICAgdGhpcy5tX2NvbEdyaWQ/LmdyaWQucm9vdC5kaXNwYXRjaEV2ZW50KGVlKTtcclxuICAgICAgcmV0dXJuO1xyXG4gICAgfVxyXG4qL1xyXG5cclxuICAgIGlmKHRoaXMubV9jb2xHcmlkID09PSBudWxsKVxyXG4gICAgICByZXR1cm47XHJcblxyXG4gICAgY29uc3QgZ3JpZCA9IHRoaXMubV9jb2xHcmlkPy5ncmlkO1xyXG4gICAgY29uc3Qgdmlld1RhYmxlID0gZ3JpZD8udmlldztcclxuICAgIGlmKERHLnRvRGFydChncm9rLnNoZWxsLnYpICE9PSBERy50b0RhcnQodmlld1RhYmxlKSlcclxuICAgICAgcmV0dXJuO1xyXG5cclxuICAgIGlmKGUuYnV0dG9ucyAhPT0gMSlcclxuICAgICAgcmV0dXJuO1xyXG5cclxuICAgIGxldCBlQ2FudmFzVGhpcyA9IHRoaXMubV9yb290O1xyXG4gICAgaWYoZUNhbnZhc1RoaXMgPT09IG51bGwpXHJcbiAgICAgIHJldHVybjtcclxuXHJcbiAgICB0aGlzLm1fblJlc2l6ZVJvd0dyaWRNb3ZpbmcgPSAtMTtcclxuICAgIGNvbnN0IGJBZGRUb1NlbCA6IGJvb2xlYW4gPSBlLmN0cmxLZXkgfHwgZS5zaGlmdEtleTtcclxuXHJcbiAgICBsZXQgblJvd0dyaWQgPSBiQWRkVG9TZWwgPyAtMSA6IFBpbm5lZENvbHVtbi5oaXRUZXN0Um93cyhlQ2FudmFzVGhpcywgZ3JpZCwgZSwgdHJ1ZSwgdW5kZWZpbmVkKTtcclxuICAgIGlmIChuUm93R3JpZCA+PSAwKSB7XHJcbiAgICAgIGNvbnN0IG5IUm93cyA9IEdyaWRVdGlscy5nZXRHcmlkUm93SGVpZ2h0KGdyaWQpO1xyXG4gICAgICB0aGlzLm1fblJlc2l6ZVJvd0dyaWREcmFnZ2luZyA9IG5Sb3dHcmlkO1xyXG4gICAgICB0aGlzLm1fbllSZXNpemVEcmFnZ2luZ0FuY2hvciA9IGUuY2xpZW50WTtcclxuICAgICAgdGhpcy5tX25IUmVzaXplUm93c0JlZm9yZURyYWcgPSBuSFJvd3M7XHJcbiAgICB9XHJcbiAgICBlbHNlXHJcbiAgICB7XHJcblxyXG4gICAgICBuUm93R3JpZCA9IFBpbm5lZENvbHVtbi5oaXRUZXN0Um93cyhlQ2FudmFzVGhpcywgZ3JpZCwgZSwgZmFsc2UsIHRoaXMubV9hclhZTW91c2VPbkNlbGxEb3duKTtcclxuXHJcbiAgICAgIHRoaXMubV9uUm93R3JpZERyYWdnaW5nID0gblJvd0dyaWQ7XHJcbiAgICAgIHRoaXMubV9uWURyYWdnaW5nQW5jaG9yID0gZS5jbGllbnRZO1xyXG5cclxuICAgICAgY29uc3QgY2VsbCA9IGdyaWQuY2VsbCh0aGlzLm1fY29sR3JpZC5uYW1lLCBuUm93R3JpZCk7XHJcbiAgICAgIGNvbnN0IHJlbmRlcmVyID0gZ2V0UmVuZGVyZXIoY2VsbCk7XHJcbiAgICAgIGlmKHJlbmRlcmVyIGluc3RhbmNlb2YgR3JpZENlbGxSZW5kZXJlckV4KSB7XHJcbiAgICAgICAgcmVuZGVyZXIub25Nb3VzZURvd25FeChjZWxsLCBlLCB0aGlzLm1fYXJYWU1vdXNlT25DZWxsRG93blswXSwgdGhpcy5tX2FyWFlNb3VzZU9uQ2VsbERvd25bMV0pO1xyXG4gICAgICB9XHJcbiAgICB9XHJcbiAgfVxyXG5cclxuICBwdWJsaWMgb25Nb3VzZVVwKGUgOiBNb3VzZUV2ZW50KSA6IHZvaWQge1xyXG4gICAgaWYoREVCVUcpXHJcbiAgICAgY29uc29sZS5sb2coJ01vdXNlIFVwIFBpbm5lZCBDb2x1bW46ICcgKyB0aGlzLmdldEdyaWRDb2x1bW4oKT8ubmFtZSk7XHJcbi8qXHJcbiAgICBpZihlLnZpZXcgIT0gbnVsbCkge1xyXG4gICAgICBjb25zdCBlZSA9IGRvY3VtZW50LmNyZWF0ZUV2ZW50KCBcIk1vdXNlRXZlbnRcIiApO1xyXG4gICAgICBlZS5pbml0TW91c2VFdmVudChlLnR5cGUsIGUuYnViYmxlcywgZS5jYW5jZWxhYmxlLCBlLnZpZXcsIGUuZGV0YWlsLCBlLnNjcmVlblggKyAxMDAsIGUuc2NyZWVuWSwgZS5jbGllbnRYICsgMTAwLCBlLmNsaWVudFksIGUuY3RybEtleSwgZS5hbHRLZXksIGUuc2hpZnRLZXksIGUubWV0YUtleSwgZS5idXR0b24sIGUucmVsYXRlZFRhcmdldCk7XHJcbiAgICAgIHRoaXMubV9jb2xHcmlkPy5ncmlkLnJvb3QuZGlzcGF0Y2hFdmVudChlZSk7XHJcbiAgICAgIHJldHVybjtcclxuICAgIH1cclxuKi9cclxuXHJcbiAgICBpZih0aGlzLm1fY29sR3JpZCA9PT0gbnVsbCB8fCB0aGlzLm1fcm9vdCA9PSBudWxsKVxyXG4gICAgICByZXR1cm47XHJcblxyXG4gICAgY29uc3QgZ3JpZCA9IHRoaXMubV9jb2xHcmlkPy5ncmlkO1xyXG4gICAgY29uc3Qgdmlld1RhYmxlID0gZ3JpZD8udmlldztcclxuXHJcbiAgICBpZihERy50b0RhcnQoZ3Jvay5zaGVsbC52KSAhPT0gREcudG9EYXJ0KHZpZXdUYWJsZSkpIHtcclxuICAgICAgcmV0dXJuO1xyXG4gICAgfVxyXG4vKlxyXG4gICBpZihlLmJ1dHRvbiA9PT0gMikge1xyXG5cclxuICAgICBjb25zdCBlRGl2UE9wdXAgOiBIVE1MRWxlbWVudCB8IG51bGwgPSBHcmlkVXRpbHMuZ2V0R3JpZERhcnRQb3B1cE1lbnUoKTtcclxuICAgICBlRGl2UE9wdXA/LnNldEF0dHJpYnV0ZSgnY29sdW1uX25hbWUnLCB0aGlzLm1fY29sR3JpZC5uYW1lKTtcclxuICAgICBsZXQgZCA9IDA7XHJcbiAgICAgcmV0dXJuO1xyXG4gICB9Ki9cclxuXHJcblxyXG4gICAgaWYodGhpcy5tX25SZXNpemVSb3dHcmlkRHJhZ2dpbmcgPj0gMCkge1xyXG4gICAgICBjb25zdCBuSFJvdyA9IEdyaWRVdGlscy5nZXRHcmlkUm93SGVpZ2h0KGdyaWQpO1xyXG4gICAgICBub3RpZnlBbGxQaW5uZWRDb2xzUm93c1Jlc2l6ZWQodGhpcywgbkhSb3csIGZhbHNlKTtcclxuICAgICAgbm90aWZ5QWxsQ29sc1Jvd3NSZXNpemVkKGdyaWQsIG5IUm93LCBmYWxzZSk7XHJcbiAgICB9XHJcblxyXG4gICAgdGhpcy5tX25IUmVzaXplUm93c0JlZm9yZURyYWcgPSAtMTtcclxuICAgIHRoaXMubV9uUmVzaXplUm93R3JpZERyYWdnaW5nID0gLTE7XHJcbiAgICB0aGlzLm1fbllSZXNpemVEcmFnZ2luZ0FuY2hvciA9IC0xO1xyXG4gICAgdGhpcy5tX25SZXNpemVSb3dHcmlkTW92aW5nID0gLTE7XHJcblxyXG4gICAgZG9jdW1lbnQuYm9keS5zdHlsZS5jdXJzb3IgPSBcImF1dG9cIjtcclxuXHJcbiAgICBpZih0aGlzLm1fblJvd0dyaWREcmFnZ2luZyA+PSAwKSB7XHJcbiAgICAgIGNvbnN0IGRmcmFtZSA9IGdyaWQuZGF0YUZyYW1lO1xyXG4gICAgICBjb25zdCBiQWRkVG9TZWwgPSBlLmN0cmxLZXk7XHJcbiAgICAgIGNvbnN0IGJSYW5nZVNlbCA9IGUuc2hpZnRLZXk7XHJcblxyXG4gICAgICBjb25zdCBuUm93R3JpZCA9IFBpbm5lZENvbHVtbi5oaXRUZXN0Um93cyh0aGlzLm1fcm9vdCwgZ3JpZCwgZSwgZmFsc2UsIHRoaXMubV9hclhZTW91c2VPbkNlbGxVcCk7XHJcbiAgICAgIGlmKCFiQWRkVG9TZWwgJiYgIWJSYW5nZVNlbCAmJiBuUm93R3JpZCA9PT0gdGhpcy5tX25Sb3dHcmlkRHJhZ2dpbmcpIHsgLy9jbGljayBvbiB0aGUgc2FtZSByb3cgd2hpY2ggd2lsbCBiZWNvbWUgYWN0aXZlXHJcblxyXG4gICAgICAgIGxldCBjZWxsUkggPSBudWxsO1xyXG4gICAgICAgIHRyeSB7XHJcbiAgICAgICAgICBjZWxsUkggPSBncmlkLmNlbGwoXCJcIiwgblJvd0dyaWQpO1xyXG4gICAgICAgIH1cclxuICAgICAgICBjYXRjaChlKSB7XHJcbiAgICAgICAgICBsZXQgY29sRyA9IG51bGw7XHJcbiAgICAgICAgICBjb25zdCBsc3RDb2xzID0gZ3JpZC5jb2x1bW5zO1xyXG4gICAgICAgICAgZm9yKGxldCBuQz0xOyBuQzxsc3RDb2xzLmxlbmd0aDsgKytuQykge1xyXG4gICAgICAgICAgICBjb2xHID0gbHN0Q29scy5ieUluZGV4KG5DKTtcclxuICAgICAgICAgICAgY2VsbFJIID0gY29sRyA9PT0gbnVsbCA/IG51bGwgOiBncmlkLmNlbGwoY29sRy5uYW1lLCBuUm93R3JpZCk7XHJcbiAgICAgICAgICAgIGlmKGNlbGxSSCAhPT0gbnVsbClcclxuICAgICAgICAgICAgICBicmVhaztcclxuICAgICAgICAgIH1cclxuICAgICAgICB9XHJcbiAgICAgICAgaWYoY2VsbFJIICE9PSBudWxsKSB7XHJcbiAgICAgICAgICBjb25zdCBuUm93VGFibGUgOiBhbnkgPSBjZWxsUkgudGFibGVSb3dJbmRleDtcclxuICAgICAgICAgIGlmKG5Sb3dUYWJsZSAhPT0gbnVsbClcclxuICAgICAgICAgICAgZGZyYW1lLmN1cnJlbnRSb3cgPSBuUm93VGFibGU7XHJcbiAgICAgICAgfVxyXG4gICAgICB9XHJcbiAgICAgIGVsc2VcclxuICAgICAge1xyXG4gICAgICAgIGNvbnN0IGJpdHNldFNlbCA9IGRmcmFtZS5zZWxlY3Rpb247XHJcblxyXG4gICAgICAgIGlmKCFiQWRkVG9TZWwgfHwgYlJhbmdlU2VsKVxyXG4gICAgICAgICAgYml0c2V0U2VsLnNldEFsbChmYWxzZSwgdHJ1ZSk7XHJcblxyXG4gICAgICAgIGxldCBuUm93TWluID0gdGhpcy5tX25Sb3dHcmlkRHJhZ2dpbmcgPCBuUm93R3JpZCA/IHRoaXMubV9uUm93R3JpZERyYWdnaW5nIDogblJvd0dyaWQ7XHJcbiAgICAgICAgbGV0IG5Sb3dNYXggPSB0aGlzLm1fblJvd0dyaWREcmFnZ2luZyA+IG5Sb3dHcmlkID8gdGhpcy5tX25Sb3dHcmlkRHJhZ2dpbmcgOiBuUm93R3JpZDtcclxuXHJcbiAgICAgICAgaWYoYlJhbmdlU2VsKSB7XHJcbiAgICAgICAgICBsZXQgblJvd0dyaWRBY3RpdmUgPSBHcmlkVXRpbHMuZ2V0QWN0aXZlR3JpZFJvdyhncmlkKTtcclxuICAgICAgICAgIGlmKG5Sb3dHcmlkQWN0aXZlID09PSBudWxsKVxyXG4gICAgICAgICAgICBuUm93R3JpZEFjdGl2ZSA9IDA7XHJcblxyXG4gICAgICAgICAgblJvd01pbiA9IG5Sb3dHcmlkQWN0aXZlIDwgblJvd0dyaWQgPyBuUm93R3JpZEFjdGl2ZSA6IG5Sb3dHcmlkO1xyXG4gICAgICAgICAgblJvd01heCA9IG5Sb3dHcmlkQWN0aXZlID4gblJvd0dyaWQgPyBuUm93R3JpZEFjdGl2ZSA6IG5Sb3dHcmlkO1xyXG4gICAgICAgIH1cclxuXHJcblxyXG4gICAgICAgIGxldCBjZWxsUkggPSBudWxsO1xyXG4gICAgICAgIGxldCBuUm93VGFibGUgPSAtMTtcclxuICAgICAgICBmb3IobGV0IG5Sb3c9blJvd01pbjsgblJvdzw9blJvd01heDsgKytuUm93KSB7XHJcblxyXG4gICAgICAgICAgdHJ5IHtcclxuICAgICAgICAgICAgY2VsbFJIID0gZ3JpZC5jZWxsKFwiXCIsIG5Sb3cpO1xyXG4gICAgICAgICAgfVxyXG4gICAgICAgICAgY2F0Y2goZSkge1xyXG4gICAgICAgICAgICBsZXQgY29sRyA9IG51bGw7XHJcbiAgICAgICAgICAgIGNvbnN0IGxzdENvbHMgPSBncmlkLmNvbHVtbnM7XHJcbiAgICAgICAgICAgIGZvcihsZXQgbkM9MTsgbkM8bHN0Q29scy5sZW5ndGg7ICsrbkMpIHtcclxuICAgICAgICAgICAgICBjb2xHID0gbHN0Q29scy5ieUluZGV4KG5DKTtcclxuICAgICAgICAgICAgICBjZWxsUkggPSBjb2xHID09PSBudWxsID8gbnVsbCA6IGdyaWQuY2VsbChjb2xHLm5hbWUsIG5Sb3dHcmlkKTtcclxuICAgICAgICAgICAgICBpZihjZWxsUkggIT09IG51bGwpXHJcbiAgICAgICAgICAgICAgICBicmVhaztcclxuICAgICAgICAgICAgfVxyXG4gICAgICAgICAgfVxyXG5cclxuICAgICAgICAgIGlmKGNlbGxSSCAhPT0gbnVsbCAmJiBjZWxsUkgudGFibGVSb3dJbmRleCAhPT0gbnVsbCkge1xyXG4gICAgICAgICAgICBuUm93VGFibGUgPSBjZWxsUkgudGFibGVSb3dJbmRleDtcclxuICAgICAgICAgICAgYml0c2V0U2VsLnNldChuUm93VGFibGUsIHRydWUsIHRydWUpO1xyXG4gICAgICAgICAgfVxyXG4gICAgICAgIH1cclxuICAgICAgfVxyXG5cclxuICAgICAgY29uc3QgY2VsbCA9IGdyaWQuY2VsbCh0aGlzLm1fY29sR3JpZC5uYW1lLCBuUm93R3JpZCk7XHJcbiAgICAgIGNvbnN0IHJlbmRlcmVyID0gZ2V0UmVuZGVyZXIoY2VsbCk7XHJcbiAgICAgIGlmKHJlbmRlcmVyIGluc3RhbmNlb2YgR3JpZENlbGxSZW5kZXJlckV4KSB7XHJcbiAgICAgICAgcmVuZGVyZXIub25Nb3VzZVVwRXgoY2VsbCwgZSwgdGhpcy5tX2FyWFlNb3VzZU9uQ2VsbFVwWzBdLCB0aGlzLm1fYXJYWU1vdXNlT25DZWxsVXBbMV0pO1xyXG4gICAgICB9XHJcblxyXG4gICAgICBpZih0aGlzLm1fYXJYWU1vdXNlT25DZWxsVXBbMF0gPT09IHRoaXMubV9hclhZTW91c2VPbkNlbGxEb3duWzBdICYmIHRoaXMubV9hclhZTW91c2VPbkNlbGxEb3duWzFdID09PSB0aGlzLm1fYXJYWU1vdXNlT25DZWxsVXBbMV0pIHtcclxuICAgICAgICBpZihyZW5kZXJlciBpbnN0YW5jZW9mIEdyaWRDZWxsUmVuZGVyZXJFeCkge1xyXG4gICAgICAgICAgcmVuZGVyZXIub25DbGlja0V4KGNlbGwsIGUsIHRoaXMubV9hclhZTW91c2VPbkNlbGxVcFswXSwgdGhpcy5tX2FyWFlNb3VzZU9uQ2VsbFVwWzFdKTtcclxuICAgICAgICB9XHJcbiAgICAgIH1cclxuXHJcbiAgICAgIHRoaXMubV9uUm93R3JpZERyYWdnaW5nID0gLTE7XHJcbiAgICAgIHRoaXMubV9uWURyYWdnaW5nQW5jaG9yID0gLTE7XHJcbiAgICAgIHRoaXMubV9hclhZTW91c2VPbkNlbGxEb3duWzBdID0gLTI7XHJcbiAgICAgIHRoaXMubV9hclhZTW91c2VPbkNlbGxEb3duWzFdID0gLTI7XHJcbiAgICAgIHRoaXMubV9hclhZTW91c2VPbkNlbGxVcFswXSA9IC0xO1xyXG4gICAgICB0aGlzLm1fYXJYWU1vdXNlT25DZWxsVXBbMV0gPSAtMTtcclxuICAgIH1cclxuICB9XHJcblxyXG4gIHB1YmxpYyBvbkNvbnRleHRNZW51KGUgOiBNb3VzZUV2ZW50KSA6IHZvaWQge1xyXG4gICBpZihERUJVRylcclxuICAgIGNvbnNvbGUubG9nKCdDb250ZXh0IG1lbnUgUGlubmVkIENvbHVtbjogJyArIHRoaXMuZ2V0R3JpZENvbHVtbigpPy5uYW1lKTtcclxuICB9XHJcblxyXG4gIHB1YmxpYyBvbk1vdXNlV2hlZWwoZSA6IFdoZWVsRXZlbnQpIDogdm9pZCB7XHJcblxyXG4gICAgaWYodGhpcy5tX2NvbEdyaWQgPT09IG51bGwgfHwgdGhpcy5tX3Jvb3QgPT0gbnVsbClcclxuICAgICAgcmV0dXJuO1xyXG5cclxuICAgIGNvbnN0IGdyaWQgPSB0aGlzLm1fY29sR3JpZD8uZ3JpZDtcclxuICAgIGNvbnN0IHZpZXdUYWJsZSA9IGdyaWQ/LnZpZXc7XHJcblxyXG4gICAgaWYgKERHLnRvRGFydChncm9rLnNoZWxsLnYpICE9PSBERy50b0RhcnQodmlld1RhYmxlKSkge1xyXG4gICAgICByZXR1cm47XHJcbiAgICB9XHJcblxyXG4gICAgaWYoZS5kZWx0YVggIT09IDAgfHwgZS5kZWx0YVogIT09IDApIHtcclxuICAgICAgcmV0dXJuO1xyXG4gICAgfVxyXG5cclxuICAgIHNldFRpbWVvdXQoKCkgPT57XHJcbiAgICAgIGNvbnN0IGVlID0gbmV3IFdoZWVsRXZlbnQoZS50eXBlLCBlKTtcclxuICAgICAgdHJ5e2dyaWQub3ZlcmxheS5kaXNwYXRjaEV2ZW50KGVlKTt9XHJcbiAgICAgIGNhdGNoKGV4KSB7XHJcbiAgICAgICAgLy9jb25zb2xlLmVycm9yKGV4Lm1lc3NhZ2UpO1xyXG4gICAgICB9XHJcbiAgICB9LCAxKTtcclxuXHJcblxyXG4gICAgaWYodHJ1ZSlcclxuICAgICAgcmV0dXJuO1xyXG4gICAgLy9lLmNsaWVudFggPSA1O1xyXG5cclxuXHJcbiAgICBpZih0aGlzLm1fbldoZWVsQ291bnQgPT09IDEpIHtcclxuICAgICAgLy9zY3JvbGwgK1xyXG4gICAgICBjb25zdCBuUm93Q291bnQgPSBHcmlkVXRpbHMuZ2V0R3JpZFZpc2libGVSb3dDb3VudChncmlkKTtcclxuICAgICAgY29uc3Qgc2Nyb2xsWSA9IGdyaWQudmVydFNjcm9sbDtcclxuICAgICAgaWYoblJvd0NvdW50IC0xID4gc2Nyb2xsWS5tYXgpIHtcclxuICAgICAgICBzY3JvbGxZLnNldFZhbHVlcyhzY3JvbGxZLm1pblJhbmdlLCBzY3JvbGxZLm1heFJhbmdlLCBzY3JvbGxZLm1pbiArIDEsIHNjcm9sbFkubWF4ICsgMSk7XHJcbiAgICAgIH1cclxuICAgICAgdGhpcy5tX25XaGVlbENvdW50ID0gMDtcclxuICAgIH1cclxuICAgIGVsc2UgaWYodGhpcy5tX25XaGVlbENvdW50ID09PSAtMSlcclxuICAgIHtcclxuICAgICAgLy9zY3JvbGwgLVxyXG4gICAgICBjb25zdCBzY3JvbGxZID0gZ3JpZC52ZXJ0U2Nyb2xsO1xyXG4gICAgICBpZihzY3JvbGxZLm1pbiA+PTEpIHtcclxuICAgICAgICBzY3JvbGxZLnNldFZhbHVlcyhzY3JvbGxZLm1pblJhbmdlLCBzY3JvbGxZLm1heFJhbmdlLCBzY3JvbGxZLm1pbiAtIDEsIHNjcm9sbFkubWF4IC0gMSk7XHJcbiAgICAgIH1cclxuICAgICAgdGhpcy5tX25XaGVlbENvdW50ID0gMDtcclxuICAgIH1cclxuICAgIGVsc2Uge1xyXG4gICAgICB0aGlzLm1fbldoZWVsQ291bnQgPSBlLmRlbHRhWSA+IDAgPyAxIDogLTE7XHJcbiAgICB9XHJcbiAgfVxyXG5cclxuXHJcbiAgcHJpdmF0ZSBwYWludChnIDogQ2FudmFzUmVuZGVyaW5nQ29udGV4dDJEIHwgbnVsbCwgZ3JpZCA6IERHLkdyaWQpIDogdm9pZCB7XHJcbiAgICAvL2NvbnN0IG5XRGl2ID0gZW50cnkuY29udGVudEJveFNpemUgPyBlbnRyeS5jb250ZW50Qm94U2l6ZVswXS5pbmxpbmVTaXplIDogZW50cnkuY29udGVudFJlY3Qud2lkdGg7XHJcblxyXG4gICAgaWYoZyA9PT0gbnVsbCkge1xyXG4gICAgICByZXR1cm47XHJcbiAgICB9XHJcblxyXG4gICAgaWYodGhpcy5tX3Jvb3QgPT09IG51bGwpIHtcclxuICAgICAgdGhyb3cgbmV3IEVycm9yKCdSb290IGNhbm5vdCBiZSBudWxsLicpO1xyXG4gICAgfVxyXG5cclxuICAgIGlmKHRoaXMubV9jb2xHcmlkID09PSBudWxsKSB7XHJcbiAgICAgIHRocm93IG5ldyBFcnJvcignQ29sdW1uIGdyaWQgY2Fubm90IGJlIG51bGwuJyk7XHJcbiAgICB9XHJcbiAgICBjb25zdCBkZnJhbWUgPSBncmlkLmRhdGFGcmFtZTtcclxuICAgIGNvbnN0IG5XID0gdGhpcy5tX3Jvb3Qub2Zmc2V0V2lkdGg7XHJcbiAgICBjb25zdCBuSCA9IHRoaXMubV9yb290Lm9mZnNldEhlaWdodDtcclxuXHJcbiAgICBnLmZpbGxTdHlsZSA9IFwid2hpdGVcIjtcclxuICAgIGcuZmlsbFJlY3QoMCwwLCBuVyp3aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbywgbkgqd2luZG93LmRldmljZVBpeGVsUmF0aW8pO1xyXG5cclxuICAgIGlmKHRoaXMubV9jb2xHcmlkLm5hbWUgPT09IG51bGwpXHJcbiAgICAgIHJldHVybjtcclxuXHJcbiAgICBjb25zdCBiaXRzZXRGaWx0ZXIgPSBkZnJhbWUuZmlsdGVyO1xyXG4gICAgaWYoYml0c2V0RmlsdGVyLmZhbHNlQ291bnQgPT09IGRmcmFtZS5yb3dDb3VudClcclxuICAgICAgcmV0dXJuOyAvL2V2ZXJ5dGhpbmcgaXMgZmlsdGVyZWRcclxuXHJcbiAgICAvL2NvbHVtbiBIZWFkZXJcclxuICAgIGNvbnN0IG9wdGlvbnMgOiBhbnkgPSBncmlkLmdldE9wdGlvbnModHJ1ZSk7XHJcblxyXG4gICAgY29uc3QgZm9udENlbGxEZWZhdWx0ID0gb3B0aW9ucy5sb29rLmRlZmF1bHRDZWxsRm9udDtcclxuXHJcbiAgICBsZXQgZm9udCA9IG9wdGlvbnMubG9vay5jb2xIZWFkZXJGb250ID09IG51bGwgfHwgb3B0aW9ucy5sb29rLmNvbEhlYWRlckZvbnQgPT09IHVuZGVmaW5lZCA/IFwiYm9sZCAxNHB4IFZvbHRhIFRleHQsIEFyaWFsXCIgOiBvcHRpb25zLmxvb2suY29sSGVhZGVyRm9udDtcclxuICAgIGxldCBmb250U2NhbGVkID0gR3JpZFV0aWxzLnNjYWxlRm9udChmb250LCB3aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbyk7XHJcbiAgICBnLmZvbnQgPSBmb250U2NhbGVkO1xyXG5cclxuICAgIGxldCBzdHIgPSBUZXh0VXRpbHMudHJpbVRleHQodGhpcy5tX2NvbEdyaWQubmFtZSwgZywgblcpO1xyXG5cclxuICAgIGNvbnN0IHRtID0gZy5tZWFzdXJlVGV4dChzdHIpO1xyXG4gICAgY29uc3QgbldMYWJlbCA9IHRtLndpZHRoO1xyXG5cclxuICAgIGNvbnN0IG5Bc2NlbnQgPSBNYXRoLmFicyh0bS5hY3R1YWxCb3VuZGluZ0JveEFzY2VudCk7XHJcbiAgICBjb25zdCBuRGVzY2VudCA9IHRtLmFjdHVhbEJvdW5kaW5nQm94RGVzY2VudDtcclxuICAgIGNvbnN0IG5IRm9udCA9ICBuQXNjZW50ICsgbkRlc2NlbnQ7Ly8gKyAyKm5ZSW5zZXQ7XHJcblxyXG4gICAgLy9sZXQgY2VsbENIID0gZ3JpZC5jZWxsKHRoaXMubV9jb2xHcmlkLm5hbWUsIC0xKTtcclxuICAgIC8vbGV0IHJlbmRlcmVyID0gY2VsbENILnJlbmRlcmVyO1xyXG5cclxuICAgIGxldCBuWCA9IDA7XHJcbiAgICBsZXQgblkgPSAwO1xyXG4gICAgY29uc3QgbkhDSCA9IEdyaWRVdGlscy5nZXRHcmlkQ29sdW1uSGVhZGVySGVpZ2h0KGdyaWQpKndpbmRvdy5kZXZpY2VQaXhlbFJhdGlvO1xyXG4gICAgZy50ZXh0QWxpZ24gPSAnc3RhcnQnO1xyXG4gICAgZy5maWxsU3R5bGUgPSBcIkJsYWNrXCI7XHJcbiAgICBsZXQgbllPZmZzZXQgPSBNYXRoLmZsb29yKChuSENIIC0gbkhGb250KS8yKTtcclxuICAgIGNvbnN0IG5YWCA9IG5YICsgKChuVyp3aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbyAtIG5XTGFiZWwpID4+IDEpO1xyXG4gICAgbGV0IG5ZWSA9IChuWSArIG5IQ0ggLSBNYXRoLmNlaWwoMyp3aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbykpOy8vLTIqd2luZG93LmRldmljZVBpeGVsUmF0aW8pO1xyXG4gICAgLy9vbnNvbGUubG9nKFwiblhYIFwiICsgblhYICsgXCIgbllZID0gXCIgKyBuWVkgKyBcIiBDSEggXCIgKyBuSENIKTtcclxuICAgIGcuZmlsbFRleHQoc3RyLCBuWFgsIG5ZWSk7XHJcblxyXG4gICAgLy9pZihvcHRpb25zLmxvb2suc2hvd1Jvd0dyaWRsaW5lcykge1xyXG5cclxuXHJcbiAgICAvL31cclxuXHJcblxyXG5cclxuICAgIC8vUmVndWxhciBjZWxsc1xyXG4gICAgY29uc3QgblJvd0N1cnJlbnQgPSAgZGZyYW1lLmN1cnJlbnRSb3cuaWR4O1xyXG4gICAgY29uc3QgYml0c2V0U2VsID0gZGZyYW1lLnNlbGVjdGlvbjtcclxuXHJcbiAgICBjb25zdCBhclJvd3NNaW5NYXggPSBbLTEsLTFdO1xyXG4gICAgR3JpZFV0aWxzLmZpbGxWaXNpYmxlVmlld3BvcnRSb3dzKGFyUm93c01pbk1heCwgZ3JpZCk7XHJcbiAgICBjb25zdCBuUm93TWluID0gYXJSb3dzTWluTWF4WzBdO1xyXG4gICAgY29uc3QgblJvd01heCA9IGFyUm93c01pbk1heFsxXTtcclxuXHJcbiAgICAvL2NvbnNvbGUubG9nKG5Sb3dNaW4gKyBcIiBcIiArIG5Sb3dNYXgpO1xyXG4gICAgY29uc3QgbkhSb3cgPSBHcmlkVXRpbHMuZ2V0R3JpZFJvd0hlaWdodChncmlkKTtcclxuICAgIG5ZT2Zmc2V0ID0gbkhDSDtcclxuICAgIGNvbnN0IG5IUm93R3JpZCA9IG5IUm93KndpbmRvdy5kZXZpY2VQaXhlbFJhdGlvO1xyXG4gICAgbGV0IGNlbGxSSCA9IG51bGw7XHJcblxyXG4gICAgbGV0IG5XVyA9IG5XKndpbmRvdy5kZXZpY2VQaXhlbFJhdGlvO1xyXG4gICAgLy9jb25zdCBuSEggPSBuSFJvd0dyaWQ7XHJcblxyXG4gICAgY29uc3QgYXJUYWJsZVJvd3MgPSBuZXcgQXJyYXkoblJvd01heCAtIG5Sb3dNaW4gKzEpO1xyXG4gICAgbGV0IG5Sb3dUYWJsZSA9IC0xO1xyXG4gICAgbGV0IGJTZWwgPSBmYWxzZTtcclxuICAgIGZvcihsZXQgblJHPW5Sb3dNaW47IG5SRzw9blJvd01heDsgKytuUkcpIHtcclxuICAgICAgdHJ5IHtcclxuICAgICAgICBjZWxsUkggPSBncmlkLmNlbGwodGhpcy5tX2NvbEdyaWQubmFtZSwgblJHKTtcclxuICAgICAgfSBjYXRjaCAoZSkgLy90byBhZGRyZXNzIERHIGJ1ZyB3aGVuIGV2ZXJ5dGhpbmcgaXMgZmlsdGVyZWRcclxuICAgICAge1xyXG4gICAgICAgIGNvbnRpbnVlO1xyXG4gICAgICB9XHJcblxyXG4gICAgICBpZiAoY2VsbFJILnRhYmxlUm93SW5kZXggPT09IHVuZGVmaW5lZCkvL0RHIGJ1Z1xyXG4gICAgICAgIGNvbnRpbnVlO1xyXG5cclxuICAgICAgblJvd1RhYmxlID0gY2VsbFJILnRhYmxlUm93SW5kZXggPT09IG51bGwgPyAtMSA6IGNlbGxSSC50YWJsZVJvd0luZGV4O1xyXG4gICAgICBhclRhYmxlUm93c1tuUkcgLSBuUm93TWluXSA9IG5Sb3dUYWJsZTtcclxuXHJcbiAgICAgIG5ZWSA9IG5ZT2Zmc2V0ICsgKG5SRyAtIG5Sb3dNaW4pICogbkhSb3dHcmlkO1xyXG5cclxuICAgICAgbGV0IHJlbmRlcmVyOiBhbnkgPSBHcmlkVXRpbHMuZ2V0R3JpZENvbHVtblJlbmRlcmVyKGNlbGxSSC5ncmlkQ29sdW1uKTtcclxuICAgICAgaWYgKHJlbmRlcmVyID09PSBudWxsKSB7XHJcbiAgICAgICAgdHJ5IHtcclxuICAgICAgICAgIHJlbmRlcmVyID0gY2VsbFJILnJlbmRlcmVyO1xyXG4gICAgICAgIH0gY2F0Y2ggKGUpIHtcclxuICAgICAgICAgIGNvbnNvbGUuZXJyb3IoXCJDb3VsZCBub3Qgb2J0YWluIHJlbmRlcmVyIGZvciBERyBjZWxsLiBERyBidWcgXCIgKyB0aGlzLm1fY29sR3JpZC5uYW1lICsgXCIgcm93IFwiICsgblJHKTtcclxuICAgICAgICAgIGNvbnRpbnVlO1xyXG4gICAgICAgIH1cclxuICAgICAgfVxyXG5cclxuICAgICAgaWYgKHJlbmRlcmVyID09PSBudWxsIHx8IHJlbmRlcmVyID09PSB1bmRlZmluZWQpIHtcclxuICAgICAgICBjb25zb2xlLmVycm9yKFwiQ291bGRuJ3QgZmluZCByZW5kZXJlciBmb3IgcGlubmVkIGNvbHVtbiBcIiArIHRoaXMubV9jb2xHcmlkLm5hbWUgKyBcIiByb3cgXCIgKyBuUkcpO1xyXG4gICAgICAgIGNvbnRpbnVlO1xyXG4gICAgICB9XHJcblxyXG4gICAgICAvL2xldCBuWVkgPSBuWTsvLyp3aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbztcclxuXHJcblxyXG4gICAgICBmb250ID0gY2VsbFJILnN0eWxlLmZvbnQ7XHJcbiAgICAgIGZvbnRTY2FsZWQgPSBHcmlkVXRpbHMuc2NhbGVGb250KGZvbnQsIHdpbmRvdy5kZXZpY2VQaXhlbFJhdGlvKTtcclxuICAgICAgaWYgKGZvbnRTY2FsZWQgIT09IG51bGwpIHtcclxuICAgICAgICBjZWxsUkguc3R5bGUuZm9udCA9IGZvbnRTY2FsZWQ7XHJcbiAgICAgIH1cclxuXHJcbiAgICAgIGlmIChuVyA+IDAgJiYgbkhSb3dHcmlkID4gMCkgeyAvL3RvIGFkZHJlc3MgYSBidWcgY2F1c2VkIGVpdGhlciBERyBvciBjbGllbnQgYXBwXHJcbiAgICAgICAgdHJ5IHtcclxuICAgICAgICAgIGlmIChyZW5kZXJlci5uYW1lID09PSAnTW9sZWN1bGUnKSB7XHJcbiAgICAgICAgICAgIHJlbmRlcmVyLnJlbmRlcihnLCAwLCBuWVkvd2luZG93LmRldmljZVBpeGVsUmF0aW8sIG5XVy93aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbywgbkhSb3dHcmlkL3dpbmRvdy5kZXZpY2VQaXhlbFJhdGlvLCBjZWxsUkgsIGNlbGxSSC5zdHlsZSk7XHJcbiAgICAgICAgICB9XHJcbiAgICAgICAgICBlbHNlIHJlbmRlcmVyLnJlbmRlcihnLCAwLCBuWVksIG5XVywgbkhSb3dHcmlkLCBjZWxsUkgsIGNlbGxSSC5zdHlsZSk7XHJcblxyXG4gICAgICAgIH0gY2F0Y2ggKGUpIHtcclxuICAgICAgICAgIGNvbnNvbGUuZXJyb3IoXCJDb3VsZCBub3QgcGFpbnQgY2VsbCBmb3IgcGlubmVkIGNvbHVtbiBcIiArIHRoaXMubV9jb2xHcmlkLm5hbWUgKyBcIiByb3cgXCIgKyBuUkcpO1xyXG4gICAgICAgICAgY29udGludWU7XHJcbiAgICAgICAgICAvL3Rocm93IGU7XHJcbiAgICAgICAgfVxyXG4gICAgICB9XHJcbiAgICB9XHJcblxyXG5cclxuICAgIC8vUGFpbnQgR3JpZFxyXG4gICAgZy5zdHJva2VTdHlsZSA9IFwiR2FpbnNib3JvXCI7XHJcbiAgICBnLmJlZ2luUGF0aCgpO1xyXG4gICAgZy5tb3ZlVG8oMCwgblkqd2luZG93LmRldmljZVBpeGVsUmF0aW8pO1xyXG4gICAgZy5saW5lVG8oMCwgKG5ZICsgbkhDSC0xKndpbmRvdy5kZXZpY2VQaXhlbFJhdGlvKSk7XHJcbiAgICBnLnN0cm9rZSgpO1xyXG5cclxuICAgIGcuYmVnaW5QYXRoKCk7XHJcbiAgICBnLm1vdmVUbygwLCBuWU9mZnNldCArIDEpO1xyXG4gICAgZy5saW5lVG8obldXLCBuWU9mZnNldCArIDEpO1xyXG4gICAgZy5zdHJva2UoKTtcclxuXHJcbiAgICBmb3IobGV0IG5SRz1uUm93TWluOyBuUkc8PW5Sb3dNYXg7ICsrblJHKVxyXG4gICAge1xyXG4gICAgICBuWVkgPSBuWU9mZnNldCArIChuUkcgLSBuUm93TWluKSAqIG5IUm93R3JpZDtcclxuXHJcbiAgICAgIC8vaWYob3B0aW9ucy5sb29rLnNob3dSb3dHcmlkbGluZXMpIHtcclxuXHJcbiAgICAgICAgZy5iZWdpblBhdGgoKTtcclxuICAgICAgICBnLm1vdmVUbygwLCBuWVkgKyBuSFJvd0dyaWQrMSk7XHJcbiAgICAgICAgZy5saW5lVG8obldXLCBuWVkgKyBuSFJvd0dyaWQrMSk7XHJcbiAgICAgICAgZy5zdHJva2UoKTtcclxuXHJcbiAgICAgICAgZy5iZWdpblBhdGgoKTtcclxuICAgICAgICBnLm1vdmVUbygwLCBuWVkpO1xyXG4gICAgICAgIGcubGluZVRvKDAsIG5ZWSArIG5IUm93R3JpZCsxKTtcclxuICAgICAgICBnLnN0cm9rZSgpO1xyXG4gICAgICAvL31cclxuICAgICAgblJvd1RhYmxlID0gYXJUYWJsZVJvd3NbblJHIC0gblJvd01pbl07XHJcbiAgICAgIHRyeXtiU2VsID0gblJvd1RhYmxlID09PSB1bmRlZmluZWQgfHwgblJvd1RhYmxlIDwgMCA/IGZhbHNlIDogYml0c2V0U2VsLmdldChuUm93VGFibGUpO31cclxuICAgICAgY2F0Y2ggKGUpe1xyXG4gICAgICAgIGNvbnNvbGUuZXJyb3IoJ1BhaW50RXJyb3I6IHJvd19taW46ICcgKyBuUm93TWluICsgJyByb3dfbWF4OiAnICsgblJvd01heCArICcgblIgJyArIG5SRyArICcgJyArIG5Sb3dUYWJsZSk7XHJcbiAgICAgICAgdGhyb3cgZTtcclxuICAgICAgfVxyXG4gICAgICBpZihiU2VsKVxyXG4gICAgICB7XHJcbiAgICAgICAgZy5nbG9iYWxBbHBoYSA9IDAuMjtcclxuICAgICAgICBnLmZpbGxTdHlsZSA9IFBpbm5lZENvbHVtbi5TRUxFQ1RJT05fQ09MT1I7XHJcbiAgICAgICAgZy5maWxsUmVjdCgwLCBuWVksIG5XVywgbkhSb3dHcmlkKTtcclxuICAgICAgICBnLmdsb2JhbEFscGhhID0gMTtcclxuICAgICAgfVxyXG5cclxuICAgICAgaWYoblJvd0N1cnJlbnQgPT09IG5Sb3dUYWJsZSlcclxuICAgICAge1xyXG4gICAgICAgIGcuZ2xvYmFsQWxwaGEgPSAwLjI7XHJcbiAgICAgICAgZy5maWxsU3R5bGUgPSBQaW5uZWRDb2x1bW4uQUNUSVZFX0NFTExfQ09MT1I7XHJcbiAgICAgICAgZy5maWxsUmVjdCgwLCBuWVksIG5XVywgbkhSb3dHcmlkKTtcclxuICAgICAgICBnLmdsb2JhbEFscGhhID0gMTtcclxuICAgICAgfVxyXG4gICAgfVxyXG4gIH1cclxuXHJcblxyXG4gIHByaXZhdGUgc3RhdGljIGhpdFRlc3RSb3dzKGVDYW52YXNQaW5uZWQgOiBIVE1MQ2FudmFzRWxlbWVudCwgZ3JpZCA6IERHLkdyaWQsIGUgOiBNb3VzZUV2ZW50LCBiQm9yZGVyIDogYm9vbGVhbiwgYXJYWU9uQ2VsbCA6IEFycmF5PG51bWJlcj4gfCB1bmRlZmluZWQpXHJcbiAge1xyXG4gICAgY29uc3QgcmVjdCA9IGVDYW52YXNQaW5uZWQuZ2V0Qm91bmRpbmdDbGllbnRSZWN0KCk7XHJcbiAgICBjb25zdCBzY3JvbGxMZWZ0PSB3aW5kb3cucGFnZVhPZmZzZXQgfHwgZG9jdW1lbnQuZG9jdW1lbnRFbGVtZW50LnNjcm9sbExlZnQ7XHJcbiAgICBjb25zdCBzY3JvbGxUb3AgPSB3aW5kb3cucGFnZVlPZmZzZXQgfHwgZG9jdW1lbnQuZG9jdW1lbnRFbGVtZW50LnNjcm9sbFRvcDtcclxuICAgIGNvbnN0IG5ZID0gcmVjdC50b3AgICsgc2Nyb2xsVG9wO1xyXG4gICAgY29uc3QgblggPSByZWN0LmxlZnQgKyBzY3JvbGxMZWZ0O1xyXG5cclxuICAgIGlmKG5YIDw9IGUuY2xpZW50WCAmJiBlLmNsaWVudFggPD0gblggKyBlQ2FudmFzUGlubmVkLm9mZnNldFdpZHRoKSAgIC8vb24gdGhlIHJvd3MgaGVhZGVyXHJcbiAgICB7XHJcbiAgICAgIGNvbnN0IG5ISGVhZGVyQ29scyA9IEdyaWRVdGlscy5nZXRHcmlkQ29sdW1uSGVhZGVySGVpZ2h0KGdyaWQpO1xyXG4gICAgICBjb25zdCBuSFJvd0dyaWQgPSBHcmlkVXRpbHMuZ2V0R3JpZFJvd0hlaWdodChncmlkKTtcclxuXHJcbiAgICAgIGNvbnN0IGFyTWluTWF4Um93cyA9IFstMSwtMV07XHJcbiAgICAgIEdyaWRVdGlscy5maWxsVmlzaWJsZVZpZXdwb3J0Um93cyhhck1pbk1heFJvd3MsIGdyaWQpO1xyXG4gICAgICBjb25zdCBuUm93TWluID0gYXJNaW5NYXhSb3dzWzBdO1xyXG4gICAgICBjb25zdCBuUm93TWF4ID0gYXJNaW5NYXhSb3dzWzFdO1xyXG5cclxuICAgICAgY29uc3QgbllNb3VzZU9uSGVhZGVyID0gZS5jbGllbnRZIC0gblk7XHJcblxyXG4gICAgICBsZXQgbllCb3JkZXIgPSAtMTtcclxuICAgICAgbGV0IG5ZRGlmZiA9IC0xO1xyXG5cclxuICAgICAgZm9yKGxldCBuUm93PW5Sb3dNaW47IG5Sb3c8PSBuUm93TWF4OyArK25Sb3cpXHJcbiAgICAgIHtcclxuICAgICAgICBuWUJvcmRlciA9IG5ISGVhZGVyQ29scyArIChuUm93IC0gblJvd01pbisxKSpuSFJvd0dyaWQ7XHJcbiAgICAgICAgbllEaWZmID0gbllNb3VzZU9uSGVhZGVyIC0gbllCb3JkZXI7XHJcblxyXG4gICAgICAgIGlmKGJCb3JkZXIgJiYgTWF0aC5hYnMobllEaWZmKSA8PSBQaW5uZWRDb2x1bW4uWV9SRVNJWkVfU0VOU0lUSVZJVFkpXHJcbiAgICAgICAge1xyXG4gICAgICAgICAgcmV0dXJuIG5Sb3c7XHJcbiAgICAgICAgfVxyXG5cclxuICAgICAgICBpZighYkJvcmRlciAmJiBuWUJvcmRlciAtIG5IUm93R3JpZCA8PSBuWU1vdXNlT25IZWFkZXIgJiYgbllNb3VzZU9uSGVhZGVyIDw9IG5ZQm9yZGVyKSB7XHJcblxyXG4gICAgICAgICAgaWYoYXJYWU9uQ2VsbCAhPT0gdW5kZWZpbmVkKSB7XHJcbiAgICAgICAgICAgIGFyWFlPbkNlbGxbMF0gPSBlLmNsaWVudFggLSBuWDtcclxuICAgICAgICAgICAgYXJYWU9uQ2VsbFsxXSA9IG5ZTW91c2VPbkhlYWRlciAtIG5ZQm9yZGVyICsgbkhSb3dHcmlkO1xyXG4gICAgICAgICAgfVxyXG5cclxuICAgICAgICAgIHJldHVybiBuUm93O1xyXG4gICAgICAgIH1cclxuICAgICAgfVxyXG4gICAgfVxyXG5cclxuICAgIHJldHVybiAtMTtcclxuICB9XHJcbn1cclxuIl19