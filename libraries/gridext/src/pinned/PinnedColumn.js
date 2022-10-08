import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as GridUtils from '../utils/GridUtils';
import * as TextUtils from '../utils/TextUtils';
import { ColorUtils } from '../utils/ColorUtils';
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
        this.m_handlerColsRemoved.unsubscribe();
        this.m_handlerColsRemoved = null;
        this.m_handlerColNameChanged.unsubscribe();
        this.m_handlerColNameChanged = null;
        this.m_handlerVScroll.unsubscribe();
        this.m_handlerVScroll = null;
        this.m_handlerRowsResized.unsubscribe();
        this.m_handlerRowsResized = null;
        this.m_handlerRowsSorted.unsubscribe();
        this.m_handlerRowsSorted = null;
        this.m_handlerRowsFiltering.unsubscribe();
        this.m_handlerRowsFiltering = null;
        this.m_handlerCurrRow.unsubscribe();
        this.m_handlerCurrRow = null;
        this.m_handlerSel.unsubscribe();
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
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJmaWxlIjoiUGlubmVkQ29sdW1uLmpzIiwic291cmNlUm9vdCI6IiIsInNvdXJjZXMiOlsiUGlubmVkQ29sdW1uLnRzIl0sIm5hbWVzIjpbXSwibWFwcGluZ3MiOiJBQUFBLE9BQU8sS0FBSyxJQUFJLE1BQU0sbUJBQW1CLENBQUM7QUFDMUMsT0FBTyxLQUFLLEVBQUUsTUFBTSxpQkFBaUIsQ0FBQztBQUN0QyxPQUFPLEtBQUssRUFBRSxNQUFNLGlCQUFpQixDQUFDO0FBQ3RDLE9BQU8sS0FBSyxTQUFTLE1BQU0sb0JBQW9CLENBQUM7QUFDaEQsT0FBTyxLQUFLLFNBQVMsTUFBTSxvQkFBb0IsQ0FBQztBQUNoRCxPQUFPLEVBQUMsVUFBVSxFQUFDLE1BQU0scUJBQXFCLENBQUM7QUFFL0MsT0FBTyxFQUFFLGtCQUFrQixFQUFDLE1BQU0sZ0NBQWdDLENBQUM7QUFDbkUsT0FBTyxLQUFLLFdBQVcsTUFBTSxlQUFlLENBQUM7QUFFN0MsT0FBTyxFQUFDLGVBQWUsRUFBQyxNQUFNLHVCQUF1QixDQUFDO0FBQ3RELE9BQU8sRUFBYyxNQUFNLEVBQUMsTUFBTSxpQkFBaUIsQ0FBQztBQUNwRCw0Q0FBNEM7QUFHNUM7Ozs7Ozs7Ozs7Ozs7Ozs7RUFnQkU7QUFFRixTQUFTLFdBQVcsQ0FBQyxJQUFrQjtJQUNyQyxNQUFNLE9BQU8sR0FBRyxJQUFJLENBQUMsVUFBVSxDQUFDO0lBQ2hDLElBQUksT0FBTyxLQUFLLElBQUksSUFBSSxPQUFPLEtBQUssU0FBUyxFQUFFO1FBQzdDLE1BQU0sSUFBSSxLQUFLLENBQUMsNENBQTRDLENBQUMsQ0FBQztLQUMvRDtJQUVELElBQUksUUFBUSxHQUFHLFNBQVMsQ0FBQyxxQkFBcUIsQ0FBQyxPQUFPLENBQUMsQ0FBQztJQUN4RCxJQUFHLFFBQVEsWUFBWSxrQkFBa0IsRUFBRTtRQUN6QyxPQUFPLFFBQVEsQ0FBQztLQUNqQjtJQUVELE9BQU8sSUFBSSxDQUFDLFFBQVEsQ0FBQztBQUN2QixDQUFDO0FBR0QsU0FBUyxPQUFPLENBQUMsT0FBdUI7SUFD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QUFDLFVBQVUsQ0FBQyxVQUFVLENBQUMsQ0FBQyxDQUFDLDZCQUE2QjtBQUMxRixpQ0FBb0IsR0FBRyxDQUFDLENBQUMiLCJzb3VyY2VzQ29udGVudCI6WyJpbXBvcnQgKiBhcyBncm9rIGZyb20gJ2RhdGFncm9rLWFwaS9ncm9rJztcclxuaW1wb3J0ICogYXMgREcgZnJvbSAnZGF0YWdyb2stYXBpL2RnJztcclxuaW1wb3J0ICogYXMgdWkgZnJvbSAnZGF0YWdyb2stYXBpL3VpJztcclxuaW1wb3J0ICogYXMgR3JpZFV0aWxzIGZyb20gJy4uL3V0aWxzL0dyaWRVdGlscyc7XHJcbmltcG9ydCAqIGFzIFRleHRVdGlscyBmcm9tICcuLi91dGlscy9UZXh0VXRpbHMnO1xyXG5pbXBvcnQge0NvbG9yVXRpbHN9IGZyb20gJy4uL3V0aWxzL0NvbG9yVXRpbHMnO1xyXG5pbXBvcnQgKiBhcyByeGpzIGZyb20gJ3J4anMnO1xyXG5pbXBvcnQgeyBHcmlkQ2VsbFJlbmRlcmVyRXh9IGZyb20gXCIuLi9yZW5kZXJlci9HcmlkQ2VsbFJlbmRlcmVyRXhcIjtcclxuaW1wb3J0ICogYXMgUGlubmVkVXRpbHMgZnJvbSBcIi4vUGlubmVkVXRpbHNcIjtcclxuaW1wb3J0IHtnZXRHcmlkRGFydFBvcHVwTWVudSwgaXNIaXRUZXN0T25FbGVtZW50fSBmcm9tIFwiLi4vdXRpbHMvR3JpZFV0aWxzXCI7XHJcbmltcG9ydCB7TW91c2VEaXNwYXRjaGVyfSBmcm9tIFwiLi4vdWkvTW91c2VEaXNwYXRjaGVyXCI7XHJcbmltcG9ydCB7Q29sdW1uc0FyZ3MsIHRvRGFydH0gZnJvbSBcImRhdGFncm9rLWFwaS9kZ1wiO1xyXG4vL2ltcG9ydCB7VGFibGVWaWV3fSBmcm9tIFwiZGF0YWdyb2stYXBpL2RnXCI7XHJcblxyXG5cclxuLypcclxuY29uc3QgaFN1YnNjcmliZXIgID0gZ3Jvay5ldmVudHMub25WaWV3TGF5b3V0QXBwbGllZC5zdWJzY3JpYmUoKGxheW91dCA6IERHLlZpZXdMYXlvdXQpID0+IHtcclxuICBjb25zdCB2aWV3IDogREcuVGFibGVWaWV3ID0gbGF5b3V0LnZpZXcgYXMgVGFibGVWaWV3O1xyXG4gIGNvbnN0IGl0Vmlld2VycyA9IHZpZXcudmlld2VycztcclxuICBjb25zdCBhclZpZXdlcnMgPSBBcnJheS5mcm9tKGl0Vmlld2Vycyk7XHJcblxyXG4gIGxldCB2aWV3ZXIgPSBudWxsO1xyXG4gIGNvbnN0IG5WaWV3ZXJDb3VudCA9IGFyVmlld2Vycy5sZW5ndGg7XHJcbiAgZm9yIChsZXQgbiA9IDA7IG4gPCBuVmlld2VyQ291bnQ7ICsrbikge1xyXG4gICAgdmlld2VyID0gYXJWaWV3ZXJzW25dO1xyXG4gICAgaWYgKHZpZXdlci50eXBlICE9PSBcIkdyaWRcIilcclxuICAgICAgY29udGludWU7XHJcblxyXG4gICAgUGlubmVkVXRpbHMuaW5zdGFsbFBpbm5lZENvbHVtbnModmlld2VyIGFzIERHLkdyaWQpO1xyXG4gIH1cclxufSk7XHJcbiovXHJcblxyXG5mdW5jdGlvbiBnZXRSZW5kZXJlcihjZWxsIDogREcuR3JpZENlbGwpIDogR3JpZENlbGxSZW5kZXJlckV4IHwgREcuR3JpZENlbGxSZW5kZXJlciB7XHJcbiAgY29uc3QgY29sR3JpZCA9IGNlbGwuZ3JpZENvbHVtbjtcclxuICBpZiAoY29sR3JpZCA9PT0gbnVsbCB8fCBjb2xHcmlkID09PSB1bmRlZmluZWQpIHtcclxuICAgIHRocm93IG5ldyBFcnJvcignR3JpZCBjZWxsIGlzIGRldGFjaGVkIGZyb20gdGhlIEdyaWQgY29sdW1uJyk7XHJcbiAgfVxyXG5cclxuICBsZXQgcmVuZGVyZXIgPSBHcmlkVXRpbHMuZ2V0R3JpZENvbHVtblJlbmRlcmVyKGNvbEdyaWQpO1xyXG4gIGlmKHJlbmRlcmVyIGluc3RhbmNlb2YgR3JpZENlbGxSZW5kZXJlckV4KSB7XHJcbiAgICByZXR1cm4gcmVuZGVyZXI7XHJcbiAgfVxyXG5cclxuICByZXR1cm4gY2VsbC5yZW5kZXJlcjtcclxufVxyXG5cclxuXHJcbmZ1bmN0aW9uIGdldEdyaWQoY29sR3JpZCA6IERHLkdyaWRDb2x1bW4pIDogREcuR3JpZCB8IG51bGwge1xyXG4gIGxldCBncmlkIDogREcuR3JpZCB8IG51bGwgPSBjb2xHcmlkLmdyaWQ7XHJcbiAgaWYoIGdyaWQgPT09IG51bGwpIHtcclxuICAgIGdyaWQgPSBHcmlkVXRpbHMuZ2V0SW5zdGFsbGVkR3JpZEZvckNvbHVtbihjb2xHcmlkKTtcclxuICAgIGlmKGdyaWQgaW5zdGFuY2VvZiBERy5HcmlkKVxyXG4gICAgICByZXR1cm4gZ3JpZDtcclxuICB9XHJcblxyXG4gIHJldHVybiBncmlkO1xyXG59XHJcblxyXG5cclxuZnVuY3Rpb24gbm90aWZ5QWxsQ29sc1Jvd3NSZXNpemVkKGdyaWQgOiBERy5HcmlkLCBuSFJvd3MgOiBudW1iZXIsIGJBZGp1c3RpbmcgOiBib29sZWFuKSA6IHZvaWQge1xyXG5cclxuICBsZXQgcmVuZGVyZXIgOiBHcmlkQ2VsbFJlbmRlcmVyRXggfCBudWxsID0gbnVsbFxyXG4gIGxldCBjb2xHcmlkID0gbnVsbDtcclxuICBjb25zdCBsc3RDb2xzR3JpZCA9IGdyaWQuY29sdW1ucztcclxuICBjb25zdCBuQ29sQ291bnQgPSBsc3RDb2xzR3JpZC5sZW5ndGg7XHJcbiAgZm9yKGxldCBuQ29sPTA7IG5Db2w8bkNvbENvdW50OyArK25Db2wpIHtcclxuICAgIGNvbEdyaWQgPSBsc3RDb2xzR3JpZC5ieUluZGV4KG5Db2wpO1xyXG4gICAgaWYoY29sR3JpZCA9PT0gbnVsbCB8fCAhY29sR3JpZC52aXNpYmxlKXtcclxuICAgICAgY29udGludWVcclxuICAgIH1cclxuXHJcbiAgICByZW5kZXJlciA9IEdyaWRVdGlscy5nZXRHcmlkQ29sdW1uUmVuZGVyZXIoY29sR3JpZCk7XHJcbiAgICBpZiAocmVuZGVyZXIgaW5zdGFuY2VvZiBHcmlkQ2VsbFJlbmRlcmVyRXgpIHtcclxuICAgICAgcmVuZGVyZXIub25SZXNpemVIZWlnaHQoY29sR3JpZCwgZ3JpZCwgbkhSb3dzLCBiQWRqdXN0aW5nKTtcclxuICAgIH1cclxuICB9XHJcbn1cclxuXHJcblxyXG5mdW5jdGlvbiBub3RpZnlBbGxQaW5uZWRDb2xzUm93c1Jlc2l6ZWQoY29sUGlubmVkU291cmNlIDogUGlubmVkQ29sdW1uLCBuSFJvd3MgOiBudW1iZXIsIGJBZGp1c3RpbmcgOiBib29sZWFuKSA6IHZvaWQge1xyXG5cclxuICBjb25zdCBjb2xHcmlkU291cmNlICA9IGNvbFBpbm5lZFNvdXJjZS5nZXRHcmlkQ29sdW1uKCk7XHJcbiAgaWYoY29sR3JpZFNvdXJjZSA9PT0gbnVsbCl7XHJcbiAgICByZXR1cm47XHJcbiAgfVxyXG5cclxuICBjb25zdCBncmlkID0gZ2V0R3JpZChjb2xHcmlkU291cmNlKTtcclxuICBjb25zdCBkYXJ0ID0gREcudG9EYXJ0KGdyaWQpO1xyXG4gIGlmKGRhcnQubV9hclBpbm5lZENvbHMgPT09IHVuZGVmaW5lZCkge1xyXG4gICAgdGhyb3cgbmV3IEVycm9yKCdQaW5uZWQgQ29sdW1ucyBhcmUgbm90IGluc3RhbGxlZC4nKTtcclxuICB9XHJcblxyXG4gIGxldCByZW5kZXJlciA6IEdyaWRDZWxsUmVuZGVyZXJFeCB8IG51bGwgPSBudWxsXHJcbiAgbGV0IGNvbFBpbm5lZCA9IG51bGw7XHJcbiAgbGV0IGNvbEdyaWQgPSBudWxsO1xyXG4gIGNvbnN0IG5QaW5uZWRDb2xDb3VudCA9IGRhcnQubV9hclBpbm5lZENvbHMubGVuZ3RoO1xyXG4gIGZvcihsZXQgbkNvbFBpbj0wOyBuQ29sUGluPG5QaW5uZWRDb2xDb3VudDsgKytuQ29sUGluKSB7XHJcbiAgICBjb2xQaW5uZWQgPSBkYXJ0Lm1fYXJQaW5uZWRDb2xzW25Db2xQaW5dO1xyXG4gICAgY29sR3JpZCA9IGNvbFBpbm5lZC5tX2NvbEdyaWQ7XHJcbiAgICBpZihjb2xHcmlkID09PSBudWxsKSB7XHJcbiAgICAgIHRocm93IG5ldyBFcnJvcignUGlubmVkIENvbHVtbiBpcyBkZXRhY2hlZC4nKTtcclxuICAgIH1cclxuXHJcbiAgICByZW5kZXJlciA9IEdyaWRVdGlscy5nZXRHcmlkQ29sdW1uUmVuZGVyZXIoY29sR3JpZCk7XHJcbiAgICBpZiAocmVuZGVyZXIgaW5zdGFuY2VvZiBHcmlkQ2VsbFJlbmRlcmVyRXggICYmIGNvbFBpbm5lZC5tX3Jvb3QgIT09IG51bGwgJiYgZ3JpZCAhPT0gbnVsbCkge1xyXG4gICAgICByZW5kZXJlci5vblJlc2l6ZUhlaWdodChjb2xQaW5uZWQsIGdyaWQsIG5IUm93cywgYkFkanVzdGluZyk7XHJcbiAgICB9XHJcbiAgfVxyXG59XHJcblxyXG5cclxuY29uc3QgREVCVUcgOiBib29sZWFuID0gZmFsc2U7XHJcblxyXG5cclxuZXhwb3J0IGNsYXNzIFBpbm5lZENvbHVtbiB7XHJcblxyXG4gIHByaXZhdGUgc3RhdGljIE1JTl9ST1dfSEVJR0hUID0gMjA7XHJcbiAgcHJpdmF0ZSBzdGF0aWMgTUFYX1JPV19IRUlHSFQgPSA1MDA7XHJcbiAgcHJpdmF0ZSBzdGF0aWMgU0VMRUNUSU9OX0NPTE9SID0gQ29sb3JVdGlscy50b1JnYihDb2xvclV0aWxzLmNvbFNlbGVjdGlvbik7IC8vXCJyZ2JhKDIzNywgMjIwLCA4OCwgMC4xNSlcIjtcclxuICBwcml2YXRlIHN0YXRpYyBBQ1RJVkVfQ0VMTF9DT0xPUiA9IENvbG9yVXRpbHMudG9SZ2IoQ29sb3JVdGlscy5jdXJyZW50Um93KTsgLy9cInJnYmEoMTUzLCAyMzcsIDgyLCAwLjI1KVwiO1xyXG4gIHByaXZhdGUgc3RhdGljIFlfUkVTSVpFX1NFTlNJVElWSVRZID0gMjtcclxuXHJcbiAgcHJpdmF0ZSBtX2ZEZXZpY2VQaXhlbFJhdGlvIDogbnVtYmVyO1xyXG4gIHByaXZhdGUgbV9jb2xHcmlkIDogREcuR3JpZENvbHVtbiB8IG51bGw7XHJcbiAgcHJpdmF0ZSBtX3Jvb3QgOiBIVE1MQ2FudmFzRWxlbWVudCB8IG51bGw7XHJcbiAgcHJpdmF0ZSBtX25XaWR0aEJ1ZyA6IG51bWJlcjtcclxuICAvL3ByaXZhdGUgbV9vYnNlcnZlclJlc2l6ZSA6IFJlc2l6ZU9ic2VydmVyIHwgbnVsbDtcclxuICBwcml2YXRlIG1fb2JzZXJ2ZXJSZXNpemVHcmlkIDogUmVzaXplT2JzZXJ2ZXIgfCBudWxsO1xyXG4gIHByaXZhdGUgbV9oYW5kbGVyQ29sc1JlbW92ZWQgOiBhbnk7XHJcbiAgcHJpdmF0ZSBtX2hhbmRsZXJDb2xOYW1lQ2hhbmdlZCA6IGFueTtcclxuICBwcml2YXRlIG1faGFuZGxlclZTY3JvbGwgOiBhbnk7XHJcbiAgcHJpdmF0ZSBtX2hhbmRsZXJSb3dzRmlsdGVyaW5nIDogYW55O1xyXG4gIHByaXZhdGUgbV9oYW5kbGVyQ3VyclJvdyA6IGFueTtcclxuICBwcml2YXRlIG1faGFuZGxlclNlbCA6IGFueTtcclxuICAvL3ByaXZhdGUgbV9oYW5kbGVyRmlsdGVyIDogYW55O1xyXG4gIHByaXZhdGUgbV9oYW5kbGVyUm93c1Jlc2l6ZWQgOiBhbnk7XHJcbiAgcHJpdmF0ZSBtX2hhbmRsZXJSb3dzU29ydGVkIDogYW55O1xyXG5cclxuICBwcml2YXRlIG1fbkhSZXNpemVSb3dzQmVmb3JlRHJhZyA9IC0xO1xyXG4gIHByaXZhdGUgbV9uUmVzaXplUm93R3JpZERyYWdnaW5nID0gLTE7XHJcbiAgcHJpdmF0ZSBtX25ZUmVzaXplRHJhZ2dpbmdBbmNob3IgPSAtMTtcclxuICBwcml2YXRlIG1fblJlc2l6ZVJvd0dyaWRNb3ZpbmcgPSAtMTtcclxuXHJcbiAgcHJpdmF0ZSBtX25ZRHJhZ2dpbmdBbmNob3IgPSAtMTtcclxuICBwcml2YXRlIG1fblJvd0dyaWREcmFnZ2luZyA9IC0xO1xyXG5cclxuICBwcml2YXRlIG1fbldoZWVsQ291bnQgOiBudW1iZXIgPSAwO1xyXG5cclxuXHJcbiAgcHJpdmF0ZSBtX2FyWFlNb3VzZU9uQ2VsbERvd24gPSBbLTIsIC0yXTtcclxuICBwcml2YXRlIG1fYXJYWU1vdXNlT25DZWxsVXAgPSBbLTEsIC0xXTtcclxuICBwcml2YXRlIG1fYlNvcnRlZEFzY2VuZGluZyA6IGJvb2xlYW4gfCBudWxsID0gbnVsbDtcclxuXHJcbiAgcHJpdmF0ZSBtX2NlbGxDdXJyZW50IDogREcuR3JpZENlbGwgfCBudWxsID0gbnVsbDtcclxuXHJcbiAgY29uc3RydWN0b3IoY29sR3JpZCA6IERHLkdyaWRDb2x1bW4pIHtcclxuXHJcbiAgICBNb3VzZURpc3BhdGNoZXIuY3JlYXRlKCk7XHJcblxyXG4gICAgY29uc3QgZ3JpZCA9IGdldEdyaWQoY29sR3JpZCk7XHJcbiAgICBpZihncmlkID09PSBudWxsKSB7XHJcbiAgICAgIHRocm93IG5ldyBFcnJvcihcIkNvbHVtbiAnXCIgKyBjb2xHcmlkLm5hbWUgKyBcIicgaXMgbm90IGF0dGFjaGVkIHRvIHRoZSBncmlkLlwiKTtcclxuICAgIH1cclxuXHJcbiAgICBpZighUGlubmVkVXRpbHMuaXNQaW5uYWJsZUNvbHVtbihjb2xHcmlkKSkge1xyXG4gICAgICB0aHJvdyBuZXcgRXJyb3IoXCJDb2x1bW4gJ1wiICsgY29sR3JpZC5uYW1lICsgXCInIGNhbm5vdCBiZSBwaW5uZWQuIEl0IGVpdGhlciBwaW5uZWQgb3IgSFRNTC5cIik7XHJcbiAgICB9XHJcblxyXG4gICAgLy9sZXQgblJvd01pbiA9IGdyaWQubWluVmlzaWJsZVJvdztcclxuICAgIC8vbGV0IG5Sb3dNYXggPSBncmlkLm1heFZpc2libGVSb3c7XHJcbiAgICAvL2xldCBuQ29sTWluID0gZ3JpZC5taW5WaXNpYmxlQ29sdW1uO1xyXG4gICAgLy9sZXQgbkNvbE1heCA9IGdyaWQubWF4VmlzaWJsZUNvbHVtbjtcclxuXHJcblxyXG4gICAgdGhpcy5tX2ZEZXZpY2VQaXhlbFJhdGlvID0gd2luZG93LmRldmljZVBpeGVsUmF0aW87XHJcblxyXG4gICAgY29uc3QgZGFydCA9IERHLnRvRGFydChncmlkKTtcclxuXHJcbiAgICBpZihkYXJ0Lm1fYXJQaW5uZWRDb2xzID09PSB1bmRlZmluZWQpXHJcbiAgICAgIGRhcnQubV9hclBpbm5lZENvbHMgPSBbXTtcclxuXHJcbiAgICBpZihkYXJ0Lm1fYXJQaW5uZWRDb2xzLmxlbmd0aCA9PT0gMCAmJiAhR3JpZFV0aWxzLmlzUm93SGVhZGVyKGNvbEdyaWQpKSB7XHJcbiAgICAgIGNvbnN0IGNvbEdyaWQwID0gZ3JpZC5jb2x1bW5zLmJ5SW5kZXgoMCk7XHJcbiAgICAgIGlmKGNvbEdyaWQwICE9PSBudWxsICYmIGNvbEdyaWQwICE9PSB1bmRlZmluZWQpXHJcbiAgICAgIG5ldyBQaW5uZWRDb2x1bW4oY29sR3JpZDApO1xyXG4gICAgfVxyXG5cclxuICAgIGNvbnN0IG5XVG90YWxQaW5uZWRDb2xzID0gUGlubmVkVXRpbHMuZ2V0VG90YWxQaW5uZWRDb2xzV2lkdGgoZ3JpZCk7XHJcbiAgICBkYXJ0Lm1fYXJQaW5uZWRDb2xzLnB1c2godGhpcyk7XHJcblxyXG4gICAgY29uc3Qgdmlld1RhYmxlID0gZ3JpZC52aWV3O1xyXG4gICAgY29uc3QgZGZyYW1lID0gZ3JpZC5kYXRhRnJhbWU7XHJcblxyXG4gICAgY29uc3QgblcgPSBjb2xHcmlkLndpZHRoO1xyXG4gICAgdGhpcy5tX2NvbEdyaWQgPSBjb2xHcmlkO1xyXG4gICAgdGhpcy5tX25XaWR0aEJ1ZyA9IC0xO1xyXG4gICAgdHJ5IHtcclxuICAgICAgY29sR3JpZC52aXNpYmxlID0gZmFsc2U7XHJcbiAgICB9XHJcbiAgICBjYXRjaChlKSB7XHJcbiAgICAgIC8vREcgYnVnXHJcbiAgICAgIGNvbnNvbGUuZXJyb3IoXCJFUlJPUjogQ291bGRuJ3QgaGlkZSBjb2x1bW4gJ1wiICsgY29sR3JpZC5uYW1lICsgXCInIGR1ZSB0byBhIERHIGJ1Zy4gQXR0ZW1wdCB0byBzZXQgdGhlIHdpZHRoIHRvIDBcIik7XHJcbiAgICAgIHRyeSB7XHJcbiAgICAgICAgdGhpcy5tX25XaWR0aEJ1ZyA9IGNvbEdyaWQud2lkdGg7XHJcbiAgICAgICAgY29sR3JpZC53aWR0aCA9IDA7XHJcbiAgICAgIH0gY2F0Y2ggKGUpIHtcclxuICAgICAgICAvL0RHIGJ1Z1xyXG4gICAgICAgIGNvbnNvbGUuZXJyb3IoXCJFUlJPUjogQ291bGRuJ3Qgc2V0IHRoZSB3aWR0aCB0byAwIGZvciBjb2x1bW4gJ1wiICsgY29sR3JpZC5uYW1lICsgXCInIGR1ZSB0byBhIERHIGJ1Zy4gVGhpcyBjb3VsZCBiZSBpZ25vcmVkIGlmIHRoZSBjb2x1bW4gdmlzdWFsbHkgbG9va3Mgb2suXCIpO1xyXG4gICAgICB9XHJcbiAgICB9XHJcblxyXG4gICAgaWYoIUdyaWRVdGlscy5pc1Jvd0hlYWRlcihjb2xHcmlkKSkge1xyXG4gICAgICBpZiAoY29sR3JpZC5zZXR0aW5ncyA9PT0gbnVsbCB8fCBjb2xHcmlkLnNldHRpbmdzID09PSB1bmRlZmluZWQpXHJcbiAgICAgICAgY29sR3JpZC5zZXR0aW5ncyA9IHt9O1xyXG5cclxuICAgICAgY29sR3JpZC5zZXR0aW5ncy5pc1Bpbm5lZCA9IHRydWU7IC8vdGhpcyB3aWxsIGJlIHNhdmVkIHdpdGggdGhlIGxheW91dFxyXG4gICAgICBjb2xHcmlkLnNldHRpbmdzLmlkeFBpbm5lZCA9IGRhcnQubV9hclBpbm5lZENvbHMubGVuZ3RoIC0gMTtcclxuICAgIH1cclxuXHJcbiAgICBncmlkLmNhbnZhcy5zdHlsZS5sZWZ0ID0gKGdyaWQuY2FudmFzLm9mZnNldExlZnQgKyBuVykudG9TdHJpbmcoKSArIFwicHhcIjtcclxuICAgIGdyaWQub3ZlcmxheS5zdHlsZS5sZWZ0PSAoZ3JpZC5vdmVybGF5Lm9mZnNldExlZnQgKyBuVykudG9TdHJpbmcoKSArIFwicHhcIjtcclxuXHJcbiAgICBncmlkLmNhbnZhcy5zdHlsZS53aWR0aCA9IChncmlkLmNhbnZhcy5vZmZzZXRXaWR0aCAtIG5XKS50b1N0cmluZygpICsgXCJweFwiO1xyXG4gICAgZ3JpZC5vdmVybGF5LnN0eWxlLndpZHRoPSAoZ3JpZC5vdmVybGF5Lm9mZnNldFdpZHRoIC0gblcpLnRvU3RyaW5nKCkgKyBcInB4XCI7XHJcblxyXG4gICAgY29uc3QgbkhlaWdodCA9IGdyaWQuY2FudmFzLmhlaWdodDsvL2NhbnZhcyBwaXhlbCBoZWlnaHRcclxuICAgIGNvbnN0IGVDYW52YXNUaGlzID0gdWkuY2FudmFzKG5XKndpbmRvdy5kZXZpY2VQaXhlbFJhdGlvLCBuSGVpZ2h0KTtcclxuICAgIGNvbnN0IHRhYkluZGV4ID0gIGdyaWQuY2FudmFzLmdldEF0dHJpYnV0ZShcInRhYkluZGV4XCIpO1xyXG4gICAgaWYodGFiSW5kZXggIT09IG51bGwpXHJcbiAgICAgZUNhbnZhc1RoaXMuc2V0QXR0cmlidXRlKFwidGFiSW5kZXhcIiwgdGFiSW5kZXgpO1xyXG5cclxuICAgIGVDYW52YXNUaGlzLnN0eWxlLnBvc2l0aW9uID0gXCJhYnNvbHV0ZVwiO1xyXG4gICAgZUNhbnZhc1RoaXMuc3R5bGUubGVmdCA9IG5XVG90YWxQaW5uZWRDb2xzICsgXCJweFwiO1xyXG4gICAgZUNhbnZhc1RoaXMuc3R5bGUudG9wID0gZ3JpZC5jYW52YXMub2Zmc2V0VG9wICsgXCJweFwiO1xyXG4gICAgZUNhbnZhc1RoaXMuc3R5bGUud2lkdGggPSBuVyArIFwicHhcIjtcclxuICAgIGVDYW52YXNUaGlzLnN0eWxlLmhlaWdodCA9IE1hdGgucm91bmQobkhlaWdodC93aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbykgKyBcInB4XCI7XHJcblxyXG4gICAgLy9jb25zb2xlLmxvZyhcImggXCIgKyBncmlkLmNhbnZhcy5oZWlnaHQgKyBcIiBvZmZzZXQgXCIgKyBncmlkLmNhbnZhcy5vZmZzZXRIZWlnaHQpO1xyXG5cclxuICAgIGlmKGdyaWQuY2FudmFzLnBhcmVudE5vZGUgPT09IG51bGwpXHJcbiAgICAgIHRocm93IG5ldyBFcnJvcihcIlBhcmVudCBub2RlIGZvciBjYW52YXMgY2Fubm90IGJlIG51bGwuXCIpO1xyXG5cclxuICAgIGdyaWQuY2FudmFzLnBhcmVudE5vZGUuaW5zZXJ0QmVmb3JlKGVDYW52YXNUaGlzLCBncmlkLmNhbnZhcyk7XHJcbiAgICB0aGlzLm1fcm9vdCA9IGVDYW52YXNUaGlzO1xyXG5cclxuXHJcbiAgICBjb25zdCBjb2xHcmlkMCA9IGdyaWQuY29sdW1ucy5ieUluZGV4KDApO1xyXG4gICAgaWYoY29sR3JpZDAgIT09IG51bGwgJiYgY29sR3JpZDAgIT09IHVuZGVmaW5lZCkgey8vREcgQnVnIGZyb20gcmVhZGluZyBsYXlvdXRcclxuICAgIHRyeXtcclxuICAgICAgICBjb2xHcmlkMC52aXNpYmxlID0gZmFsc2U7XHJcbiAgICAgIH1cclxuICAgICAgY2F0Y2goZSkge1xyXG4gICAgICAgIGNvbnNvbGUuZXJyb3IoXCJFUlJPUjogQ291bGRuJ3QgaGlkZSByb3cgaGVhZGVyLlwiKTtcclxuICAgICAgfVxyXG4gICAgfVxyXG5cclxuXHJcbiAgICAvL09uUmVzaXplIFJvdyBoZWFkZXJcclxuICAgIGNvbnN0IGhlYWRlclRoaXMgPSB0aGlzOy8qXHJcbiAgICB0aGlzLm1fb2JzZXJ2ZXJSZXNpemUgPSBuZXcgUmVzaXplT2JzZXJ2ZXIoZW50cmllcyA9PiB7XHJcbiAgICAgIGNvbnN0IGcgPSBoZWFkZXJUaGlzLm1fcm9vdC5nZXRDb250ZXh0KCcyZCcpO1xyXG4gICAgICBmb3IgKGxldCBlbnRyeSBvZiBlbnRyaWVzKSB7XHJcbiAgICAgICAgaGVhZGVyVGhpcy5wYWludChnLCBncmlkKTtcclxuICAgICAgfVxyXG4gICAgfSk7XHJcbiAgICB0aGlzLm1fb2JzZXJ2ZXJSZXNpemUub2JzZXJ2ZShoZWFkZXJUaGlzLm1fcm9vdCk7Ki9cclxuXHJcblxyXG5cclxuICAgIC8vT25SZXNpemUgR3JpZFxyXG4gICAgdGhpcy5tX29ic2VydmVyUmVzaXplR3JpZCA9IG5ldyBSZXNpemVPYnNlcnZlcihmdW5jdGlvbiAoZW50cmllcyA6IGFueSkge1xyXG5cclxuICAgICAgY29uc3QgYkN1cnJlbnQgPSAgREcudG9EYXJ0KGdyb2suc2hlbGwudikgPT09IERHLnRvRGFydCh2aWV3VGFibGUpO1xyXG4gICAgICBpZighYkN1cnJlbnQpXHJcbiAgICAgICAgcmV0dXJuO1xyXG5cclxuICAgICAgaWYoaGVhZGVyVGhpcy5tX2ZEZXZpY2VQaXhlbFJhdGlvICE9PSB3aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbyB8fCBncmlkLmNhbnZhcy5oZWlnaHQgIT09IGVDYW52YXNUaGlzLmhlaWdodCkge1xyXG4gICAgICAgIGVDYW52YXNUaGlzLndpZHRoID0gblcqd2luZG93LmRldmljZVBpeGVsUmF0aW87XHJcbiAgICAgICAgZUNhbnZhc1RoaXMuaGVpZ2h0ID0gZ3JpZC5jYW52YXMuaGVpZ2h0O1xyXG4gICAgICAgIGVDYW52YXNUaGlzLnN0eWxlLnRvcCA9IGdyaWQuY2FudmFzLm9mZnNldFRvcCArIFwicHhcIjtcclxuICAgICAgICBlQ2FudmFzVGhpcy5zdHlsZS53aWR0aCA9IG5XICsgXCJweFwiO1xyXG4gICAgICAgIGVDYW52YXNUaGlzLnN0eWxlLmhlaWdodCA9IE1hdGgucm91bmQoZ3JpZC5jYW52YXMuaGVpZ2h0L3dpbmRvdy5kZXZpY2VQaXhlbFJhdGlvKSArIFwicHhcIjtcclxuXHJcbiAgICAgICAgaGVhZGVyVGhpcy5tX2ZEZXZpY2VQaXhlbFJhdGlvID0gd2luZG93LmRldmljZVBpeGVsUmF0aW87XHJcbiAgICAgIH1cclxuXHJcbiAgICAgIC8vY29uc29sZS5sb2coXCJHcmlkIFJlc2l6ZTogXCIgKyBncmlkLmNhbnZhcy5oZWlnaHQgKyBcIiBcIiArIHdpbmRvdy5kZXZpY2VQaXhlbFJhdGlvKTtcclxuICAgICAgLy9lQ2FudmFzVGhpcy5zdHlsZS5oZWlnaHQgPSBncmlkLnJvb3Quc3R5bGUuaGVpZ2h0O1xyXG4vKlxyXG4gICAgICBjb25zdCBlQ2FudmFzTmV3ID0gdWkuY2FudmFzKG5XLCBncmlkLnJvb3Qub2Zmc2V0SGVpZ2h0KTtcclxuICAgICAgaWYoaGVhZGVyVGhpcy5tX3Jvb3QucGFyZW50Tm9kZSAhPT0gbnVsbCkge1xyXG4gICAgICAgIGhlYWRlclRoaXMubV9yb290LnBhcmVudE5vZGUucmVwbGFjZUNoaWxkKGVDYW52YXNOZXcsIGhlYWRlclRoaXMubV9yb290KTtcclxuICAgICAgICBoZWFkZXJUaGlzLm1fcm9vdCA9IGVDYW52YXNOZXc7XHJcbiAgICAgIH0qL1xyXG4gICAgICAvL2hlYWRlclRoaXMubV9yb290LmhlaWdodCA9IGdyaWQucm9vdC5vZmZzZXRIZWlnaHQ7XHJcbiAgICAgIGNvbnN0IGcgPSBlQ2FudmFzVGhpcy5nZXRDb250ZXh0KCcyZCcpO1xyXG4gICAgICBmb3IgKGxldCBlbnRyeSBvZiBlbnRyaWVzKSB7XHJcbiAgICAgICAgc2V0VGltZW91dCgoKT0+IHtoZWFkZXJUaGlzLnBhaW50KGcsIGdyaWQpO30sIDEwMCk7XHJcbiAgICAgIH1cclxuICAgIH0pO1xyXG5cclxuICAgIHRoaXMubV9vYnNlcnZlclJlc2l6ZUdyaWQ/Lm9ic2VydmUoZ3JpZC5jYW52YXMpO1xyXG5cclxuICAgIGNvbnN0IHNjcm9sbFZlcnQgPSBncmlkLnZlcnRTY3JvbGw7XHJcbiAgICB0aGlzLm1faGFuZGxlclZTY3JvbGwgPSBzY3JvbGxWZXJ0Lm9uVmFsdWVzQ2hhbmdlZC5zdWJzY3JpYmUoKCkgPT4ge1xyXG4gICAgICBjb25zdCBnID0gZUNhbnZhc1RoaXMuZ2V0Q29udGV4dCgnMmQnKTtcclxuICAgICAgaGVhZGVyVGhpcy5wYWludChnLCBncmlkKTtcclxuICAgIH0pO1xyXG5cclxuICAgIHRoaXMubV9oYW5kbGVyUm93c0ZpbHRlcmluZyA9IGRmcmFtZS5vblJvd3NGaWx0ZXJpbmcuc3Vic2NyaWJlKCgpID0+IHtcclxuICAgICAgc2V0VGltZW91dCgoKSA9PiB7XHJcbiAgICAgICAgY29uc3QgZyA9IGVDYW52YXNUaGlzLmdldENvbnRleHQoJzJkJyk7XHJcbiAgICAgICAgaGVhZGVyVGhpcy5wYWludChnLCBncmlkKTtcclxuICAgICAgfSwgMTAwKTtcclxuXHJcbiAgICB9KTtcclxuXHJcbiAgICB0aGlzLm1faGFuZGxlckN1cnJSb3cgPSBkZnJhbWUub25DdXJyZW50Um93Q2hhbmdlZC5zdWJzY3JpYmUoKCkgPT4ge1xyXG4gICAgICAgIGNvbnN0IGcgPSBlQ2FudmFzVGhpcy5nZXRDb250ZXh0KCcyZCcpO1xyXG4gICAgICAgIGhlYWRlclRoaXMucGFpbnQoZywgZ3JpZCk7XHJcbiAgICAgIH1cclxuICAgICk7XHJcblxyXG4gICAgdGhpcy5tX2hhbmRsZXJTZWwgPSBkZnJhbWUub25TZWxlY3Rpb25DaGFuZ2VkLnN1YnNjcmliZSgoZSA6IGFueSkgPT4ge1xyXG4gICAgICAgIGNvbnN0IGcgPSBlQ2FudmFzVGhpcy5nZXRDb250ZXh0KCcyZCcpO1xyXG4gICAgICAgIGhlYWRlclRoaXMucGFpbnQoZywgZ3JpZCk7XHJcbiAgICAgIH1cclxuICAgICk7XHJcblxyXG4gICAgdGhpcy5tX2hhbmRsZXJDb2xzUmVtb3ZlZCA9IGRmcmFtZS5vbkNvbHVtbnNSZW1vdmVkLnN1YnNjcmliZSgoZSA6IENvbHVtbnNBcmdzKSA9PiB7XHJcblxyXG4gICAgICAgICAgaWYoaGVhZGVyVGhpcy5tX2NvbEdyaWQgPT09IG51bGwpXHJcbiAgICAgICAgICAgIHJldHVybjtcclxuICAgICAgICAgIGZvcihsZXQgbkM9MDsgbkM8ZS5jb2x1bW5zLmxlbmd0aDsgKytuQykge1xyXG4gICAgICAgICAgICBpZihlLmNvbHVtbnNbbkNdLm5hbWUgPT09IGhlYWRlclRoaXMubV9jb2xHcmlkLm5hbWUpXHJcbiAgICAgICAgICAgICAgaGVhZGVyVGhpcy5jbG9zZSgpO1xyXG4gICAgICAgICAgfVxyXG4gICAgICAgIH1cclxuICAgICk7XHJcblxyXG4gICAgdGhpcy5tX2hhbmRsZXJDb2xOYW1lQ2hhbmdlZCA9IGRmcmFtZS5vbkNvbHVtbk5hbWVDaGFuZ2VkLnN1YnNjcmliZSgoZSA6IGFueSkgPT4ge1xyXG5cclxuICAgICAgICAgIGNvbnN0IGRhcnQgPSB0b0RhcnQoZSk7XHJcbiAgICAgICAgICBjb25zdCBzdHJDb2xOYW1lT2xkID0gZGFydC5uZXdOYW1lO1xyXG4gICAgICAgICAgaWYoc3RyQ29sTmFtZU9sZCA9PT0gaGVhZGVyVGhpcy5tX2NvbEdyaWQ/Lm5hbWUpIHtcclxuICAgICAgICAgICAgY29uc3QgZyA9IGVDYW52YXNUaGlzLmdldENvbnRleHQoJzJkJyk7XHJcbiAgICAgICAgICAgIGhlYWRlclRoaXMucGFpbnQoZywgZ3JpZCk7XHJcbiAgICAgICAgICB9XHJcbiAgICAgICAgfVxyXG4gICAgKTtcclxuXHJcblxyXG4vKlxyXG4gICAgdGhpcy5tX2hhbmRsZXJGaWx0ZXIgPSBkZnJhbWUub25Sb3dzRmlsdGVyZWQuc3Vic2NyaWJlKChlIDogYW55KSA9PiB7XHJcbiAgICAgICAgY29uc3QgZyA9IGVDYW52YXNUaGlzLmdldENvbnRleHQoJzJkJyk7XHJcbiAgICAgICAgaGVhZGVyVGhpcy5wYWludChnLCBncmlkKTtcclxuICAgICAgfVxyXG4gICAgKTtcclxuKi9cclxuXHJcbiAgICB0aGlzLm1faGFuZGxlclJvd3NSZXNpemVkID0gZ3JpZC5vblJvd3NSZXNpemVkLnN1YnNjcmliZSgoZSA6IGFueSkgPT4ge1xyXG4gICAgICAgIGNvbnN0IGcgPSBlQ2FudmFzVGhpcy5nZXRDb250ZXh0KCcyZCcpO1xyXG4gICAgICAgIGhlYWRlclRoaXMucGFpbnQoZywgZ3JpZCk7XHJcbiAgICAgIH1cclxuICAgICk7XHJcblxyXG4gICAgdGhpcy5tX2hhbmRsZXJSb3dzU29ydGVkID0gZ3JpZC5vblJvd3NTb3J0ZWQuc3Vic2NyaWJlKChlIDogYW55KSA9PiB7XHJcbiAgICAgICAgY29uc3QgZyA9IGVDYW52YXNUaGlzLmdldENvbnRleHQoJzJkJyk7XHJcbiAgICAgICAgaGVhZGVyVGhpcy5wYWludChnLCBncmlkKTtcclxuICAgICAgfVxyXG4gICAgKTtcclxuICB9XHJcblxyXG4gIGlzUGlubmVkKCkgOiBib29sZWFuIHtcclxuICAgIHJldHVybiB0aGlzLm1fY29sR3JpZCAhPT0gbnVsbDtcclxuICB9XHJcblxyXG4gIGdldEdyaWRDb2x1bW4oKSA6IERHLkdyaWRDb2x1bW4gfCBudWxse1xyXG4gICAgcmV0dXJuIHRoaXMubV9jb2xHcmlkO1xyXG4gIH1cclxuXHJcbiAgZ2V0V2lkdGgoKSA6IG51bWJlciB7XHJcbiAgICByZXR1cm4gdGhpcy5tX3Jvb3QgPT09IG51bGwgPyAtMSA6IHRoaXMubV9yb290Lm9mZnNldFdpZHRoO1xyXG4gIH1cclxuXHJcbiAgZ2V0Um9vdCgpIDogSFRNTENhbnZhc0VsZW1lbnQgfCBudWxsIHtcclxuICAgIHJldHVybiB0aGlzLm1fcm9vdDtcclxuICB9XHJcblxyXG4gIHB1YmxpYyBjbG9zZSgpIDogdm9pZCB7XHJcblxyXG4gICAgaWYodGhpcy5tX2NvbEdyaWQgPT09IG51bGwpIHtcclxuICAgICAgdGhyb3cgbmV3IEVycm9yKFwiQ29sdW1uIGhhcyBhbHJlYWR5IGJlZW4gdW5waW5uZWRcIik7XHJcbiAgICB9XHJcblxyXG4gICAgaWYodGhpcy5tX29ic2VydmVyUmVzaXplR3JpZCAhPT0gbnVsbCkge1xyXG4gICAgICB0aGlzLm1fb2JzZXJ2ZXJSZXNpemVHcmlkLmRpc2Nvbm5lY3QoKTtcclxuICAgICAgdGhpcy5tX29ic2VydmVyUmVzaXplR3JpZCA9IG51bGw7XHJcbiAgICB9XHJcbi8qbXkgY2hhbmdlc1xyXG4gICAgaWYodGhpcy5tX29ic2VydmVyUmVzaXplICE9PSBudWxsKSB7XHJcbiAgICAgIHRoaXMubV9vYnNlcnZlclJlc2l6ZS5kaXNjb25uZWN0KCk7XHJcbiAgICAgIHRoaXMubV9vYnNlcnZlclJlc2l6ZSA9IG51bGw7XHJcbiAgICB9XHJcbiAgICAqL1xyXG5cclxuICAgIHRoaXMubV9oYW5kbGVyQ29sc1JlbW92ZWQudW5zdWJzY3JpYmUoKTtcclxuICAgIHRoaXMubV9oYW5kbGVyQ29sc1JlbW92ZWQgPSBudWxsO1xyXG5cclxuICAgIHRoaXMubV9oYW5kbGVyQ29sTmFtZUNoYW5nZWQudW5zdWJzY3JpYmUoKTtcclxuICAgIHRoaXMubV9oYW5kbGVyQ29sTmFtZUNoYW5nZWQgPSBudWxsO1xyXG5cclxuICAgIHRoaXMubV9oYW5kbGVyVlNjcm9sbC51bnN1YnNjcmliZSgpO1xyXG4gICAgdGhpcy5tX2hhbmRsZXJWU2Nyb2xsID0gbnVsbDtcclxuXHJcbiAgICB0aGlzLm1faGFuZGxlclJvd3NSZXNpemVkLnVuc3Vic2NyaWJlKCk7XHJcbiAgICB0aGlzLm1faGFuZGxlclJvd3NSZXNpemVkID0gbnVsbDtcclxuXHJcbiAgICB0aGlzLm1faGFuZGxlclJvd3NTb3J0ZWQudW5zdWJzY3JpYmUoKTtcclxuICAgIHRoaXMubV9oYW5kbGVyUm93c1NvcnRlZCA9IG51bGw7XHJcblxyXG4gICAgdGhpcy5tX2hhbmRsZXJSb3dzRmlsdGVyaW5nLnVuc3Vic2NyaWJlKCk7XHJcbiAgICB0aGlzLm1faGFuZGxlclJvd3NGaWx0ZXJpbmcgPSBudWxsO1xyXG5cclxuICAgIHRoaXMubV9oYW5kbGVyQ3VyclJvdy51bnN1YnNjcmliZSgpO1xyXG4gICAgdGhpcy5tX2hhbmRsZXJDdXJyUm93ID0gbnVsbDtcclxuXHJcbiAgICB0aGlzLm1faGFuZGxlclNlbC51bnN1YnNjcmliZSgpO1xyXG4gICAgdGhpcy5tX2hhbmRsZXJTZWwgPSBudWxsO1xyXG5cclxuICAgIGNvbnN0IGdyaWQgPSBnZXRHcmlkKHRoaXMubV9jb2xHcmlkKTtcclxuICAgIGlmKGdyaWQgPT09IG51bGwpe1xyXG4gICAgICB0aHJvdyBuZXcgRXJyb3IoXCJDb2x1bW4gJ1wiICsgdGhpcy5tX2NvbEdyaWQubmFtZSArIFwiJyBpcyBkaXNjb25uZWN0ZWQgZnJvbSBncmlkLlwiKTtcclxuICAgIH1cclxuXHJcbiAgICBjb25zdCBkYXJ0ID0gREcudG9EYXJ0KGdyaWQpO1xyXG4gICAgY29uc3QgYXIgPSBkYXJ0Lm1fYXJQaW5uZWRDb2xzO1xyXG4gICAgY29uc3QgbklkeCA9IGFyLmluZGV4T2YodGhpcyk7XHJcbiAgICBhci5zcGxpY2UobklkeCwgMSk7XHJcblxyXG4gICAgaWYodGhpcy5tX3Jvb3QgPT09IG51bGwpXHJcbiAgICAgIHRocm93IG5ldyBFcnJvcignUm9vdCBjYW5ub3QgYmUgbnVsbCcpO1xyXG5cclxuICAgIGxldCBuSWR4UGlubmVkID0gLTE7XHJcbiAgICBsZXQgY29sR3JpZFRtcD0gbnVsbDtcclxuICAgIGZvcihsZXQgbj1uSWR4OyBuPGFyLmxlbmd0aDsgKytuKSB7XHJcbiAgICAgIGNvbEdyaWRUbXAgPSBhcltuXTtcclxuICAgICAgY29sR3JpZFRtcC5tX3Jvb3Quc3R5bGUubGVmdCA9IChjb2xHcmlkVG1wLm1fcm9vdC5vZmZzZXRMZWZ0IC0gdGhpcy5tX3Jvb3Qub2Zmc2V0V2lkdGgpLnRvU3RyaW5nKCkgKyBcInB4XCI7XHJcblxyXG4gICAgICBuSWR4UGlubmVkID0gIGNvbEdyaWRUbXAubV9jb2xHcmlkLnNldHRpbmdzLmlkeFBpbm5lZDtcclxuICAgICAgY29sR3JpZFRtcC5tX2NvbEdyaWQuc2V0dGluZ3MuaWR4UGlubmVkID0gbjtcclxuICAgIH1cclxuXHJcbiAgICBpZighR3JpZFV0aWxzLmlzUm93SGVhZGVyKHRoaXMubV9jb2xHcmlkKSkge1xyXG4gICAgICB0aGlzLm1fY29sR3JpZC5zZXR0aW5ncy5pZHhQaW5uZWQgPSAtMTtcclxuICAgICAgdGhpcy5tX2NvbEdyaWQuc2V0dGluZ3MuaXNQaW5uZWQgPSBmYWxzZTtcclxuICAgIH1cclxuXHJcblxyXG4gICAgaWYodGhpcy5tX25XaWR0aEJ1ZyA+PSAwKSB7XHJcbiAgICAgIHRyeSB7XHJcbiAgICAgICAgdGhpcy5tX2NvbEdyaWQud2lkdGggPSB0aGlzLm1fbldpZHRoQnVnO1xyXG4gICAgICB9XHJcbiAgICAgIGNhdGNoKGUpIHtcclxuICAgICAgICAvL0RHIGJ1Z1xyXG4gICAgICAgIGNvbnNvbGUuZXJyb3IoXCJFUlJPUjogQ291bGRuJ3Qgc2V0IHRoZSB3aWR0aCB0byBcIiArIHRoaXMubV9uV2lkdGhCdWcgKyBcIiBmb3IgY29sdW1uICdcIiArIHRoaXMubV9jb2xHcmlkLm5hbWUgKyBcIicgZHVlIHRvIGEgREcgYnVnLiBUaGlzIGNvdWxkIGJlIGlnbm9yZWQgaWYgdGhlIGNvbHVtbiB2aXN1YWxseSBsb29rcyBvay5cIik7XHJcbiAgICAgIH1cclxuICAgIH1cclxuXHJcbiAgICB0cnkge1xyXG4gICAgICB0aGlzLm1fY29sR3JpZC52aXNpYmxlID0gdHJ1ZTtcclxuICAgIH1cclxuICAgIGNhdGNoKGUpIHtcclxuICAgICAgLy9ERyBidWdcclxuICAgICAgY29uc29sZS5lcnJvcihcIkVSUk9SOiBDb3VsZG4ndCBzaG93IGNvbHVtbiAnXCIgKyB0aGlzLm1fY29sR3JpZC5uYW1lICsgXCInIGR1ZSB0byBhIERHIGJ1Zy4gVGhpcyBjb3VsZCBiZSBpZ25vcmVkIGlmIHRoZSBjb2x1bW4gdmlzdWFsbHkgbG9va3Mgb2suXCIpO1xyXG4gICAgfVxyXG5cclxuICAgIGdyaWQuY2FudmFzLnN0eWxlLmxlZnQgPSAoZ3JpZC5jYW52YXMub2Zmc2V0TGVmdCAtIHRoaXMubV9yb290Lm9mZnNldFdpZHRoKS50b1N0cmluZygpICsgXCJweFwiO1xyXG4gICAgZ3JpZC5vdmVybGF5LnN0eWxlLmxlZnQ9IChncmlkLm92ZXJsYXkub2Zmc2V0TGVmdCAtIHRoaXMubV9yb290Lm9mZnNldFdpZHRoKS50b1N0cmluZygpICsgXCJweFwiO1xyXG4gICAgZ3JpZC5jYW52YXMuc3R5bGUud2lkdGggPSAoZ3JpZC5jYW52YXMub2Zmc2V0V2lkdGggKyB0aGlzLm1fcm9vdC5vZmZzZXRXaWR0aCkudG9TdHJpbmcoKSArIFwicHhcIjtcclxuICAgIGdyaWQub3ZlcmxheS5zdHlsZS53aWR0aD0gKGdyaWQub3ZlcmxheS5vZmZzZXRXaWR0aCArIHRoaXMubV9yb290Lm9mZnNldFdpZHRoKS50b1N0cmluZygpICsgXCJweFwiO1xyXG5cclxuICAgIGlmKHRoaXMubV9yb290LnBhcmVudE5vZGUgIT09IG51bGwpXHJcbiAgICAgdGhpcy5tX3Jvb3QucGFyZW50Tm9kZS5yZW1vdmVDaGlsZCh0aGlzLm1fcm9vdCk7XHJcblxyXG4gICAgdGhpcy5tX3Jvb3QgPSBudWxsO1xyXG5cclxuICAgIGlmIChkYXJ0Lm1fYXJQaW5uZWRDb2xzLmxlbmd0aCA9PT0gMSAmJiBkYXJ0Lm1fYXJQaW5uZWRDb2xzWzBdLm1fY29sR3JpZC5pZHggPT09IDAgJiYgdGhpcy5tX2NvbEdyaWQuaWR4ICE9PSAwKSB7XHJcblxyXG4gICAgICAgIC8vIHRyeXtjb2xHcmlkMC52aXNpYmxlID0gdHJ1ZTt9XHJcbiAgICAgICAgdHJ5IHtcclxuICAgICAgICAgIGRhcnQubV9hclBpbm5lZENvbHNbMF0uY2xvc2UoKTtcclxuICAgICAgICB9IGNhdGNoIChlKSB7XHJcbiAgICAgICAgICBjb25zb2xlLmVycm9yKFwiRVJST1I6IENvdWxkbid0IGNsb3NlIHBpbm5lZCBjb2x1bW4gJ1wiICsgZGFydC5tX2FyUGlubmVkQ29sc1swXS5tX2NvbEdyaWQubmFtZSArIFwiJyBcIik7XHJcbiAgICAgICAgfVxyXG4gICAgfVxyXG4gICAgdGhpcy5tX2NvbEdyaWQgPSBudWxsO1xyXG4gIH1cclxuXHJcblxyXG4gIHB1YmxpYyBvbk1vdXNlRW50ZXIoZSA6IE1vdXNlRXZlbnQpIDogdm9pZCB7XHJcbiAgICBpZihERUJVRylcclxuICAgICAgY29uc29sZS5sb2coJ01vdXNlIEVudGVyIFBpbm5lZCBDb2x1bW46ICcgKyB0aGlzLmdldEdyaWRDb2x1bW4oKT8ubmFtZSk7XHJcbiAgfVxyXG5cclxuICBwdWJsaWMgb25Nb3VzZU1vdmUoZSA6IE1vdXNlRXZlbnQpIDogdm9pZCB7XHJcbiAgICBpZihERUJVRylcclxuICAgICAgY29uc29sZS5sb2coJ01vdXNlIE1vdmUgUGlubmVkIENvbHVtbjogJyArIHRoaXMuZ2V0R3JpZENvbHVtbigpPy5uYW1lKTtcclxuXHJcbiAgICBpZih0aGlzLm1fY29sR3JpZCA9PT0gbnVsbCB8fCB0aGlzLm1fcm9vdCA9PT0gbnVsbClcclxuICAgICAgcmV0dXJuO1xyXG5cclxuICAgIGNvbnN0IGdyaWQgPSB0aGlzLm1fY29sR3JpZC5ncmlkO1xyXG4gICAgY29uc3Qgdmlld1RhYmxlID0gZ3JpZC52aWV3O1xyXG5cclxuICAgIGlmKERHLnRvRGFydChncm9rLnNoZWxsLnYpICE9PSBERy50b0RhcnQodmlld1RhYmxlKSkge1xyXG4gICAgICByZXR1cm47XHJcbiAgICB9XHJcblxyXG5cclxuICAgIGNvbnN0IGFyWFlPbkNlbGwgPSBbLTEsLTFdO1xyXG5cclxuICAgIGxldCBuUm93R3JpZCA9IFBpbm5lZENvbHVtbi5oaXRUZXN0Um93cyh0aGlzLm1fcm9vdCwgZ3JpZCwgZSwgZmFsc2UsIGFyWFlPbkNlbGwpO1xyXG4gICAgaWYoblJvd0dyaWQgPj0gMCkge1xyXG4gICAgICBjb25zdCBjZWxsID0gZ3JpZC5jZWxsKHRoaXMubV9jb2xHcmlkLm5hbWUsIG5Sb3dHcmlkKTtcclxuICAgICAgY29uc3QgcmVuZGVyZXIgPSBnZXRSZW5kZXJlcihjZWxsKTtcclxuXHJcbiAgICAgIGlmIChyZW5kZXJlciBpbnN0YW5jZW9mIEdyaWRDZWxsUmVuZGVyZXJFeCkge1xyXG5cclxuICAgICAgICBpZiAodGhpcy5tX2NlbGxDdXJyZW50ID09PSBudWxsKSB7XHJcbiAgICAgICAgICByZW5kZXJlci5vbk1vdXNlRW50ZXJFeChjZWxsLCBlLCBhclhZT25DZWxsWzBdLCBhclhZT25DZWxsWzFdKTtcclxuICAgICAgICB9XHJcblxyXG4gICAgICAgIGlmICh0aGlzLm1fY2VsbEN1cnJlbnQgIT09IG51bGwgJiYgblJvd0dyaWQgIT09IHRoaXMubV9jZWxsQ3VycmVudC5ncmlkUm93KSB7XHJcbiAgICAgICAgICByZW5kZXJlci5vbk1vdXNlTGVhdmVFeCh0aGlzLm1fY2VsbEN1cnJlbnQsIGUsIC0xLCAtMSk7XHJcblxyXG4gICAgICAgICAgcmVuZGVyZXIub25Nb3VzZUVudGVyRXgoY2VsbCwgZSwgYXJYWU9uQ2VsbFswXSwgYXJYWU9uQ2VsbFsxXSk7XHJcbiAgICAgICAgfVxyXG5cclxuICAgICAgICByZW5kZXJlci5vbk1vdXNlTW92ZUV4KGNlbGwsIGUsIGFyWFlPbkNlbGxbMF0sIGFyWFlPbkNlbGxbMV0pO1xyXG4gICAgICB9XHJcblxyXG4gICAgICB0aGlzLm1fY2VsbEN1cnJlbnQgPSBjZWxsO1xyXG4gICAgfVxyXG4gICAgZWxzZSBpZiAodGhpcy5tX2NlbGxDdXJyZW50ICE9PSBudWxsKSB7XHJcbiAgICAgIGNvbnN0IHJlbmRlcmVyID0gZ2V0UmVuZGVyZXIodGhpcy5tX2NlbGxDdXJyZW50KTtcclxuICAgICAgaWYgKHJlbmRlcmVyIGluc3RhbmNlb2YgR3JpZENlbGxSZW5kZXJlckV4KSB7XHJcbiAgICAgICAgcmVuZGVyZXIub25Nb3VzZUxlYXZlRXgodGhpcy5tX2NlbGxDdXJyZW50LCBlLCAtMSwgLTEpO1xyXG4gICAgICB9XHJcblxyXG4gICAgICB0aGlzLm1fY2VsbEN1cnJlbnQgPSBudWxsO1xyXG4gICAgfVxyXG5cclxuICAgIG5Sb3dHcmlkID0gUGlubmVkQ29sdW1uLmhpdFRlc3RSb3dzKHRoaXMubV9yb290LCBncmlkLCBlLCB0cnVlLCB1bmRlZmluZWQpO1xyXG4gICAgaWYgKG5Sb3dHcmlkID49IDApIHtcclxuICAgICAgdGhpcy5tX25SZXNpemVSb3dHcmlkTW92aW5nID0gblJvd0dyaWQ7XHJcbiAgICAgIGRvY3VtZW50LmJvZHkuc3R5bGUuY3Vyc29yID0gXCJyb3ctcmVzaXplXCI7XHJcbiAgICAgIHJldHVybjtcclxuICAgIH1cclxuXHJcbiAgICBpZih0aGlzLm1fblJlc2l6ZVJvd0dyaWRNb3ZpbmcgPj0gMCkge1xyXG4gICAgICB0aGlzLm1fblJlc2l6ZVJvd0dyaWRNb3ZpbmcgPSAtMTtcclxuICAgICAgZG9jdW1lbnQuYm9keS5zdHlsZS5jdXJzb3IgPSBcImF1dG9cIjtcclxuICAgIH1cclxuXHJcblxyXG4gICAgLy9IYW1idXJnZXIgTWVudVxyXG4gICAgY29uc3QgY29sR3JpZCA9IHRoaXMuZ2V0R3JpZENvbHVtbigpO1xyXG4gICAgaWYoY29sR3JpZCA9PT0gbnVsbCB8fCBjb2xHcmlkLm5hbWUgPT09ICcnKVxyXG4gICAgICByZXR1cm47XHJcblxyXG4gICAgY29uc3QgZURpdkhhbWIgPSBHcmlkVXRpbHMuZ2V0VG9vbEljb25EaXYoY29sR3JpZC5ncmlkKTtcclxuICAgIGNvbnN0IG5IQ29sSGVhZGVyID0gR3JpZFV0aWxzLmdldEdyaWRDb2x1bW5IZWFkZXJIZWlnaHQoY29sR3JpZC5ncmlkKTtcclxuICAgIGlmKDAgPD0gZS5vZmZzZXRZICYmIGUub2Zmc2V0WSA8IG5IQ29sSGVhZGVyKSB7XHJcblxyXG4gICAgICBlRGl2SGFtYj8uc3R5bGUucmVtb3ZlUHJvcGVydHkoJ3Zpc2liaWxpdHknKTtcclxuICAgICAgZURpdkhhbWI/LnNldEF0dHJpYnV0ZSgnY29sdW1uX25hbWUnLCBjb2xHcmlkLm5hbWUpO1xyXG4gICAgICAvL2NvbnNvbGUubG9nKCdUb29sc0ljb24gZm9yIGNvbHVtbiAnICsgY29sR3JpZC5uYW1lKTtcclxuICAgICAgLy8gQHRzLWlnbm9yZVxyXG4gICAgICBlRGl2SGFtYj8uc3R5bGUubGVmdCA9IChQaW5uZWRVdGlscy5nZXRQaW5uZWRDb2x1bW5MZWZ0KHRoaXMpICsgdGhpcy5nZXRXaWR0aCgpIC0gMTgpICsgJ3B4JztcclxuICAgICAgLy8gQHRzLWlnbm9yZVxyXG4gICAgICBlRGl2SGFtYj8uc3R5bGUudG9wID0gKEdyaWRVdGlscy5nZXRHcmlkQ29sdW1uSGVhZGVySGVpZ2h0KGNvbEdyaWQuZ3JpZCkgLSAxNikgKyBcInB4XCI7XHJcbiAgICB9IGVsc2Uge1xyXG4gICAgICBjb25zdCBjb2xHcmlkID0gdGhpcy5nZXRHcmlkQ29sdW1uKCk7XHJcbiAgICAgIGlmKGNvbEdyaWQgIT0gbnVsbCkge1xyXG4gICAgICAgICAgZURpdkhhbWI/LnNldEF0dHJpYnV0ZSgnY29sdW1uX25hbWUnLCAnJyk7XHJcbiAgICAgICAgICAvLyBAdHMtaWdub3JlXHJcbiAgICAgICAgICBlRGl2SGFtYj8uc3R5bGUudmlzaWJpbGl0eSA9ICdoaWRkZW4nO1xyXG4gICAgICAgIH1cclxuICAgIH1cclxuICB9XHJcblxyXG4gIHB1YmxpYyBvbk1vdXNlRHJhZyhlIDogTW91c2VFdmVudCkgOiB2b2lkIHtcclxuICAgIGlmKERFQlVHKVxyXG4gICAgIGNvbnNvbGUubG9nKCdNb3VzZSBEcmFnIFBpbm5lZCBDb2x1bW46ICcgKyB0aGlzLmdldEdyaWRDb2x1bW4oKT8ubmFtZSk7XHJcblxyXG4gICAgaWYodGhpcy5tX2NvbEdyaWQgPT09IG51bGwgfHwgdGhpcy5tX3Jvb3QgPT09IG51bGwpXHJcbiAgICByZXR1cm47XHJcblxyXG4gICAgY29uc3QgZ3JpZCA9IHRoaXMubV9jb2xHcmlkLmdyaWQ7XHJcbiAgICBjb25zdCB2aWV3VGFibGUgPSBncmlkLnZpZXc7XHJcblxyXG4gICAgaWYoREcudG9EYXJ0KGdyb2suc2hlbGwudikgIT09IERHLnRvRGFydCh2aWV3VGFibGUpKSB7XHJcbiAgICAgIHJldHVybjtcclxuICAgIH1cclxuXHJcbiAgICBjb25zdCBiUmVzaXppbmcgPSB0aGlzLm1fblJlc2l6ZVJvd0dyaWREcmFnZ2luZyA+PSAwO1xyXG4gICAgaWYgKGJSZXNpemluZykge1xyXG5cclxuICAgICAgLy9jb25zb2xlLmxvZyhcIkRyYWdnaW5nIDogXCIgKyBoZWFkZXJUaGlzLm1fc3RyQ29sTmFtZSk7XHJcbiAgICAgIGNvbnN0IG5ZRGlmZiA9IGUuY2xpZW50WSAtIHRoaXMubV9uWVJlc2l6ZURyYWdnaW5nQW5jaG9yO1xyXG4gICAgICBsZXQgbkhSb3dHcmlkID0gdGhpcy5tX25IUmVzaXplUm93c0JlZm9yZURyYWcgKyBuWURpZmY7XHJcblxyXG4gICAgICBpZiAobkhSb3dHcmlkIDwgUGlubmVkQ29sdW1uLk1JTl9ST1dfSEVJR0hUKVxyXG4gICAgICAgIG5IUm93R3JpZCA9IFBpbm5lZENvbHVtbi5NSU5fUk9XX0hFSUdIVDtcclxuICAgICAgZWxzZSBpZiAobkhSb3dHcmlkID4gUGlubmVkQ29sdW1uLk1BWF9ST1dfSEVJR0hUKVxyXG4gICAgICAgIG5IUm93R3JpZCA9IFBpbm5lZENvbHVtbi5NQVhfUk9XX0hFSUdIVDtcclxuXHJcbiAgICAgIGNvbnN0IGVDYW52YXNUaGlzID0gdGhpcy5tX3Jvb3Q7XHJcblxyXG4gICAgICBsZXQgZyA9IGVDYW52YXNUaGlzLmdldENvbnRleHQoJzJkJyk7XHJcbiAgICAgIGlmKGcgPT09IG51bGwpXHJcbiAgICAgICAgcmV0dXJuO1xyXG5cclxuICAgICAgZy5maWxsU3R5bGUgPSBcIndoaXRlXCI7XHJcbiAgICAgIGNvbnN0IG5ISGVhZGVyQ29scyA9IEdyaWRVdGlscy5nZXRHcmlkQ29sdW1uSGVhZGVySGVpZ2h0KGdyaWQpO1xyXG4gICAgICBnLmZpbGxSZWN0KDAsbkhIZWFkZXJDb2xzLCBlQ2FudmFzVGhpcy5vZmZzZXRXaWR0aCwgZUNhbnZhc1RoaXMub2Zmc2V0SGVpZ2h0KTtcclxuXHJcbiAgICAgIGdyaWQuc2V0T3B0aW9ucyh7XHJcbiAgICAgICAgcm93SGVpZ2h0OiBuSFJvd0dyaWQgLy90aGlzIHdvbid0IHRyaWdnZXIgb25Sb3dzUmV6aXplZCBldmVudCwgd2hpY2ggaXMgYSBERyBidWdcclxuICAgICAgfSk7XHJcblxyXG4gICAgICBub3RpZnlBbGxQaW5uZWRDb2xzUm93c1Jlc2l6ZWQodGhpcywgbkhSb3dHcmlkLCB0cnVlKTtcclxuICAgICAgbm90aWZ5QWxsQ29sc1Jvd3NSZXNpemVkKGdyaWQsIG5IUm93R3JpZCwgdHJ1ZSk7XHJcblxyXG4gICAgICBsZXQgaGVhZGVyID0gbnVsbDtcclxuICAgICAgY29uc3QgYXIgPSBncmlkLmRhcnQubV9hclBpbm5lZENvbHM7XHJcbiAgICAgIGZvcihsZXQgbj0wOyBuPGFyLmxlbmd0aDsgKytuKSB7XHJcbiAgICAgICAgaGVhZGVyID0gYXJbbl07XHJcbiAgICAgICAgZyA9IGhlYWRlci5tX3Jvb3QuZ2V0Q29udGV4dCgnMmQnKTtcclxuICAgICAgICBoZWFkZXIucGFpbnQoZywgZ3JpZCk7XHJcbiAgICAgIH1cclxuXHJcbiAgICAgIHRyeSB7XHJcbiAgICAgICAgY29uc3QgY29sR3JpZDAgPSBncmlkLmNvbHVtbnMuYnlJbmRleCgwKTtcclxuICAgICAgICBpZiAoY29sR3JpZDAgIT09IG51bGwpXHJcbiAgICAgICAgICBjb2xHcmlkMC52aXNpYmxlID0gZmFsc2U7Ly90ZW1wb3JhcnkgYWRkcmVzc2VkIHRoZSBERyBidWdcclxuICAgICAgfVxyXG4gICAgICBjYXRjaChlKSB7XHJcbiAgICAgICAgLy9ERyBidWdcclxuICAgICAgfVxyXG4gICAgICByZXR1cm47XHJcbiAgICB9XHJcblxyXG5cclxuICB9XHJcblxyXG4gIHB1YmxpYyBvbk1vdXNlTGVhdmUoZSA6IE1vdXNlRXZlbnQsIGJPdmVybGFwIDogYm9vbGVhbikgOiB2b2lkIHtcclxuICAgIGlmKERFQlVHKVxyXG4gICAgIGNvbnNvbGUubG9nKCdNb3VzZSBMZWZ0IFBpbm5lZCBDb2x1bW46ICcgKyB0aGlzLmdldEdyaWRDb2x1bW4oKT8ubmFtZSArICcgIG92ZXJsYXA6ICcgKyBiT3ZlcmxhcCk7XHJcblxyXG4gICAgaWYodGhpcy5tX25SZXNpemVSb3dHcmlkTW92aW5nID49IDApIHtcclxuICAgICAgdGhpcy5tX25SZXNpemVSb3dHcmlkTW92aW5nID0gLTE7XHJcbiAgICAgIGRvY3VtZW50LmJvZHkuc3R5bGUuY3Vyc29yID0gXCJhdXRvXCI7XHJcbiAgICB9XHJcblxyXG4gICAgaWYodGhpcy5tX2NlbGxDdXJyZW50ICE9PSBudWxsKSB7XHJcbiAgICAgIGNvbnN0IHJlbmRlcmVyID0gZ2V0UmVuZGVyZXIodGhpcy5tX2NlbGxDdXJyZW50KTtcclxuICAgICAgaWYgKHJlbmRlcmVyIGluc3RhbmNlb2YgR3JpZENlbGxSZW5kZXJlckV4KSB7XHJcbiAgICAgICAgY29uc3QgZU1vdXNlID0gZSBhcyBNb3VzZUV2ZW50O1xyXG4gICAgICAgIHJlbmRlcmVyLm9uTW91c2VMZWF2ZUV4KHRoaXMubV9jZWxsQ3VycmVudCwgZU1vdXNlLCAtMSwgLTEpO1xyXG4gICAgICB9XHJcbiAgICAgIHRoaXMubV9jZWxsQ3VycmVudCA9IG51bGw7XHJcbiAgICB9XHJcblxyXG4gICAgY29uc3QgY29sR3JpZCA9IHRoaXMuZ2V0R3JpZENvbHVtbigpO1xyXG4gICAgaWYoY29sR3JpZCAhPSBudWxsICYmICFiT3ZlcmxhcCkge1xyXG4gICAgICBjb25zdCBlRGl2SGFtYiA9IEdyaWRVdGlscy5nZXRUb29sSWNvbkRpdihjb2xHcmlkLmdyaWQpO1xyXG4gICAgICBlRGl2SGFtYj8uc2V0QXR0cmlidXRlKCdjb2x1bW5fbmFtZScsICcnKTtcclxuICAgICAgLy8gQHRzLWlnbm9yZVxyXG4gICAgICBlRGl2SGFtYj8uc3R5bGUudmlzaWJpbGl0eSA9ICdoaWRkZW4nO1xyXG4gICAgfVxyXG5cclxuXHJcbiAgfVxyXG5cclxuICBwdWJsaWMgb25Nb3VzZURibENsaWNrKGUgOiBNb3VzZUV2ZW50KSA6IHZvaWQge1xyXG4gICAgaWYoREVCVUcpXHJcbiAgICAgY29uc29sZS5sb2coJ01vdXNlIERibCBDbGlja2VkIFBpbm5lZCBDb2x1bW46ICcgKyB0aGlzLmdldEdyaWRDb2x1bW4oKT8ubmFtZSk7XHJcblxyXG4gICAgaWYodGhpcy5tX2NvbEdyaWQgPT09IG51bGwgfHwgdGhpcy5tX3Jvb3QgPT09IG51bGwpXHJcbiAgICAgIHJldHVybjtcclxuXHJcbiAgICBjb25zdCBncmlkID0gdGhpcy5tX2NvbEdyaWQuZ3JpZDtcclxuICAgIGNvbnN0IHZpZXdUYWJsZSA9IGdyaWQ/LnZpZXc7XHJcblxyXG4gICAgaWYgKERHLnRvRGFydChncm9rLnNoZWxsLnYpICE9PSBERy50b0RhcnQodmlld1RhYmxlKSkge1xyXG4gICAgICByZXR1cm47XHJcbiAgICB9XHJcblxyXG4gICAgaWYodGhpcy5tX2NvbEdyaWQ/Lm5hbWUgPT09ICcnKVxyXG4gICAgICByZXR1cm47XHJcblxyXG4gICAgaWYodGhpcy5tX2JTb3J0ZWRBc2NlbmRpbmcgPT0gbnVsbClcclxuICAgICAgdGhpcy5tX2JTb3J0ZWRBc2NlbmRpbmcgPSB0cnVlO1xyXG4gICAgZWxzZSBpZih0aGlzLm1fYlNvcnRlZEFzY2VuZGluZylcclxuICAgICAgdGhpcy5tX2JTb3J0ZWRBc2NlbmRpbmcgPSBmYWxzZTtcclxuICAgIGVsc2UgdGhpcy5tX2JTb3J0ZWRBc2NlbmRpbmcgPSB0cnVlO1xyXG5cclxuICAgIGNvbnN0IG5ISGVhZGVyQ29scyA9IEdyaWRVdGlscy5nZXRHcmlkQ29sdW1uSGVhZGVySGVpZ2h0KGdyaWQpO1xyXG5cclxuICAgIGlmKDAgPD0gZS5vZmZzZXRYICYmIGUub2Zmc2V0WCA8PSB0aGlzLm1fcm9vdC5vZmZzZXRXaWR0aCAmJlxyXG4gICAgICAgIDAgPD0gZS5vZmZzZXRZICYmIGUub2Zmc2V0WSA8PSBuSEhlYWRlckNvbHMpICAgLy9vbiB0aGUgcm93cyBoZWFkZXJcclxuICAgIHtcclxuICAgICAgZ3JpZD8uc29ydChbdGhpcy5tX2NvbEdyaWQ/Lm5hbWVdLCBbdGhpcy5tX2JTb3J0ZWRBc2NlbmRpbmddKTtcclxuICAgIH1cclxuICB9XHJcblxyXG4gIHB1YmxpYyBvbk1vdXNlRG93bihlIDogTW91c2VFdmVudCkgOiB2b2lkIHtcclxuICAgIGlmKERFQlVHKVxyXG4gICAgIGNvbnNvbGUubG9nKCdNb3VzZSBEb3duIFBpbm5lZCBDb2x1bW46ICcgKyB0aGlzLmdldEdyaWRDb2x1bW4oKT8ubmFtZSk7XHJcbi8qXHJcbiAgICBpZihlLnZpZXcgIT0gbnVsbCkge1xyXG4gICAgICBjb25zdCBlZSA9IGRvY3VtZW50LmNyZWF0ZUV2ZW50KCBcIk1vdXNlRXZlbnRcIiApO1xyXG4gICAgICBlZS5pbml0TW91c2VFdmVudChlLnR5cGUsIGUuYnViYmxlcywgZS5jYW5jZWxhYmxlLCBlLnZpZXcsIGUuZGV0YWlsLCBlLnNjcmVlblggKyAxMDAsIGUuc2NyZWVuWSwgZS5jbGllbnRYICsgMTAwLCBlLmNsaWVudFksIGUuY3RybEtleSwgZS5hbHRLZXksIGUuc2hpZnRLZXksIGUubWV0YUtleSwgZS5idXR0b24sIGUucmVsYXRlZFRhcmdldCk7XHJcbiAgICAgIHRoaXMubV9jb2xHcmlkPy5ncmlkLnJvb3QuZGlzcGF0Y2hFdmVudChlZSk7XHJcbiAgICAgIHJldHVybjtcclxuICAgIH1cclxuKi9cclxuXHJcbiAgICBpZih0aGlzLm1fY29sR3JpZCA9PT0gbnVsbClcclxuICAgICAgcmV0dXJuO1xyXG5cclxuICAgIGNvbnN0IGdyaWQgPSB0aGlzLm1fY29sR3JpZD8uZ3JpZDtcclxuICAgIGNvbnN0IHZpZXdUYWJsZSA9IGdyaWQ/LnZpZXc7XHJcbiAgICBpZihERy50b0RhcnQoZ3Jvay5zaGVsbC52KSAhPT0gREcudG9EYXJ0KHZpZXdUYWJsZSkpXHJcbiAgICAgIHJldHVybjtcclxuXHJcbiAgICBpZihlLmJ1dHRvbnMgIT09IDEpXHJcbiAgICAgIHJldHVybjtcclxuXHJcbiAgICBsZXQgZUNhbnZhc1RoaXMgPSB0aGlzLm1fcm9vdDtcclxuICAgIGlmKGVDYW52YXNUaGlzID09PSBudWxsKVxyXG4gICAgICByZXR1cm47XHJcblxyXG4gICAgdGhpcy5tX25SZXNpemVSb3dHcmlkTW92aW5nID0gLTE7XHJcbiAgICBjb25zdCBiQWRkVG9TZWwgOiBib29sZWFuID0gZS5jdHJsS2V5IHx8IGUuc2hpZnRLZXk7XHJcblxyXG4gICAgbGV0IG5Sb3dHcmlkID0gYkFkZFRvU2VsID8gLTEgOiBQaW5uZWRDb2x1bW4uaGl0VGVzdFJvd3MoZUNhbnZhc1RoaXMsIGdyaWQsIGUsIHRydWUsIHVuZGVmaW5lZCk7XHJcbiAgICBpZiAoblJvd0dyaWQgPj0gMCkge1xyXG4gICAgICBjb25zdCBuSFJvd3MgPSBHcmlkVXRpbHMuZ2V0R3JpZFJvd0hlaWdodChncmlkKTtcclxuICAgICAgdGhpcy5tX25SZXNpemVSb3dHcmlkRHJhZ2dpbmcgPSBuUm93R3JpZDtcclxuICAgICAgdGhpcy5tX25ZUmVzaXplRHJhZ2dpbmdBbmNob3IgPSBlLmNsaWVudFk7XHJcbiAgICAgIHRoaXMubV9uSFJlc2l6ZVJvd3NCZWZvcmVEcmFnID0gbkhSb3dzO1xyXG4gICAgfVxyXG4gICAgZWxzZVxyXG4gICAge1xyXG5cclxuICAgICAgblJvd0dyaWQgPSBQaW5uZWRDb2x1bW4uaGl0VGVzdFJvd3MoZUNhbnZhc1RoaXMsIGdyaWQsIGUsIGZhbHNlLCB0aGlzLm1fYXJYWU1vdXNlT25DZWxsRG93bik7XHJcblxyXG4gICAgICB0aGlzLm1fblJvd0dyaWREcmFnZ2luZyA9IG5Sb3dHcmlkO1xyXG4gICAgICB0aGlzLm1fbllEcmFnZ2luZ0FuY2hvciA9IGUuY2xpZW50WTtcclxuXHJcbiAgICAgIGNvbnN0IGNlbGwgPSBncmlkLmNlbGwodGhpcy5tX2NvbEdyaWQubmFtZSwgblJvd0dyaWQpO1xyXG4gICAgICBjb25zdCByZW5kZXJlciA9IGdldFJlbmRlcmVyKGNlbGwpO1xyXG4gICAgICBpZihyZW5kZXJlciBpbnN0YW5jZW9mIEdyaWRDZWxsUmVuZGVyZXJFeCkge1xyXG4gICAgICAgIHJlbmRlcmVyLm9uTW91c2VEb3duRXgoY2VsbCwgZSwgdGhpcy5tX2FyWFlNb3VzZU9uQ2VsbERvd25bMF0sIHRoaXMubV9hclhZTW91c2VPbkNlbGxEb3duWzFdKTtcclxuICAgICAgfVxyXG4gICAgfVxyXG4gIH1cclxuXHJcbiAgcHVibGljIG9uTW91c2VVcChlIDogTW91c2VFdmVudCkgOiB2b2lkIHtcclxuICAgIGlmKERFQlVHKVxyXG4gICAgIGNvbnNvbGUubG9nKCdNb3VzZSBVcCBQaW5uZWQgQ29sdW1uOiAnICsgdGhpcy5nZXRHcmlkQ29sdW1uKCk/Lm5hbWUpO1xyXG4vKlxyXG4gICAgaWYoZS52aWV3ICE9IG51bGwpIHtcclxuICAgICAgY29uc3QgZWUgPSBkb2N1bWVudC5jcmVhdGVFdmVudCggXCJNb3VzZUV2ZW50XCIgKTtcclxuICAgICAgZWUuaW5pdE1vdXNlRXZlbnQoZS50eXBlLCBlLmJ1YmJsZXMsIGUuY2FuY2VsYWJsZSwgZS52aWV3LCBlLmRldGFpbCwgZS5zY3JlZW5YICsgMTAwLCBlLnNjcmVlblksIGUuY2xpZW50WCArIDEwMCwgZS5jbGllbnRZLCBlLmN0cmxLZXksIGUuYWx0S2V5LCBlLnNoaWZ0S2V5LCBlLm1ldGFLZXksIGUuYnV0dG9uLCBlLnJlbGF0ZWRUYXJnZXQpO1xyXG4gICAgICB0aGlzLm1fY29sR3JpZD8uZ3JpZC5yb290LmRpc3BhdGNoRXZlbnQoZWUpO1xyXG4gICAgICByZXR1cm47XHJcbiAgICB9XHJcbiovXHJcblxyXG4gICAgaWYodGhpcy5tX2NvbEdyaWQgPT09IG51bGwgfHwgdGhpcy5tX3Jvb3QgPT0gbnVsbClcclxuICAgICAgcmV0dXJuO1xyXG5cclxuICAgIGNvbnN0IGdyaWQgPSB0aGlzLm1fY29sR3JpZD8uZ3JpZDtcclxuICAgIGNvbnN0IHZpZXdUYWJsZSA9IGdyaWQ/LnZpZXc7XHJcblxyXG4gICAgaWYoREcudG9EYXJ0KGdyb2suc2hlbGwudikgIT09IERHLnRvRGFydCh2aWV3VGFibGUpKSB7XHJcbiAgICAgIHJldHVybjtcclxuICAgIH1cclxuLypcclxuICAgaWYoZS5idXR0b24gPT09IDIpIHtcclxuXHJcbiAgICAgY29uc3QgZURpdlBPcHVwIDogSFRNTEVsZW1lbnQgfCBudWxsID0gR3JpZFV0aWxzLmdldEdyaWREYXJ0UG9wdXBNZW51KCk7XHJcbiAgICAgZURpdlBPcHVwPy5zZXRBdHRyaWJ1dGUoJ2NvbHVtbl9uYW1lJywgdGhpcy5tX2NvbEdyaWQubmFtZSk7XHJcbiAgICAgbGV0IGQgPSAwO1xyXG4gICAgIHJldHVybjtcclxuICAgfSovXHJcblxyXG5cclxuICAgIGlmKHRoaXMubV9uUmVzaXplUm93R3JpZERyYWdnaW5nID49IDApIHtcclxuICAgICAgY29uc3QgbkhSb3cgPSBHcmlkVXRpbHMuZ2V0R3JpZFJvd0hlaWdodChncmlkKTtcclxuICAgICAgbm90aWZ5QWxsUGlubmVkQ29sc1Jvd3NSZXNpemVkKHRoaXMsIG5IUm93LCBmYWxzZSk7XHJcbiAgICAgIG5vdGlmeUFsbENvbHNSb3dzUmVzaXplZChncmlkLCBuSFJvdywgZmFsc2UpO1xyXG4gICAgfVxyXG5cclxuICAgIHRoaXMubV9uSFJlc2l6ZVJvd3NCZWZvcmVEcmFnID0gLTE7XHJcbiAgICB0aGlzLm1fblJlc2l6ZVJvd0dyaWREcmFnZ2luZyA9IC0xO1xyXG4gICAgdGhpcy5tX25ZUmVzaXplRHJhZ2dpbmdBbmNob3IgPSAtMTtcclxuICAgIHRoaXMubV9uUmVzaXplUm93R3JpZE1vdmluZyA9IC0xO1xyXG5cclxuICAgIGRvY3VtZW50LmJvZHkuc3R5bGUuY3Vyc29yID0gXCJhdXRvXCI7XHJcblxyXG4gICAgaWYodGhpcy5tX25Sb3dHcmlkRHJhZ2dpbmcgPj0gMCkge1xyXG4gICAgICBjb25zdCBkZnJhbWUgPSBncmlkLmRhdGFGcmFtZTtcclxuICAgICAgY29uc3QgYkFkZFRvU2VsID0gZS5jdHJsS2V5O1xyXG4gICAgICBjb25zdCBiUmFuZ2VTZWwgPSBlLnNoaWZ0S2V5O1xyXG5cclxuICAgICAgY29uc3QgblJvd0dyaWQgPSBQaW5uZWRDb2x1bW4uaGl0VGVzdFJvd3ModGhpcy5tX3Jvb3QsIGdyaWQsIGUsIGZhbHNlLCB0aGlzLm1fYXJYWU1vdXNlT25DZWxsVXApO1xyXG4gICAgICBpZighYkFkZFRvU2VsICYmICFiUmFuZ2VTZWwgJiYgblJvd0dyaWQgPT09IHRoaXMubV9uUm93R3JpZERyYWdnaW5nKSB7IC8vY2xpY2sgb24gdGhlIHNhbWUgcm93IHdoaWNoIHdpbGwgYmVjb21lIGFjdGl2ZVxyXG5cclxuICAgICAgICBsZXQgY2VsbFJIID0gbnVsbDtcclxuICAgICAgICB0cnkge1xyXG4gICAgICAgICAgY2VsbFJIID0gZ3JpZC5jZWxsKFwiXCIsIG5Sb3dHcmlkKTtcclxuICAgICAgICB9XHJcbiAgICAgICAgY2F0Y2goZSkge1xyXG4gICAgICAgICAgbGV0IGNvbEcgPSBudWxsO1xyXG4gICAgICAgICAgY29uc3QgbHN0Q29scyA9IGdyaWQuY29sdW1ucztcclxuICAgICAgICAgIGZvcihsZXQgbkM9MTsgbkM8bHN0Q29scy5sZW5ndGg7ICsrbkMpIHtcclxuICAgICAgICAgICAgY29sRyA9IGxzdENvbHMuYnlJbmRleChuQyk7XHJcbiAgICAgICAgICAgIGNlbGxSSCA9IGNvbEcgPT09IG51bGwgPyBudWxsIDogZ3JpZC5jZWxsKGNvbEcubmFtZSwgblJvd0dyaWQpO1xyXG4gICAgICAgICAgICBpZihjZWxsUkggIT09IG51bGwpXHJcbiAgICAgICAgICAgICAgYnJlYWs7XHJcbiAgICAgICAgICB9XHJcbiAgICAgICAgfVxyXG4gICAgICAgIGlmKGNlbGxSSCAhPT0gbnVsbCkge1xyXG4gICAgICAgICAgY29uc3QgblJvd1RhYmxlIDogYW55ID0gY2VsbFJILnRhYmxlUm93SW5kZXg7XHJcbiAgICAgICAgICBpZihuUm93VGFibGUgIT09IG51bGwpXHJcbiAgICAgICAgICAgIGRmcmFtZS5jdXJyZW50Um93ID0gblJvd1RhYmxlO1xyXG4gICAgICAgIH1cclxuICAgICAgfVxyXG4gICAgICBlbHNlXHJcbiAgICAgIHtcclxuICAgICAgICBjb25zdCBiaXRzZXRTZWwgPSBkZnJhbWUuc2VsZWN0aW9uO1xyXG5cclxuICAgICAgICBpZighYkFkZFRvU2VsIHx8IGJSYW5nZVNlbClcclxuICAgICAgICAgIGJpdHNldFNlbC5zZXRBbGwoZmFsc2UsIHRydWUpO1xyXG5cclxuICAgICAgICBsZXQgblJvd01pbiA9IHRoaXMubV9uUm93R3JpZERyYWdnaW5nIDwgblJvd0dyaWQgPyB0aGlzLm1fblJvd0dyaWREcmFnZ2luZyA6IG5Sb3dHcmlkO1xyXG4gICAgICAgIGxldCBuUm93TWF4ID0gdGhpcy5tX25Sb3dHcmlkRHJhZ2dpbmcgPiBuUm93R3JpZCA/IHRoaXMubV9uUm93R3JpZERyYWdnaW5nIDogblJvd0dyaWQ7XHJcblxyXG4gICAgICAgIGlmKGJSYW5nZVNlbCkge1xyXG4gICAgICAgICAgbGV0IG5Sb3dHcmlkQWN0aXZlID0gR3JpZFV0aWxzLmdldEFjdGl2ZUdyaWRSb3coZ3JpZCk7XHJcbiAgICAgICAgICBpZihuUm93R3JpZEFjdGl2ZSA9PT0gbnVsbClcclxuICAgICAgICAgICAgblJvd0dyaWRBY3RpdmUgPSAwO1xyXG5cclxuICAgICAgICAgIG5Sb3dNaW4gPSBuUm93R3JpZEFjdGl2ZSA8IG5Sb3dHcmlkID8gblJvd0dyaWRBY3RpdmUgOiBuUm93R3JpZDtcclxuICAgICAgICAgIG5Sb3dNYXggPSBuUm93R3JpZEFjdGl2ZSA+IG5Sb3dHcmlkID8gblJvd0dyaWRBY3RpdmUgOiBuUm93R3JpZDtcclxuICAgICAgICB9XHJcblxyXG5cclxuICAgICAgICBsZXQgY2VsbFJIID0gbnVsbDtcclxuICAgICAgICBsZXQgblJvd1RhYmxlID0gLTE7XHJcbiAgICAgICAgZm9yKGxldCBuUm93PW5Sb3dNaW47IG5Sb3c8PW5Sb3dNYXg7ICsrblJvdykge1xyXG5cclxuICAgICAgICAgIHRyeSB7XHJcbiAgICAgICAgICAgIGNlbGxSSCA9IGdyaWQuY2VsbChcIlwiLCBuUm93KTtcclxuICAgICAgICAgIH1cclxuICAgICAgICAgIGNhdGNoKGUpIHtcclxuICAgICAgICAgICAgbGV0IGNvbEcgPSBudWxsO1xyXG4gICAgICAgICAgICBjb25zdCBsc3RDb2xzID0gZ3JpZC5jb2x1bW5zO1xyXG4gICAgICAgICAgICBmb3IobGV0IG5DPTE7IG5DPGxzdENvbHMubGVuZ3RoOyArK25DKSB7XHJcbiAgICAgICAgICAgICAgY29sRyA9IGxzdENvbHMuYnlJbmRleChuQyk7XHJcbiAgICAgICAgICAgICAgY2VsbFJIID0gY29sRyA9PT0gbnVsbCA/IG51bGwgOiBncmlkLmNlbGwoY29sRy5uYW1lLCBuUm93R3JpZCk7XHJcbiAgICAgICAgICAgICAgaWYoY2VsbFJIICE9PSBudWxsKVxyXG4gICAgICAgICAgICAgICAgYnJlYWs7XHJcbiAgICAgICAgICAgIH1cclxuICAgICAgICAgIH1cclxuXHJcbiAgICAgICAgICBpZihjZWxsUkggIT09IG51bGwgJiYgY2VsbFJILnRhYmxlUm93SW5kZXggIT09IG51bGwpIHtcclxuICAgICAgICAgICAgblJvd1RhYmxlID0gY2VsbFJILnRhYmxlUm93SW5kZXg7XHJcbiAgICAgICAgICAgIGJpdHNldFNlbC5zZXQoblJvd1RhYmxlLCB0cnVlLCB0cnVlKTtcclxuICAgICAgICAgIH1cclxuICAgICAgICB9XHJcbiAgICAgIH1cclxuXHJcbiAgICAgIGNvbnN0IGNlbGwgPSBncmlkLmNlbGwodGhpcy5tX2NvbEdyaWQubmFtZSwgblJvd0dyaWQpO1xyXG4gICAgICBjb25zdCByZW5kZXJlciA9IGdldFJlbmRlcmVyKGNlbGwpO1xyXG4gICAgICBpZihyZW5kZXJlciBpbnN0YW5jZW9mIEdyaWRDZWxsUmVuZGVyZXJFeCkge1xyXG4gICAgICAgIHJlbmRlcmVyLm9uTW91c2VVcEV4KGNlbGwsIGUsIHRoaXMubV9hclhZTW91c2VPbkNlbGxVcFswXSwgdGhpcy5tX2FyWFlNb3VzZU9uQ2VsbFVwWzFdKTtcclxuICAgICAgfVxyXG5cclxuICAgICAgaWYodGhpcy5tX2FyWFlNb3VzZU9uQ2VsbFVwWzBdID09PSB0aGlzLm1fYXJYWU1vdXNlT25DZWxsRG93blswXSAmJiB0aGlzLm1fYXJYWU1vdXNlT25DZWxsRG93blsxXSA9PT0gdGhpcy5tX2FyWFlNb3VzZU9uQ2VsbFVwWzFdKSB7XHJcbiAgICAgICAgaWYocmVuZGVyZXIgaW5zdGFuY2VvZiBHcmlkQ2VsbFJlbmRlcmVyRXgpIHtcclxuICAgICAgICAgIHJlbmRlcmVyLm9uQ2xpY2tFeChjZWxsLCBlLCB0aGlzLm1fYXJYWU1vdXNlT25DZWxsVXBbMF0sIHRoaXMubV9hclhZTW91c2VPbkNlbGxVcFsxXSk7XHJcbiAgICAgICAgfVxyXG4gICAgICB9XHJcblxyXG4gICAgICB0aGlzLm1fblJvd0dyaWREcmFnZ2luZyA9IC0xO1xyXG4gICAgICB0aGlzLm1fbllEcmFnZ2luZ0FuY2hvciA9IC0xO1xyXG4gICAgICB0aGlzLm1fYXJYWU1vdXNlT25DZWxsRG93blswXSA9IC0yO1xyXG4gICAgICB0aGlzLm1fYXJYWU1vdXNlT25DZWxsRG93blsxXSA9IC0yO1xyXG4gICAgICB0aGlzLm1fYXJYWU1vdXNlT25DZWxsVXBbMF0gPSAtMTtcclxuICAgICAgdGhpcy5tX2FyWFlNb3VzZU9uQ2VsbFVwWzFdID0gLTE7XHJcbiAgICB9XHJcbiAgfVxyXG5cclxuICBwdWJsaWMgb25Db250ZXh0TWVudShlIDogTW91c2VFdmVudCkgOiB2b2lkIHtcclxuICAgaWYoREVCVUcpXHJcbiAgICBjb25zb2xlLmxvZygnQ29udGV4dCBtZW51IFBpbm5lZCBDb2x1bW46ICcgKyB0aGlzLmdldEdyaWRDb2x1bW4oKT8ubmFtZSk7XHJcbiAgfVxyXG5cclxuICBwdWJsaWMgb25Nb3VzZVdoZWVsKGUgOiBXaGVlbEV2ZW50KSA6IHZvaWQge1xyXG5cclxuICAgIGlmKHRoaXMubV9jb2xHcmlkID09PSBudWxsIHx8IHRoaXMubV9yb290ID09IG51bGwpXHJcbiAgICAgIHJldHVybjtcclxuXHJcbiAgICBjb25zdCBncmlkID0gdGhpcy5tX2NvbEdyaWQ/LmdyaWQ7XHJcbiAgICBjb25zdCB2aWV3VGFibGUgPSBncmlkPy52aWV3O1xyXG5cclxuICAgIGlmIChERy50b0RhcnQoZ3Jvay5zaGVsbC52KSAhPT0gREcudG9EYXJ0KHZpZXdUYWJsZSkpIHtcclxuICAgICAgcmV0dXJuO1xyXG4gICAgfVxyXG5cclxuICAgIGlmKGUuZGVsdGFYICE9PSAwIHx8IGUuZGVsdGFaICE9PSAwKSB7XHJcbiAgICAgIHJldHVybjtcclxuICAgIH1cclxuXHJcbiAgICBzZXRUaW1lb3V0KCgpID0+e1xyXG4gICAgICBjb25zdCBlZSA9IG5ldyBXaGVlbEV2ZW50KGUudHlwZSwgZSk7XHJcbiAgICAgIHRyeXtncmlkLm92ZXJsYXkuZGlzcGF0Y2hFdmVudChlZSk7fVxyXG4gICAgICBjYXRjaChleCkge1xyXG4gICAgICAgIC8vY29uc29sZS5lcnJvcihleC5tZXNzYWdlKTtcclxuICAgICAgfVxyXG4gICAgfSwgMSk7XHJcblxyXG5cclxuICAgIGlmKHRydWUpXHJcbiAgICAgIHJldHVybjtcclxuICAgIC8vZS5jbGllbnRYID0gNTtcclxuXHJcblxyXG4gICAgaWYodGhpcy5tX25XaGVlbENvdW50ID09PSAxKSB7XHJcbiAgICAgIC8vc2Nyb2xsICtcclxuICAgICAgY29uc3QgblJvd0NvdW50ID0gR3JpZFV0aWxzLmdldEdyaWRWaXNpYmxlUm93Q291bnQoZ3JpZCk7XHJcbiAgICAgIGNvbnN0IHNjcm9sbFkgPSBncmlkLnZlcnRTY3JvbGw7XHJcbiAgICAgIGlmKG5Sb3dDb3VudCAtMSA+IHNjcm9sbFkubWF4KSB7XHJcbiAgICAgICAgc2Nyb2xsWS5zZXRWYWx1ZXMoc2Nyb2xsWS5taW5SYW5nZSwgc2Nyb2xsWS5tYXhSYW5nZSwgc2Nyb2xsWS5taW4gKyAxLCBzY3JvbGxZLm1heCArIDEpO1xyXG4gICAgICB9XHJcbiAgICAgIHRoaXMubV9uV2hlZWxDb3VudCA9IDA7XHJcbiAgICB9XHJcbiAgICBlbHNlIGlmKHRoaXMubV9uV2hlZWxDb3VudCA9PT0gLTEpXHJcbiAgICB7XHJcbiAgICAgIC8vc2Nyb2xsIC1cclxuICAgICAgY29uc3Qgc2Nyb2xsWSA9IGdyaWQudmVydFNjcm9sbDtcclxuICAgICAgaWYoc2Nyb2xsWS5taW4gPj0xKSB7XHJcbiAgICAgICAgc2Nyb2xsWS5zZXRWYWx1ZXMoc2Nyb2xsWS5taW5SYW5nZSwgc2Nyb2xsWS5tYXhSYW5nZSwgc2Nyb2xsWS5taW4gLSAxLCBzY3JvbGxZLm1heCAtIDEpO1xyXG4gICAgICB9XHJcbiAgICAgIHRoaXMubV9uV2hlZWxDb3VudCA9IDA7XHJcbiAgICB9XHJcbiAgICBlbHNlIHtcclxuICAgICAgdGhpcy5tX25XaGVlbENvdW50ID0gZS5kZWx0YVkgPiAwID8gMSA6IC0xO1xyXG4gICAgfVxyXG4gIH1cclxuXHJcblxyXG4gIHByaXZhdGUgcGFpbnQoZyA6IENhbnZhc1JlbmRlcmluZ0NvbnRleHQyRCB8IG51bGwsIGdyaWQgOiBERy5HcmlkKSA6IHZvaWQge1xyXG4gICAgLy9jb25zdCBuV0RpdiA9IGVudHJ5LmNvbnRlbnRCb3hTaXplID8gZW50cnkuY29udGVudEJveFNpemVbMF0uaW5saW5lU2l6ZSA6IGVudHJ5LmNvbnRlbnRSZWN0LndpZHRoO1xyXG5cclxuICAgIGlmKGcgPT09IG51bGwpIHtcclxuICAgICAgcmV0dXJuO1xyXG4gICAgfVxyXG5cclxuICAgIGlmKHRoaXMubV9yb290ID09PSBudWxsKSB7XHJcbiAgICAgIHRocm93IG5ldyBFcnJvcignUm9vdCBjYW5ub3QgYmUgbnVsbC4nKTtcclxuICAgIH1cclxuXHJcbiAgICBpZih0aGlzLm1fY29sR3JpZCA9PT0gbnVsbCkge1xyXG4gICAgICB0aHJvdyBuZXcgRXJyb3IoJ0NvbHVtbiBncmlkIGNhbm5vdCBiZSBudWxsLicpO1xyXG4gICAgfVxyXG4gICAgY29uc3QgZGZyYW1lID0gZ3JpZC5kYXRhRnJhbWU7XHJcbiAgICBjb25zdCBuVyA9IHRoaXMubV9yb290Lm9mZnNldFdpZHRoO1xyXG4gICAgY29uc3QgbkggPSB0aGlzLm1fcm9vdC5vZmZzZXRIZWlnaHQ7XHJcblxyXG4gICAgZy5maWxsU3R5bGUgPSBcIndoaXRlXCI7XHJcbiAgICBnLmZpbGxSZWN0KDAsMCwgblcqd2luZG93LmRldmljZVBpeGVsUmF0aW8sIG5IKndpbmRvdy5kZXZpY2VQaXhlbFJhdGlvKTtcclxuXHJcbiAgICBpZih0aGlzLm1fY29sR3JpZC5uYW1lID09PSBudWxsKVxyXG4gICAgICByZXR1cm47XHJcblxyXG4gICAgY29uc3QgYml0c2V0RmlsdGVyID0gZGZyYW1lLmZpbHRlcjtcclxuICAgIGlmKGJpdHNldEZpbHRlci5mYWxzZUNvdW50ID09PSBkZnJhbWUucm93Q291bnQpXHJcbiAgICAgIHJldHVybjsgLy9ldmVyeXRoaW5nIGlzIGZpbHRlcmVkXHJcblxyXG4gICAgLy9jb2x1bW4gSGVhZGVyXHJcbiAgICBjb25zdCBvcHRpb25zIDogYW55ID0gZ3JpZC5nZXRPcHRpb25zKHRydWUpO1xyXG5cclxuICAgIGNvbnN0IGZvbnRDZWxsRGVmYXVsdCA9IG9wdGlvbnMubG9vay5kZWZhdWx0Q2VsbEZvbnQ7XHJcblxyXG4gICAgbGV0IGZvbnQgPSBvcHRpb25zLmxvb2suY29sSGVhZGVyRm9udCA9PSBudWxsIHx8IG9wdGlvbnMubG9vay5jb2xIZWFkZXJGb250ID09PSB1bmRlZmluZWQgPyBcImJvbGQgMTRweCBWb2x0YSBUZXh0LCBBcmlhbFwiIDogb3B0aW9ucy5sb29rLmNvbEhlYWRlckZvbnQ7XHJcbiAgICBsZXQgZm9udFNjYWxlZCA9IEdyaWRVdGlscy5zY2FsZUZvbnQoZm9udCwgd2luZG93LmRldmljZVBpeGVsUmF0aW8pO1xyXG4gICAgZy5mb250ID0gZm9udFNjYWxlZDtcclxuXHJcbiAgICBsZXQgc3RyID0gVGV4dFV0aWxzLnRyaW1UZXh0KHRoaXMubV9jb2xHcmlkLm5hbWUsIGcsIG5XKTtcclxuXHJcbiAgICBjb25zdCB0bSA9IGcubWVhc3VyZVRleHQoc3RyKTtcclxuICAgIGNvbnN0IG5XTGFiZWwgPSB0bS53aWR0aDtcclxuXHJcbiAgICBjb25zdCBuQXNjZW50ID0gTWF0aC5hYnModG0uYWN0dWFsQm91bmRpbmdCb3hBc2NlbnQpO1xyXG4gICAgY29uc3QgbkRlc2NlbnQgPSB0bS5hY3R1YWxCb3VuZGluZ0JveERlc2NlbnQ7XHJcbiAgICBjb25zdCBuSEZvbnQgPSAgbkFzY2VudCArIG5EZXNjZW50Oy8vICsgMipuWUluc2V0O1xyXG5cclxuICAgIC8vbGV0IGNlbGxDSCA9IGdyaWQuY2VsbCh0aGlzLm1fY29sR3JpZC5uYW1lLCAtMSk7XHJcbiAgICAvL2xldCByZW5kZXJlciA9IGNlbGxDSC5yZW5kZXJlcjtcclxuXHJcbiAgICBsZXQgblggPSAwO1xyXG4gICAgbGV0IG5ZID0gMDtcclxuICAgIGNvbnN0IG5IQ0ggPSBHcmlkVXRpbHMuZ2V0R3JpZENvbHVtbkhlYWRlckhlaWdodChncmlkKSp3aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbztcclxuICAgIGcudGV4dEFsaWduID0gJ3N0YXJ0JztcclxuICAgIGcuZmlsbFN0eWxlID0gXCJCbGFja1wiO1xyXG4gICAgbGV0IG5ZT2Zmc2V0ID0gTWF0aC5mbG9vcigobkhDSCAtIG5IRm9udCkvMik7XHJcbiAgICBjb25zdCBuWFggPSBuWCArICgoblcqd2luZG93LmRldmljZVBpeGVsUmF0aW8gLSBuV0xhYmVsKSA+PiAxKTtcclxuICAgIGxldCBuWVkgPSAoblkgKyBuSENIIC0gTWF0aC5jZWlsKDMqd2luZG93LmRldmljZVBpeGVsUmF0aW8pKTsvLy0yKndpbmRvdy5kZXZpY2VQaXhlbFJhdGlvKTtcclxuICAgIC8vb25zb2xlLmxvZyhcIm5YWCBcIiArIG5YWCArIFwiIG5ZWSA9IFwiICsgbllZICsgXCIgQ0hIIFwiICsgbkhDSCk7XHJcbiAgICBnLmZpbGxUZXh0KHN0ciwgblhYLCBuWVkpO1xyXG5cclxuICAgIC8vaWYob3B0aW9ucy5sb29rLnNob3dSb3dHcmlkbGluZXMpIHtcclxuXHJcblxyXG4gICAgLy99XHJcblxyXG5cclxuXHJcbiAgICAvL1JlZ3VsYXIgY2VsbHNcclxuICAgIGNvbnN0IG5Sb3dDdXJyZW50ID0gIGRmcmFtZS5jdXJyZW50Um93LmlkeDtcclxuICAgIGNvbnN0IGJpdHNldFNlbCA9IGRmcmFtZS5zZWxlY3Rpb247XHJcblxyXG4gICAgY29uc3QgYXJSb3dzTWluTWF4ID0gWy0xLC0xXTtcclxuICAgIEdyaWRVdGlscy5maWxsVmlzaWJsZVZpZXdwb3J0Um93cyhhclJvd3NNaW5NYXgsIGdyaWQpO1xyXG4gICAgY29uc3QgblJvd01pbiA9IGFyUm93c01pbk1heFswXTtcclxuICAgIGNvbnN0IG5Sb3dNYXggPSBhclJvd3NNaW5NYXhbMV07XHJcblxyXG4gICAgLy9jb25zb2xlLmxvZyhuUm93TWluICsgXCIgXCIgKyBuUm93TWF4KTtcclxuICAgIGNvbnN0IG5IUm93ID0gR3JpZFV0aWxzLmdldEdyaWRSb3dIZWlnaHQoZ3JpZCk7XHJcbiAgICBuWU9mZnNldCA9IG5IQ0g7XHJcbiAgICBjb25zdCBuSFJvd0dyaWQgPSBuSFJvdyp3aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbztcclxuICAgIGxldCBjZWxsUkggPSBudWxsO1xyXG5cclxuICAgIGxldCBuV1cgPSBuVyp3aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbztcclxuICAgIC8vY29uc3QgbkhIID0gbkhSb3dHcmlkO1xyXG5cclxuICAgIGNvbnN0IGFyVGFibGVSb3dzID0gbmV3IEFycmF5KG5Sb3dNYXggLSBuUm93TWluICsxKTtcclxuICAgIGxldCBuUm93VGFibGUgPSAtMTtcclxuICAgIGxldCBiU2VsID0gZmFsc2U7XHJcbiAgICBmb3IobGV0IG5SRz1uUm93TWluOyBuUkc8PW5Sb3dNYXg7ICsrblJHKSB7XHJcbiAgICAgIHRyeSB7XHJcbiAgICAgICAgY2VsbFJIID0gZ3JpZC5jZWxsKHRoaXMubV9jb2xHcmlkLm5hbWUsIG5SRyk7XHJcbiAgICAgIH0gY2F0Y2ggKGUpIC8vdG8gYWRkcmVzcyBERyBidWcgd2hlbiBldmVyeXRoaW5nIGlzIGZpbHRlcmVkXHJcbiAgICAgIHtcclxuICAgICAgICBjb250aW51ZTtcclxuICAgICAgfVxyXG5cclxuICAgICAgaWYgKGNlbGxSSC50YWJsZVJvd0luZGV4ID09PSB1bmRlZmluZWQpLy9ERyBidWdcclxuICAgICAgICBjb250aW51ZTtcclxuXHJcbiAgICAgIG5Sb3dUYWJsZSA9IGNlbGxSSC50YWJsZVJvd0luZGV4ID09PSBudWxsID8gLTEgOiBjZWxsUkgudGFibGVSb3dJbmRleDtcclxuICAgICAgYXJUYWJsZVJvd3NbblJHIC0gblJvd01pbl0gPSBuUm93VGFibGU7XHJcblxyXG4gICAgICBuWVkgPSBuWU9mZnNldCArIChuUkcgLSBuUm93TWluKSAqIG5IUm93R3JpZDtcclxuXHJcbiAgICAgIGxldCByZW5kZXJlcjogYW55ID0gR3JpZFV0aWxzLmdldEdyaWRDb2x1bW5SZW5kZXJlcihjZWxsUkguZ3JpZENvbHVtbik7XHJcbiAgICAgIGlmIChyZW5kZXJlciA9PT0gbnVsbCkge1xyXG4gICAgICAgIHRyeSB7XHJcbiAgICAgICAgICByZW5kZXJlciA9IGNlbGxSSC5yZW5kZXJlcjtcclxuICAgICAgICB9IGNhdGNoIChlKSB7XHJcbiAgICAgICAgICBjb25zb2xlLmVycm9yKFwiQ291bGQgbm90IG9idGFpbiByZW5kZXJlciBmb3IgREcgY2VsbC4gREcgYnVnIFwiICsgdGhpcy5tX2NvbEdyaWQubmFtZSArIFwiIHJvdyBcIiArIG5SRyk7XHJcbiAgICAgICAgICBjb250aW51ZTtcclxuICAgICAgICB9XHJcbiAgICAgIH1cclxuXHJcbiAgICAgIGlmIChyZW5kZXJlciA9PT0gbnVsbCB8fCByZW5kZXJlciA9PT0gdW5kZWZpbmVkKSB7XHJcbiAgICAgICAgY29uc29sZS5lcnJvcihcIkNvdWxkbid0IGZpbmQgcmVuZGVyZXIgZm9yIHBpbm5lZCBjb2x1bW4gXCIgKyB0aGlzLm1fY29sR3JpZC5uYW1lICsgXCIgcm93IFwiICsgblJHKTtcclxuICAgICAgICBjb250aW51ZTtcclxuICAgICAgfVxyXG5cclxuICAgICAgLy9sZXQgbllZID0gblk7Ly8qd2luZG93LmRldmljZVBpeGVsUmF0aW87XHJcblxyXG5cclxuICAgICAgZm9udCA9IGNlbGxSSC5zdHlsZS5mb250O1xyXG4gICAgICBmb250U2NhbGVkID0gR3JpZFV0aWxzLnNjYWxlRm9udChmb250LCB3aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbyk7XHJcbiAgICAgIGlmIChmb250U2NhbGVkICE9PSBudWxsKSB7XHJcbiAgICAgICAgY2VsbFJILnN0eWxlLmZvbnQgPSBmb250U2NhbGVkO1xyXG4gICAgICB9XHJcblxyXG4gICAgICBpZiAoblcgPiAwICYmIG5IUm93R3JpZCA+IDApIHsgLy90byBhZGRyZXNzIGEgYnVnIGNhdXNlZCBlaXRoZXIgREcgb3IgY2xpZW50IGFwcFxyXG4gICAgICAgIHRyeSB7XHJcbiAgICAgICAgICBpZiAocmVuZGVyZXIubmFtZSA9PT0gJ01vbGVjdWxlJykge1xyXG4gICAgICAgICAgICByZW5kZXJlci5yZW5kZXIoZywgMCwgbllZL3dpbmRvdy5kZXZpY2VQaXhlbFJhdGlvLCBuV1cvd2luZG93LmRldmljZVBpeGVsUmF0aW8sIG5IUm93R3JpZC93aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbywgY2VsbFJILCBjZWxsUkguc3R5bGUpO1xyXG4gICAgICAgICAgfVxyXG4gICAgICAgICAgZWxzZSByZW5kZXJlci5yZW5kZXIoZywgMCwgbllZLCBuV1csIG5IUm93R3JpZCwgY2VsbFJILCBjZWxsUkguc3R5bGUpO1xyXG5cclxuICAgICAgICB9IGNhdGNoIChlKSB7XHJcbiAgICAgICAgICBjb25zb2xlLmVycm9yKFwiQ291bGQgbm90IHBhaW50IGNlbGwgZm9yIHBpbm5lZCBjb2x1bW4gXCIgKyB0aGlzLm1fY29sR3JpZC5uYW1lICsgXCIgcm93IFwiICsgblJHKTtcclxuICAgICAgICAgIGNvbnRpbnVlO1xyXG4gICAgICAgICAgLy90aHJvdyBlO1xyXG4gICAgICAgIH1cclxuICAgICAgfVxyXG4gICAgfVxyXG5cclxuXHJcbiAgICAvL1BhaW50IEdyaWRcclxuICAgIGcuc3Ryb2tlU3R5bGUgPSBcIkdhaW5zYm9yb1wiO1xyXG4gICAgZy5iZWdpblBhdGgoKTtcclxuICAgIGcubW92ZVRvKDAsIG5ZKndpbmRvdy5kZXZpY2VQaXhlbFJhdGlvKTtcclxuICAgIGcubGluZVRvKDAsIChuWSArIG5IQ0gtMSp3aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbykpO1xyXG4gICAgZy5zdHJva2UoKTtcclxuXHJcbiAgICBnLmJlZ2luUGF0aCgpO1xyXG4gICAgZy5tb3ZlVG8oMCwgbllPZmZzZXQgKyAxKTtcclxuICAgIGcubGluZVRvKG5XVywgbllPZmZzZXQgKyAxKTtcclxuICAgIGcuc3Ryb2tlKCk7XHJcblxyXG4gICAgZm9yKGxldCBuUkc9blJvd01pbjsgblJHPD1uUm93TWF4OyArK25SRylcclxuICAgIHtcclxuICAgICAgbllZID0gbllPZmZzZXQgKyAoblJHIC0gblJvd01pbikgKiBuSFJvd0dyaWQ7XHJcblxyXG4gICAgICAvL2lmKG9wdGlvbnMubG9vay5zaG93Um93R3JpZGxpbmVzKSB7XHJcblxyXG4gICAgICAgIGcuYmVnaW5QYXRoKCk7XHJcbiAgICAgICAgZy5tb3ZlVG8oMCwgbllZICsgbkhSb3dHcmlkKzEpO1xyXG4gICAgICAgIGcubGluZVRvKG5XVywgbllZICsgbkhSb3dHcmlkKzEpO1xyXG4gICAgICAgIGcuc3Ryb2tlKCk7XHJcblxyXG4gICAgICAgIGcuYmVnaW5QYXRoKCk7XHJcbiAgICAgICAgZy5tb3ZlVG8oMCwgbllZKTtcclxuICAgICAgICBnLmxpbmVUbygwLCBuWVkgKyBuSFJvd0dyaWQrMSk7XHJcbiAgICAgICAgZy5zdHJva2UoKTtcclxuICAgICAgLy99XHJcbiAgICAgIG5Sb3dUYWJsZSA9IGFyVGFibGVSb3dzW25SRyAtIG5Sb3dNaW5dO1xyXG4gICAgICB0cnl7YlNlbCA9IG5Sb3dUYWJsZSA9PT0gdW5kZWZpbmVkIHx8IG5Sb3dUYWJsZSA8IDAgPyBmYWxzZSA6IGJpdHNldFNlbC5nZXQoblJvd1RhYmxlKTt9XHJcbiAgICAgIGNhdGNoIChlKXtcclxuICAgICAgICBjb25zb2xlLmVycm9yKCdQYWludEVycm9yOiByb3dfbWluOiAnICsgblJvd01pbiArICcgcm93X21heDogJyArIG5Sb3dNYXggKyAnIG5SICcgKyBuUkcgKyAnICcgKyBuUm93VGFibGUpO1xyXG4gICAgICAgIHRocm93IGU7XHJcbiAgICAgIH1cclxuICAgICAgaWYoYlNlbClcclxuICAgICAge1xyXG4gICAgICAgIGcuZ2xvYmFsQWxwaGEgPSAwLjI7XHJcbiAgICAgICAgZy5maWxsU3R5bGUgPSBQaW5uZWRDb2x1bW4uU0VMRUNUSU9OX0NPTE9SO1xyXG4gICAgICAgIGcuZmlsbFJlY3QoMCwgbllZLCBuV1csIG5IUm93R3JpZCk7XHJcbiAgICAgICAgZy5nbG9iYWxBbHBoYSA9IDE7XHJcbiAgICAgIH1cclxuXHJcbiAgICAgIGlmKG5Sb3dDdXJyZW50ID09PSBuUm93VGFibGUpXHJcbiAgICAgIHtcclxuICAgICAgICBnLmdsb2JhbEFscGhhID0gMC4yO1xyXG4gICAgICAgIGcuZmlsbFN0eWxlID0gUGlubmVkQ29sdW1uLkFDVElWRV9DRUxMX0NPTE9SO1xyXG4gICAgICAgIGcuZmlsbFJlY3QoMCwgbllZLCBuV1csIG5IUm93R3JpZCk7XHJcbiAgICAgICAgZy5nbG9iYWxBbHBoYSA9IDE7XHJcbiAgICAgIH1cclxuICAgIH1cclxuICB9XHJcblxyXG5cclxuICBwcml2YXRlIHN0YXRpYyBoaXRUZXN0Um93cyhlQ2FudmFzUGlubmVkIDogSFRNTENhbnZhc0VsZW1lbnQsIGdyaWQgOiBERy5HcmlkLCBlIDogTW91c2VFdmVudCwgYkJvcmRlciA6IGJvb2xlYW4sIGFyWFlPbkNlbGwgOiBBcnJheTxudW1iZXI+IHwgdW5kZWZpbmVkKVxyXG4gIHtcclxuICAgIGNvbnN0IHJlY3QgPSBlQ2FudmFzUGlubmVkLmdldEJvdW5kaW5nQ2xpZW50UmVjdCgpO1xyXG4gICAgY29uc3Qgc2Nyb2xsTGVmdD0gd2luZG93LnBhZ2VYT2Zmc2V0IHx8IGRvY3VtZW50LmRvY3VtZW50RWxlbWVudC5zY3JvbGxMZWZ0O1xyXG4gICAgY29uc3Qgc2Nyb2xsVG9wID0gd2luZG93LnBhZ2VZT2Zmc2V0IHx8IGRvY3VtZW50LmRvY3VtZW50RWxlbWVudC5zY3JvbGxUb3A7XHJcbiAgICBjb25zdCBuWSA9IHJlY3QudG9wICArIHNjcm9sbFRvcDtcclxuICAgIGNvbnN0IG5YID0gcmVjdC5sZWZ0ICsgc2Nyb2xsTGVmdDtcclxuXHJcbiAgICBpZihuWCA8PSBlLmNsaWVudFggJiYgZS5jbGllbnRYIDw9IG5YICsgZUNhbnZhc1Bpbm5lZC5vZmZzZXRXaWR0aCkgICAvL29uIHRoZSByb3dzIGhlYWRlclxyXG4gICAge1xyXG4gICAgICBjb25zdCBuSEhlYWRlckNvbHMgPSBHcmlkVXRpbHMuZ2V0R3JpZENvbHVtbkhlYWRlckhlaWdodChncmlkKTtcclxuICAgICAgY29uc3QgbkhSb3dHcmlkID0gR3JpZFV0aWxzLmdldEdyaWRSb3dIZWlnaHQoZ3JpZCk7XHJcblxyXG4gICAgICBjb25zdCBhck1pbk1heFJvd3MgPSBbLTEsLTFdO1xyXG4gICAgICBHcmlkVXRpbHMuZmlsbFZpc2libGVWaWV3cG9ydFJvd3MoYXJNaW5NYXhSb3dzLCBncmlkKTtcclxuICAgICAgY29uc3QgblJvd01pbiA9IGFyTWluTWF4Um93c1swXTtcclxuICAgICAgY29uc3QgblJvd01heCA9IGFyTWluTWF4Um93c1sxXTtcclxuXHJcbiAgICAgIGNvbnN0IG5ZTW91c2VPbkhlYWRlciA9IGUuY2xpZW50WSAtIG5ZO1xyXG5cclxuICAgICAgbGV0IG5ZQm9yZGVyID0gLTE7XHJcbiAgICAgIGxldCBuWURpZmYgPSAtMTtcclxuXHJcbiAgICAgIGZvcihsZXQgblJvdz1uUm93TWluOyBuUm93PD0gblJvd01heDsgKytuUm93KVxyXG4gICAgICB7XHJcbiAgICAgICAgbllCb3JkZXIgPSBuSEhlYWRlckNvbHMgKyAoblJvdyAtIG5Sb3dNaW4rMSkqbkhSb3dHcmlkO1xyXG4gICAgICAgIG5ZRGlmZiA9IG5ZTW91c2VPbkhlYWRlciAtIG5ZQm9yZGVyO1xyXG5cclxuICAgICAgICBpZihiQm9yZGVyICYmIE1hdGguYWJzKG5ZRGlmZikgPD0gUGlubmVkQ29sdW1uLllfUkVTSVpFX1NFTlNJVElWSVRZKVxyXG4gICAgICAgIHtcclxuICAgICAgICAgIHJldHVybiBuUm93O1xyXG4gICAgICAgIH1cclxuXHJcbiAgICAgICAgaWYoIWJCb3JkZXIgJiYgbllCb3JkZXIgLSBuSFJvd0dyaWQgPD0gbllNb3VzZU9uSGVhZGVyICYmIG5ZTW91c2VPbkhlYWRlciA8PSBuWUJvcmRlcikge1xyXG5cclxuICAgICAgICAgIGlmKGFyWFlPbkNlbGwgIT09IHVuZGVmaW5lZCkge1xyXG4gICAgICAgICAgICBhclhZT25DZWxsWzBdID0gZS5jbGllbnRYIC0gblg7XHJcbiAgICAgICAgICAgIGFyWFlPbkNlbGxbMV0gPSBuWU1vdXNlT25IZWFkZXIgLSBuWUJvcmRlciArIG5IUm93R3JpZDtcclxuICAgICAgICAgIH1cclxuXHJcbiAgICAgICAgICByZXR1cm4gblJvdztcclxuICAgICAgICB9XHJcbiAgICAgIH1cclxuICAgIH1cclxuXHJcbiAgICByZXR1cm4gLTE7XHJcbiAgfVxyXG59XHJcbiJdfQ==