import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as GridUtils from '../utils/GridUtils';
import * as TextUtils from '../utils/TextUtils';
import { ColorUtils } from '../utils/ColorUtils';
import * as rxjs from 'rxjs';
import { GridCellRendererEx } from "../renderer/GridCellRendererEx";
import * as PinnedUtils from "./PinnedUtils";
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
export class PinnedColumn {
    constructor(colGrid) {
        const grid = getGrid(colGrid);
        if (grid === null) {
            throw new Error("Column '" + colGrid.name + "' is not attached to the grid.");
        }
        if (!PinnedUtils.isPinnableColumn(colGrid)) {
            throw new Error("Column '" + colGrid.name + "' cannot be pinned. It either pinned or HTML.");
        }
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
        this.m_observerResizeGrid = new ResizeObserver(entries => {
            const bCurrent = DG.toDart(grok.shell.v) === DG.toDart(viewTable);
            if (!bCurrent)
                return;
            if (this.m_fDevicePixelRatio !== window.devicePixelRatio || grid.canvas.height !== eCanvasThis.height) {
                eCanvasThis.width = nW * window.devicePixelRatio;
                eCanvasThis.height = grid.canvas.height;
                eCanvasThis.style.top = grid.canvas.offsetTop + "px";
                eCanvasThis.style.width = nW + "px";
                eCanvasThis.style.height = Math.round(grid.canvas.height / window.devicePixelRatio) + "px";
                this.m_fDevicePixelRatio = window.devicePixelRatio;
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
        this.m_observerResizeGrid.observe(grid.canvas);
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
        let nHResizeRowsBeforeDrag = -1;
        let nResizeRowGridDragging = -1;
        let nYResizeDraggingAnchor = -1;
        let nResizeRowGridMoving = -1;
        let nYDraggingAnchor = -1;
        let nRowGridDragging = -1;
        let arXYMouseOnCellDown = [-2, -2];
        let arXYMouseOnCellUp = [-1, -1];
        this.m_handlerMouseDown = rxjs.fromEvent(document, 'mousedown').subscribe((e) => {
            if (DG.toDart(grok.shell.v) !== DG.toDart(viewTable))
                return;
            const eMouse = e;
            if (eMouse.buttons !== 1)
                return;
            nResizeRowGridMoving = -1;
            const bAddToSel = eMouse.ctrlKey || eMouse.shiftKey;
            let nRowGrid = bAddToSel ? -1 : PinnedColumn.hitTestRows(eCanvasThis, grid, eMouse, true, undefined);
            if (nRowGrid >= 0) {
                const nHRows = GridUtils.getGridRowHeight(grid);
                nResizeRowGridDragging = nRowGrid;
                nYResizeDraggingAnchor = eMouse.clientY;
                nHResizeRowsBeforeDrag = nHRows;
            }
            else {
                nRowGrid = PinnedColumn.hitTestRows(eCanvasThis, grid, eMouse, false, arXYMouseOnCellDown);
                nRowGridDragging = nRowGrid;
                nYDraggingAnchor = eMouse.clientY;
                const cell = grid.cell(colGrid.name, nRowGrid);
                const renderer = getRenderer(cell);
                if (renderer instanceof GridCellRendererEx) {
                    renderer.onMouseDownEx(cell, eMouse, arXYMouseOnCellDown[0], arXYMouseOnCellDown[1]);
                }
            }
            e.preventDefault();
            e.stopPropagation();
        });
        this.m_handlerMouseUp = rxjs.fromEvent(document, 'mouseup').subscribe((e) => {
            if (DG.toDart(grok.shell.v) !== DG.toDart(viewTable)) {
                return;
            }
            const eMouse = e;
            if (nResizeRowGridDragging >= 0) {
                const nHRow = GridUtils.getGridRowHeight(grid);
                notifyAllPinnedColsRowsResized(headerThis, nHRow, false);
                notifyAllColsRowsResized(grid, nHRow, false);
            }
            nHResizeRowsBeforeDrag = -1;
            nResizeRowGridDragging = -1;
            nYResizeDraggingAnchor = -1;
            nResizeRowGridMoving = -1;
            document.body.style.cursor = "auto";
            if (nRowGridDragging >= 0) {
                const bAddToSel = eMouse.ctrlKey || eMouse.shiftKey;
                const nRowGrid = PinnedColumn.hitTestRows(eCanvasThis, grid, eMouse, false, arXYMouseOnCellUp);
                if (!bAddToSel && nRowGrid === nRowGridDragging) {
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
                    if (!bAddToSel)
                        bitsetSel.setAll(false, true);
                    const nRowMin = nRowGridDragging < nRowGrid ? nRowGridDragging : nRowGrid;
                    const nRowMax = nRowGridDragging > nRowGrid ? nRowGridDragging : nRowGrid;
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
                const cell = grid.cell(colGrid.name, nRowGrid);
                const renderer = getRenderer(cell);
                if (renderer instanceof GridCellRendererEx) {
                    renderer.onMouseUpEx(cell, eMouse, arXYMouseOnCellUp[0], arXYMouseOnCellUp[1]);
                }
                if (arXYMouseOnCellUp[0] === arXYMouseOnCellDown[0] && arXYMouseOnCellDown[1] === arXYMouseOnCellUp[1]) {
                    if (renderer instanceof GridCellRendererEx) {
                        renderer.onClickEx(cell, eMouse, arXYMouseOnCellUp[0], arXYMouseOnCellUp[1]);
                    }
                }
                nRowGridDragging = -1;
                nYDraggingAnchor = -1;
                arXYMouseOnCellDown[0] = -2;
                arXYMouseOnCellDown[1] = -2;
                arXYMouseOnCellUp[0] = -1;
                arXYMouseOnCellUp[1] = -1;
            }
            //e.preventDefault();
            //e.stopPropagation();
        });
        let cellCurrent = null;
        this.m_handlerMouseLeave = rxjs.fromEvent(document, 'mouseleave').subscribe((e) => {
            if (cellCurrent !== null) {
                const renderer = getRenderer(cellCurrent);
                if (renderer instanceof GridCellRendererEx) {
                    const eMouse = e;
                    renderer.onMouseLeaveEx(cellCurrent, eMouse, -1, -1);
                }
                cellCurrent = null;
            }
        });
        this.m_handlerMouseMove = rxjs.fromEvent(document, 'mousemove').subscribe((e) => {
            if (DG.toDart(grok.shell.v) !== DG.toDart(viewTable)) {
                return;
            }
            const bResizing = nResizeRowGridDragging >= 0;
            if (bResizing) {
                //console.log("Dragging : " + headerThis.m_strColName);
                const eMouse = e;
                const nYDiff = eMouse.clientY - nYResizeDraggingAnchor;
                let nHRowGrid = nHResizeRowsBeforeDrag + nYDiff;
                if (nHRowGrid < PinnedColumn.MIN_ROW_HEIGHT)
                    nHRowGrid = PinnedColumn.MIN_ROW_HEIGHT;
                else if (nHRowGrid > PinnedColumn.MAX_ROW_HEIGHT)
                    nHRowGrid = PinnedColumn.MAX_ROW_HEIGHT;
                let g = eCanvasThis.getContext('2d');
                if (g === null)
                    return;
                g.fillStyle = "white";
                const nHHeaderCols = GridUtils.getGridColumnHeaderHeight(grid);
                g.fillRect(0, nHHeaderCols, eCanvasThis.offsetWidth, eCanvasThis.offsetHeight);
                grid.setOptions({
                    rowHeight: nHRowGrid //this won't trigger onRowsRezized event, which is a DG bug
                });
                notifyAllPinnedColsRowsResized(headerThis, nHRowGrid, true);
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
            const arXYOnCell = [-1, -1];
            let nRowGrid = PinnedColumn.hitTestRows(eCanvasThis, grid, e, false, arXYOnCell);
            if (nRowGrid >= 0) {
                const cell = grid.cell(colGrid.name, nRowGrid);
                const renderer = getRenderer(cell);
                if (renderer instanceof GridCellRendererEx) {
                    if (cellCurrent === null) {
                        renderer.onMouseEnterEx(cell, e, arXYOnCell[0], arXYOnCell[1]);
                    }
                    if (cellCurrent !== null && nRowGrid !== cellCurrent.gridRow) {
                        renderer.onMouseLeaveEx(cellCurrent, e, -1, -1);
                        renderer.onMouseEnterEx(cell, e, arXYOnCell[0], arXYOnCell[1]);
                    }
                    renderer.onMouseMoveEx(cell, e, arXYOnCell[0], arXYOnCell[1]);
                }
                cellCurrent = cell;
            }
            else if (cellCurrent !== null) {
                const renderer = getRenderer(cellCurrent);
                if (renderer instanceof GridCellRendererEx) {
                    renderer.onMouseLeaveEx(cellCurrent, e, -1, -1);
                }
                cellCurrent = null;
            }
            nRowGrid = PinnedColumn.hitTestRows(eCanvasThis, grid, e, true, undefined);
            if (nRowGrid >= 0) {
                nResizeRowGridMoving = nRowGrid;
                document.body.style.cursor = "row-resize";
                return;
            }
            if (nResizeRowGridMoving >= 0) {
                nResizeRowGridMoving = -1;
                document.body.style.cursor = "auto";
            }
        });
        let nCount = 0;
        this.m_handlerMouseWheel = rxjs.fromEvent(this.m_root, 'wheel').subscribe((e) => {
            if (DG.toDart(grok.shell.v) !== DG.toDart(viewTable)) {
                return;
            }
            const eWheel = e;
            if (eWheel.deltaX !== 0 || eWheel.deltaZ !== 0) {
                return;
            }
            if (nCount === 1) {
                //scroll +
                const nRowCount = GridUtils.getGridVisibleRowCount(grid);
                const scrollY = grid.vertScroll;
                if (nRowCount - 1 > scrollY.max) {
                    scrollY.setValues(scrollY.minRange, scrollY.maxRange, scrollY.min + 1, scrollY.max + 1);
                }
                nCount = 0;
            }
            else if (nCount === -1) {
                //scroll -
                const scrollY = grid.vertScroll;
                if (scrollY.min >= 1) {
                    scrollY.setValues(scrollY.minRange, scrollY.maxRange, scrollY.min - 1, scrollY.max - 1);
                }
                nCount = 0;
            }
            else {
                nCount = eWheel.deltaY > 0 ? 1 : -1;
            }
            //console.log(eWheel.deltaX + " " + eWheel.deltaY);
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
        this.m_handlerMouseDown.unsubscribe();
        this.m_handlerMouseDown = null;
        this.m_handlerMouseUp.unsubscribe();
        this.m_handlerMouseUp = null;
        this.m_handlerMouseLeave.unsubscribe();
        this.m_handlerMouseLeave = null;
        this.m_handlerMouseMove.unsubscribe();
        this.m_handlerMouseMove = null;
        this.m_handlerMouseWheel.unsubscribe();
        this.m_handlerMouseWheel = null;
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
            bSel = nRowTable < 0 ? false : bitsetSel.get(nRowTable);
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
//# sourceMappingURL=data:application/json;base64,eyJ2ZXJzaW9uIjozLCJmaWxlIjoiUGlubmVkQ29sdW1uLmpzIiwic291cmNlUm9vdCI6IiIsInNvdXJjZXMiOlsiUGlubmVkQ29sdW1uLnRzIl0sIm5hbWVzIjpbXSwibWFwcGluZ3MiOiJBQUFBLE9BQU8sS0FBSyxJQUFJLE1BQU0sbUJBQW1CLENBQUM7QUFDMUMsT0FBTyxLQUFLLEVBQUUsTUFBTSxpQkFBaUIsQ0FBQztBQUN0QyxPQUFPLEtBQUssRUFBRSxNQUFNLGlCQUFpQixDQUFDO0FBQ3RDLE9BQU8sS0FBSyxTQUFTLE1BQU0sb0JBQW9CLENBQUM7QUFDaEQsT0FBTyxLQUFLLFNBQVMsTUFBTSxvQkFBb0IsQ0FBQztBQUNoRCxPQUFPLEVBQUMsVUFBVSxFQUFDLE1BQU0scUJBQXFCLENBQUM7QUFDL0MsT0FBTyxLQUFLLElBQUksTUFBTSxNQUFNLENBQUM7QUFDN0IsT0FBTyxFQUFFLGtCQUFrQixFQUFDLE1BQU0sZ0NBQWdDLENBQUM7QUFDbkUsT0FBTyxLQUFLLFdBQVcsTUFBTSxlQUFlLENBQUM7QUFDN0MsNENBQTRDO0FBRzVDOzs7Ozs7Ozs7Ozs7Ozs7O0VBZ0JFO0FBRUYsU0FBUyxXQUFXLENBQUMsSUFBa0I7SUFDckMsTUFBTSxPQUFPLEdBQUcsSUFBSSxDQUFDLFVBQVUsQ0FBQztJQUNoQyxJQUFJLE9BQU8sS0FBSyxJQUFJLElBQUksT0FBTyxLQUFLLFNBQVMsRUFBRTtRQUM3QyxNQUFNLElBQUksS0FBSyxDQUFDLDRDQUE0QyxDQUFDLENBQUM7S0FDL0Q7SUFFRCxJQUFJLFFBQVEsR0FBRyxTQUFTLENBQUMscUJBQXFCLENBQUMsT0FBTyxDQUFDLENBQUM7SUFDeEQsSUFBRyxRQUFRLFlBQVksa0JBQWtCLEVBQUU7UUFDekMsT0FBTyxRQUFRLENBQUM7S0FDakI7SUFFRCxPQUFPLElBQUksQ0FBQyxRQUFRLENBQUM7QUFDdkIsQ0FBQztBQUdELFNBQVMsT0FBTyxDQUFDLE9BQXVCO0lBQ3RDLElBQUksSUFBSSxHQUFvQixPQUFPLENBQUMsSUFBSSxDQUFDO0lBQ3pDLElBQUks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QixVQUFVLENBQUMsQ0FBQyxDQUFDLEdBQUcsQ0FBQyxDQUFDLE9BQU8sR0FBRyxFQUFFLENBQUM7d0JBQy9CLFVBQVUsQ0FBQyxDQUFDLENBQUMsR0FBRyxlQUFlLEdBQUcsUUFBUSxHQUFHLFNBQVMsQ0FBQztxQkFDeEQ7b0JBRUQsT0FBTyxJQUFJLENBQUM7aUJBQ2I7YUFDRjtTQUNGO1FBRUQsT0FBTyxDQUFDLENBQUMsQ0FBQztJQUNaLENBQUM7O0FBaDRCYywyQkFBYyxHQUFHLEVBQUUsQ0FBQztBQUNwQiwyQkFBYyxHQUFHLEdBQUcsQ0FBQztBQUNyQiw0QkFBZSxHQUFHLFVBQVUsQ0FBQyxLQUFLLENBQUMsVUFBVSxDQUFDLFlBQVksQ0FBQyxDQUFDLENBQUMsNkJBQTZCO0FBQzFGLDhCQUFpQixHQUFHLFVBQVUsQ0FBQyxLQUFLLENBQUMsVUFBVSxDQUFDLFVBQVUsQ0FBQyxDQUFDLENBQUMsNkJBQTZCO0FBQzFGLGlDQUFvQixHQUFHLENBQUMsQ0FBQyIsInNvdXJjZXNDb250ZW50IjpbImltcG9ydCAqIGFzIGdyb2sgZnJvbSAnZGF0YWdyb2stYXBpL2dyb2snO1xyXG5pbXBvcnQgKiBhcyBERyBmcm9tICdkYXRhZ3Jvay1hcGkvZGcnO1xyXG5pbXBvcnQgKiBhcyB1aSBmcm9tICdkYXRhZ3Jvay1hcGkvdWknO1xyXG5pbXBvcnQgKiBhcyBHcmlkVXRpbHMgZnJvbSAnLi4vdXRpbHMvR3JpZFV0aWxzJztcclxuaW1wb3J0ICogYXMgVGV4dFV0aWxzIGZyb20gJy4uL3V0aWxzL1RleHRVdGlscyc7XHJcbmltcG9ydCB7Q29sb3JVdGlsc30gZnJvbSAnLi4vdXRpbHMvQ29sb3JVdGlscyc7XHJcbmltcG9ydCAqIGFzIHJ4anMgZnJvbSAncnhqcyc7XHJcbmltcG9ydCB7IEdyaWRDZWxsUmVuZGVyZXJFeH0gZnJvbSBcIi4uL3JlbmRlcmVyL0dyaWRDZWxsUmVuZGVyZXJFeFwiO1xyXG5pbXBvcnQgKiBhcyBQaW5uZWRVdGlscyBmcm9tIFwiLi9QaW5uZWRVdGlsc1wiO1xyXG4vL2ltcG9ydCB7VGFibGVWaWV3fSBmcm9tIFwiZGF0YWdyb2stYXBpL2RnXCI7XHJcblxyXG5cclxuLypcclxuY29uc3QgaFN1YnNjcmliZXIgID0gZ3Jvay5ldmVudHMub25WaWV3TGF5b3V0QXBwbGllZC5zdWJzY3JpYmUoKGxheW91dCA6IERHLlZpZXdMYXlvdXQpID0+IHtcclxuICBjb25zdCB2aWV3IDogREcuVGFibGVWaWV3ID0gbGF5b3V0LnZpZXcgYXMgVGFibGVWaWV3O1xyXG4gIGNvbnN0IGl0Vmlld2VycyA9IHZpZXcudmlld2VycztcclxuICBjb25zdCBhclZpZXdlcnMgPSBBcnJheS5mcm9tKGl0Vmlld2Vycyk7XHJcblxyXG4gIGxldCB2aWV3ZXIgPSBudWxsO1xyXG4gIGNvbnN0IG5WaWV3ZXJDb3VudCA9IGFyVmlld2Vycy5sZW5ndGg7XHJcbiAgZm9yIChsZXQgbiA9IDA7IG4gPCBuVmlld2VyQ291bnQ7ICsrbikge1xyXG4gICAgdmlld2VyID0gYXJWaWV3ZXJzW25dO1xyXG4gICAgaWYgKHZpZXdlci50eXBlICE9PSBcIkdyaWRcIilcclxuICAgICAgY29udGludWU7XHJcblxyXG4gICAgUGlubmVkVXRpbHMuaW5zdGFsbFBpbm5lZENvbHVtbnModmlld2VyIGFzIERHLkdyaWQpO1xyXG4gIH1cclxufSk7XHJcbiovXHJcblxyXG5mdW5jdGlvbiBnZXRSZW5kZXJlcihjZWxsIDogREcuR3JpZENlbGwpIDogR3JpZENlbGxSZW5kZXJlckV4IHwgREcuR3JpZENlbGxSZW5kZXJlciB7XHJcbiAgY29uc3QgY29sR3JpZCA9IGNlbGwuZ3JpZENvbHVtbjtcclxuICBpZiAoY29sR3JpZCA9PT0gbnVsbCB8fCBjb2xHcmlkID09PSB1bmRlZmluZWQpIHtcclxuICAgIHRocm93IG5ldyBFcnJvcignR3JpZCBjZWxsIGlzIGRldGFjaGVkIGZyb20gdGhlIEdyaWQgY29sdW1uJyk7XHJcbiAgfVxyXG5cclxuICBsZXQgcmVuZGVyZXIgPSBHcmlkVXRpbHMuZ2V0R3JpZENvbHVtblJlbmRlcmVyKGNvbEdyaWQpO1xyXG4gIGlmKHJlbmRlcmVyIGluc3RhbmNlb2YgR3JpZENlbGxSZW5kZXJlckV4KSB7XHJcbiAgICByZXR1cm4gcmVuZGVyZXI7XHJcbiAgfVxyXG5cclxuICByZXR1cm4gY2VsbC5yZW5kZXJlcjtcclxufVxyXG5cclxuXHJcbmZ1bmN0aW9uIGdldEdyaWQoY29sR3JpZCA6IERHLkdyaWRDb2x1bW4pIDogREcuR3JpZCB8IG51bGwge1xyXG4gIGxldCBncmlkIDogREcuR3JpZCB8IG51bGwgPSBjb2xHcmlkLmdyaWQ7XHJcbiAgaWYoIGdyaWQgPT09IG51bGwpIHtcclxuICAgIGdyaWQgPSBHcmlkVXRpbHMuZ2V0SW5zdGFsbGVkR3JpZEZvckNvbHVtbihjb2xHcmlkKTtcclxuICAgIGlmKGdyaWQgaW5zdGFuY2VvZiBERy5HcmlkKVxyXG4gICAgICByZXR1cm4gZ3JpZDtcclxuICB9XHJcblxyXG4gIHJldHVybiBncmlkO1xyXG59XHJcblxyXG5cclxuZnVuY3Rpb24gbm90aWZ5QWxsQ29sc1Jvd3NSZXNpemVkKGdyaWQgOiBERy5HcmlkLCBuSFJvd3MgOiBudW1iZXIsIGJBZGp1c3RpbmcgOiBib29sZWFuKSA6IHZvaWQge1xyXG5cclxuICBsZXQgcmVuZGVyZXIgOiBHcmlkQ2VsbFJlbmRlcmVyRXggfCBudWxsID0gbnVsbFxyXG4gIGxldCBjb2xHcmlkID0gbnVsbDtcclxuICBjb25zdCBsc3RDb2xzR3JpZCA9IGdyaWQuY29sdW1ucztcclxuICBjb25zdCBuQ29sQ291bnQgPSBsc3RDb2xzR3JpZC5sZW5ndGg7XHJcbiAgZm9yKGxldCBuQ29sPTA7IG5Db2w8bkNvbENvdW50OyArK25Db2wpIHtcclxuICAgIGNvbEdyaWQgPSBsc3RDb2xzR3JpZC5ieUluZGV4KG5Db2wpO1xyXG4gICAgaWYoY29sR3JpZCA9PT0gbnVsbCB8fCAhY29sR3JpZC52aXNpYmxlKXtcclxuICAgICAgY29udGludWVcclxuICAgIH1cclxuXHJcbiAgICByZW5kZXJlciA9IEdyaWRVdGlscy5nZXRHcmlkQ29sdW1uUmVuZGVyZXIoY29sR3JpZCk7XHJcbiAgICBpZiAocmVuZGVyZXIgaW5zdGFuY2VvZiBHcmlkQ2VsbFJlbmRlcmVyRXgpIHtcclxuICAgICAgcmVuZGVyZXIub25SZXNpemVIZWlnaHQoY29sR3JpZCwgZ3JpZCwgbkhSb3dzLCBiQWRqdXN0aW5nKTtcclxuICAgIH1cclxuICB9XHJcbn1cclxuXHJcblxyXG5mdW5jdGlvbiBub3RpZnlBbGxQaW5uZWRDb2xzUm93c1Jlc2l6ZWQoY29sUGlubmVkU291cmNlIDogUGlubmVkQ29sdW1uLCBuSFJvd3MgOiBudW1iZXIsIGJBZGp1c3RpbmcgOiBib29sZWFuKSA6IHZvaWQge1xyXG5cclxuICBjb25zdCBjb2xHcmlkU291cmNlICA9IGNvbFBpbm5lZFNvdXJjZS5nZXRHcmlkQ29sdW1uKCk7XHJcbiAgaWYoY29sR3JpZFNvdXJjZSA9PT0gbnVsbCl7XHJcbiAgICByZXR1cm47XHJcbiAgfVxyXG5cclxuICBjb25zdCBncmlkID0gZ2V0R3JpZChjb2xHcmlkU291cmNlKTtcclxuICBjb25zdCBkYXJ0ID0gREcudG9EYXJ0KGdyaWQpO1xyXG4gIGlmKGRhcnQubV9hclBpbm5lZENvbHMgPT09IHVuZGVmaW5lZCkge1xyXG4gICAgdGhyb3cgbmV3IEVycm9yKCdQaW5uZWQgQ29sdW1ucyBhcmUgbm90IGluc3RhbGxlZC4nKTtcclxuICB9XHJcblxyXG4gIGxldCByZW5kZXJlciA6IEdyaWRDZWxsUmVuZGVyZXJFeCB8IG51bGwgPSBudWxsXHJcbiAgbGV0IGNvbFBpbm5lZCA9IG51bGw7XHJcbiAgbGV0IGNvbEdyaWQgPSBudWxsO1xyXG4gIGNvbnN0IG5QaW5uZWRDb2xDb3VudCA9IGRhcnQubV9hclBpbm5lZENvbHMubGVuZ3RoO1xyXG4gIGZvcihsZXQgbkNvbFBpbj0wOyBuQ29sUGluPG5QaW5uZWRDb2xDb3VudDsgKytuQ29sUGluKSB7XHJcbiAgICBjb2xQaW5uZWQgPSBkYXJ0Lm1fYXJQaW5uZWRDb2xzW25Db2xQaW5dO1xyXG4gICAgY29sR3JpZCA9IGNvbFBpbm5lZC5tX2NvbEdyaWQ7XHJcbiAgICBpZihjb2xHcmlkID09PSBudWxsKSB7XHJcbiAgICAgIHRocm93IG5ldyBFcnJvcignUGlubmVkIENvbHVtbiBpcyBkZXRhY2hlZC4nKTtcclxuICAgIH1cclxuXHJcbiAgICByZW5kZXJlciA9IEdyaWRVdGlscy5nZXRHcmlkQ29sdW1uUmVuZGVyZXIoY29sR3JpZCk7XHJcbiAgICBpZiAocmVuZGVyZXIgaW5zdGFuY2VvZiBHcmlkQ2VsbFJlbmRlcmVyRXggICYmIGNvbFBpbm5lZC5tX3Jvb3QgIT09IG51bGwgJiYgZ3JpZCAhPT0gbnVsbCkge1xyXG4gICAgICByZW5kZXJlci5vblJlc2l6ZUhlaWdodChjb2xQaW5uZWQsIGdyaWQsIG5IUm93cywgYkFkanVzdGluZyk7XHJcbiAgICB9XHJcbiAgfVxyXG59XHJcblxyXG5cclxuXHJcblxyXG5leHBvcnQgY2xhc3MgUGlubmVkQ29sdW1uIHtcclxuXHJcbiAgcHJpdmF0ZSBzdGF0aWMgTUlOX1JPV19IRUlHSFQgPSAyMDtcclxuICBwcml2YXRlIHN0YXRpYyBNQVhfUk9XX0hFSUdIVCA9IDUwMDtcclxuICBwcml2YXRlIHN0YXRpYyBTRUxFQ1RJT05fQ09MT1IgPSBDb2xvclV0aWxzLnRvUmdiKENvbG9yVXRpbHMuY29sU2VsZWN0aW9uKTsgLy9cInJnYmEoMjM3LCAyMjAsIDg4LCAwLjE1KVwiO1xyXG4gIHByaXZhdGUgc3RhdGljIEFDVElWRV9DRUxMX0NPTE9SID0gQ29sb3JVdGlscy50b1JnYihDb2xvclV0aWxzLmN1cnJlbnRSb3cpOyAvL1wicmdiYSgxNTMsIDIzNywgODIsIDAuMjUpXCI7XHJcbiAgcHJpdmF0ZSBzdGF0aWMgWV9SRVNJWkVfU0VOU0lUSVZJVFkgPSAyO1xyXG5cclxuICBwcml2YXRlIG1fZkRldmljZVBpeGVsUmF0aW8gOiBudW1iZXI7XHJcbiAgcHJpdmF0ZSBtX2NvbEdyaWQgOiBERy5HcmlkQ29sdW1uIHwgbnVsbDtcclxuICBwcml2YXRlIG1fcm9vdCA6IEhUTUxDYW52YXNFbGVtZW50IHwgbnVsbDtcclxuICBwcml2YXRlIG1fbldpZHRoQnVnIDogbnVtYmVyO1xyXG4gIC8vcHJpdmF0ZSBtX29ic2VydmVyUmVzaXplIDogUmVzaXplT2JzZXJ2ZXIgfCBudWxsO1xyXG4gIHByaXZhdGUgbV9vYnNlcnZlclJlc2l6ZUdyaWQgOiBSZXNpemVPYnNlcnZlciB8IG51bGw7XHJcbiAgcHJpdmF0ZSBtX2hhbmRsZXJWU2Nyb2xsIDogYW55O1xyXG4gIHByaXZhdGUgbV9oYW5kbGVyUm93c0ZpbHRlcmluZyA6IGFueTtcclxuICBwcml2YXRlIG1faGFuZGxlckN1cnJSb3cgOiBhbnk7XHJcbiAgcHJpdmF0ZSBtX2hhbmRsZXJTZWwgOiBhbnk7XHJcbiAgLy9wcml2YXRlIG1faGFuZGxlckZpbHRlciA6IGFueTtcclxuICBwcml2YXRlIG1faGFuZGxlclJvd3NSZXNpemVkIDogYW55O1xyXG4gIHByaXZhdGUgbV9oYW5kbGVyUm93c1NvcnRlZCA6IGFueTtcclxuICBwcml2YXRlIG1faGFuZGxlck1vdXNlRG93biA6IGFueTtcclxuICBwcml2YXRlIG1faGFuZGxlck1vdXNlVXAgOiBhbnk7XHJcbiAgcHJpdmF0ZSBtX2hhbmRsZXJNb3VzZUxlYXZlIDogYW55O1xyXG4gIHByaXZhdGUgbV9oYW5kbGVyTW91c2VNb3ZlIDogYW55O1xyXG4gIHByaXZhdGUgbV9oYW5kbGVyTW91c2VXaGVlbCA6IGFueTtcclxuXHJcbiAgY29uc3RydWN0b3IoY29sR3JpZCA6IERHLkdyaWRDb2x1bW4pIHtcclxuXHJcbiAgICBjb25zdCBncmlkID0gZ2V0R3JpZChjb2xHcmlkKTtcclxuICAgIGlmKGdyaWQgPT09IG51bGwpIHtcclxuICAgICAgdGhyb3cgbmV3IEVycm9yKFwiQ29sdW1uICdcIiArIGNvbEdyaWQubmFtZSArIFwiJyBpcyBub3QgYXR0YWNoZWQgdG8gdGhlIGdyaWQuXCIpO1xyXG4gICAgfVxyXG5cclxuICAgIGlmKCFQaW5uZWRVdGlscy5pc1Bpbm5hYmxlQ29sdW1uKGNvbEdyaWQpKSB7XHJcbiAgICAgIHRocm93IG5ldyBFcnJvcihcIkNvbHVtbiAnXCIgKyBjb2xHcmlkLm5hbWUgKyBcIicgY2Fubm90IGJlIHBpbm5lZC4gSXQgZWl0aGVyIHBpbm5lZCBvciBIVE1MLlwiKTtcclxuICAgIH1cclxuXHJcbiAgICB0aGlzLm1fZkRldmljZVBpeGVsUmF0aW8gPSB3aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbztcclxuXHJcbiAgICBjb25zdCBkYXJ0ID0gREcudG9EYXJ0KGdyaWQpO1xyXG5cclxuICAgIGlmKGRhcnQubV9hclBpbm5lZENvbHMgPT09IHVuZGVmaW5lZClcclxuICAgICAgZGFydC5tX2FyUGlubmVkQ29scyA9IFtdO1xyXG5cclxuICAgIGlmKGRhcnQubV9hclBpbm5lZENvbHMubGVuZ3RoID09PSAwICYmICFHcmlkVXRpbHMuaXNSb3dIZWFkZXIoY29sR3JpZCkpIHtcclxuICAgICAgY29uc3QgY29sR3JpZDAgPSBncmlkLmNvbHVtbnMuYnlJbmRleCgwKTtcclxuICAgICAgaWYoY29sR3JpZDAgIT09IG51bGwgJiYgY29sR3JpZDAgIT09IHVuZGVmaW5lZClcclxuICAgICAgbmV3IFBpbm5lZENvbHVtbihjb2xHcmlkMCk7XHJcbiAgICB9XHJcblxyXG4gICAgY29uc3QgbldUb3RhbFBpbm5lZENvbHMgPSBQaW5uZWRVdGlscy5nZXRUb3RhbFBpbm5lZENvbHNXaWR0aChncmlkKTtcclxuICAgIGRhcnQubV9hclBpbm5lZENvbHMucHVzaCh0aGlzKTtcclxuXHJcbiAgICBjb25zdCB2aWV3VGFibGUgPSBncmlkLnZpZXc7XHJcbiAgICBjb25zdCBkZnJhbWUgPSBncmlkLmRhdGFGcmFtZTtcclxuXHJcbiAgICBjb25zdCBuVyA9IGNvbEdyaWQud2lkdGg7XHJcbiAgICB0aGlzLm1fY29sR3JpZCA9IGNvbEdyaWQ7XHJcbiAgICB0aGlzLm1fbldpZHRoQnVnID0gLTE7XHJcbiAgICB0cnkge1xyXG4gICAgICBjb2xHcmlkLnZpc2libGUgPSBmYWxzZTtcclxuICAgIH1cclxuICAgIGNhdGNoKGUpIHtcclxuICAgICAgLy9ERyBidWdcclxuICAgICAgY29uc29sZS5lcnJvcihcIkVSUk9SOiBDb3VsZG4ndCBoaWRlIGNvbHVtbiAnXCIgKyBjb2xHcmlkLm5hbWUgKyBcIicgZHVlIHRvIGEgREcgYnVnLiBBdHRlbXB0IHRvIHNldCB0aGUgd2lkdGggdG8gMFwiKTtcclxuICAgICAgdHJ5IHtcclxuICAgICAgICB0aGlzLm1fbldpZHRoQnVnID0gY29sR3JpZC53aWR0aDtcclxuICAgICAgICBjb2xHcmlkLndpZHRoID0gMDtcclxuICAgICAgfSBjYXRjaCAoZSkge1xyXG4gICAgICAgIC8vREcgYnVnXHJcbiAgICAgICAgY29uc29sZS5lcnJvcihcIkVSUk9SOiBDb3VsZG4ndCBzZXQgdGhlIHdpZHRoIHRvIDAgZm9yIGNvbHVtbiAnXCIgKyBjb2xHcmlkLm5hbWUgKyBcIicgZHVlIHRvIGEgREcgYnVnLiBUaGlzIGNvdWxkIGJlIGlnbm9yZWQgaWYgdGhlIGNvbHVtbiB2aXN1YWxseSBsb29rcyBvay5cIik7XHJcbiAgICAgIH1cclxuICAgIH1cclxuXHJcbiAgICBpZighR3JpZFV0aWxzLmlzUm93SGVhZGVyKGNvbEdyaWQpKSB7XHJcbiAgICAgIGlmIChjb2xHcmlkLnNldHRpbmdzID09PSBudWxsIHx8IGNvbEdyaWQuc2V0dGluZ3MgPT09IHVuZGVmaW5lZClcclxuICAgICAgICBjb2xHcmlkLnNldHRpbmdzID0ge307XHJcblxyXG4gICAgICBjb2xHcmlkLnNldHRpbmdzLmlzUGlubmVkID0gdHJ1ZTsgLy90aGlzIHdpbGwgYmUgc2F2ZWQgd2l0aCB0aGUgbGF5b3V0XHJcbiAgICAgIGNvbEdyaWQuc2V0dGluZ3MuaWR4UGlubmVkID0gZGFydC5tX2FyUGlubmVkQ29scy5sZW5ndGggLSAxO1xyXG4gICAgfVxyXG5cclxuICAgIGdyaWQuY2FudmFzLnN0eWxlLmxlZnQgPSAoZ3JpZC5jYW52YXMub2Zmc2V0TGVmdCArIG5XKS50b1N0cmluZygpICsgXCJweFwiO1xyXG4gICAgZ3JpZC5vdmVybGF5LnN0eWxlLmxlZnQ9IChncmlkLm92ZXJsYXkub2Zmc2V0TGVmdCArIG5XKS50b1N0cmluZygpICsgXCJweFwiO1xyXG5cclxuICAgIGdyaWQuY2FudmFzLnN0eWxlLndpZHRoID0gKGdyaWQuY2FudmFzLm9mZnNldFdpZHRoIC0gblcpLnRvU3RyaW5nKCkgKyBcInB4XCI7XHJcbiAgICBncmlkLm92ZXJsYXkuc3R5bGUud2lkdGg9IChncmlkLm92ZXJsYXkub2Zmc2V0V2lkdGggLSBuVykudG9TdHJpbmcoKSArIFwicHhcIjtcclxuXHJcbiAgICBjb25zdCBuSGVpZ2h0ID0gZ3JpZC5jYW52YXMuaGVpZ2h0Oy8vY2FudmFzIHBpeGVsIGhlaWdodFxyXG4gICAgY29uc3QgZUNhbnZhc1RoaXMgPSB1aS5jYW52YXMoblcqd2luZG93LmRldmljZVBpeGVsUmF0aW8sIG5IZWlnaHQpO1xyXG4gICAgY29uc3QgdGFiSW5kZXggPSAgZ3JpZC5jYW52YXMuZ2V0QXR0cmlidXRlKFwidGFiSW5kZXhcIik7XHJcbiAgICBpZih0YWJJbmRleCAhPT0gbnVsbClcclxuICAgICBlQ2FudmFzVGhpcy5zZXRBdHRyaWJ1dGUoXCJ0YWJJbmRleFwiLCB0YWJJbmRleCk7XHJcblxyXG4gICAgZUNhbnZhc1RoaXMuc3R5bGUucG9zaXRpb24gPSBcImFic29sdXRlXCI7XHJcbiAgICBlQ2FudmFzVGhpcy5zdHlsZS5sZWZ0ID0gbldUb3RhbFBpbm5lZENvbHMgKyBcInB4XCI7XHJcbiAgICBlQ2FudmFzVGhpcy5zdHlsZS50b3AgPSBncmlkLmNhbnZhcy5vZmZzZXRUb3AgKyBcInB4XCI7XHJcbiAgICBlQ2FudmFzVGhpcy5zdHlsZS53aWR0aCA9IG5XICsgXCJweFwiO1xyXG4gICAgZUNhbnZhc1RoaXMuc3R5bGUuaGVpZ2h0ID0gTWF0aC5yb3VuZChuSGVpZ2h0L3dpbmRvdy5kZXZpY2VQaXhlbFJhdGlvKSArIFwicHhcIjtcclxuXHJcbiAgICAvL2NvbnNvbGUubG9nKFwiaCBcIiArIGdyaWQuY2FudmFzLmhlaWdodCArIFwiIG9mZnNldCBcIiArIGdyaWQuY2FudmFzLm9mZnNldEhlaWdodCk7XHJcblxyXG4gICAgaWYoZ3JpZC5jYW52YXMucGFyZW50Tm9kZSA9PT0gbnVsbClcclxuICAgICAgdGhyb3cgbmV3IEVycm9yKFwiUGFyZW50IG5vZGUgZm9yIGNhbnZhcyBjYW5ub3QgYmUgbnVsbC5cIik7XHJcblxyXG4gICAgZ3JpZC5jYW52YXMucGFyZW50Tm9kZS5pbnNlcnRCZWZvcmUoZUNhbnZhc1RoaXMsIGdyaWQuY2FudmFzKTtcclxuICAgIHRoaXMubV9yb290ID0gZUNhbnZhc1RoaXM7XHJcblxyXG5cclxuICAgIGNvbnN0IGNvbEdyaWQwID0gZ3JpZC5jb2x1bW5zLmJ5SW5kZXgoMCk7XHJcbiAgICBpZihjb2xHcmlkMCAhPT0gbnVsbCAmJiBjb2xHcmlkMCAhPT0gdW5kZWZpbmVkKSB7Ly9ERyBCdWcgZnJvbSByZWFkaW5nIGxheW91dFxyXG4gICAgdHJ5e1xyXG4gICAgICAgIGNvbEdyaWQwLnZpc2libGUgPSBmYWxzZTtcclxuICAgICAgfVxyXG4gICAgICBjYXRjaChlKSB7XHJcbiAgICAgICAgY29uc29sZS5lcnJvcihcIkVSUk9SOiBDb3VsZG4ndCBoaWRlIHJvdyBoZWFkZXIuXCIpO1xyXG4gICAgICB9XHJcbiAgICB9XHJcblxyXG5cclxuICAgIC8vT25SZXNpemUgUm93IGhlYWRlclxyXG4gICAgY29uc3QgaGVhZGVyVGhpcyA9IHRoaXM7LypcclxuICAgIHRoaXMubV9vYnNlcnZlclJlc2l6ZSA9IG5ldyBSZXNpemVPYnNlcnZlcihlbnRyaWVzID0+IHtcclxuICAgICAgY29uc3QgZyA9IGhlYWRlclRoaXMubV9yb290LmdldENvbnRleHQoJzJkJyk7XHJcbiAgICAgIGZvciAobGV0IGVudHJ5IG9mIGVudHJpZXMpIHtcclxuICAgICAgICBoZWFkZXJUaGlzLnBhaW50KGcsIGdyaWQpO1xyXG4gICAgICB9XHJcbiAgICB9KTtcclxuICAgIHRoaXMubV9vYnNlcnZlclJlc2l6ZS5vYnNlcnZlKGhlYWRlclRoaXMubV9yb290KTsqL1xyXG5cclxuXHJcblxyXG4gICAgLy9PblJlc2l6ZSBHcmlkXHJcbiAgICB0aGlzLm1fb2JzZXJ2ZXJSZXNpemVHcmlkID0gbmV3IFJlc2l6ZU9ic2VydmVyKGVudHJpZXMgPT4ge1xyXG5cclxuICAgICAgY29uc3QgYkN1cnJlbnQgPSAgREcudG9EYXJ0KGdyb2suc2hlbGwudikgPT09IERHLnRvRGFydCh2aWV3VGFibGUpO1xyXG4gICAgICBpZighYkN1cnJlbnQpXHJcbiAgICAgICAgcmV0dXJuO1xyXG5cclxuICAgICAgaWYodGhpcy5tX2ZEZXZpY2VQaXhlbFJhdGlvICE9PSB3aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbyB8fCBncmlkLmNhbnZhcy5oZWlnaHQgIT09IGVDYW52YXNUaGlzLmhlaWdodCkge1xyXG4gICAgICAgIGVDYW52YXNUaGlzLndpZHRoID0gblcqd2luZG93LmRldmljZVBpeGVsUmF0aW87XHJcbiAgICAgICAgZUNhbnZhc1RoaXMuaGVpZ2h0ID0gZ3JpZC5jYW52YXMuaGVpZ2h0O1xyXG4gICAgICAgIGVDYW52YXNUaGlzLnN0eWxlLnRvcCA9IGdyaWQuY2FudmFzLm9mZnNldFRvcCArIFwicHhcIjtcclxuICAgICAgICBlQ2FudmFzVGhpcy5zdHlsZS53aWR0aCA9IG5XICsgXCJweFwiO1xyXG4gICAgICAgIGVDYW52YXNUaGlzLnN0eWxlLmhlaWdodCA9IE1hdGgucm91bmQoZ3JpZC5jYW52YXMuaGVpZ2h0L3dpbmRvdy5kZXZpY2VQaXhlbFJhdGlvKSArIFwicHhcIjtcclxuXHJcbiAgICAgICAgdGhpcy5tX2ZEZXZpY2VQaXhlbFJhdGlvID0gd2luZG93LmRldmljZVBpeGVsUmF0aW87XHJcbiAgICAgIH1cclxuXHJcbiAgICAgIC8vY29uc29sZS5sb2coXCJHcmlkIFJlc2l6ZTogXCIgKyBncmlkLmNhbnZhcy5oZWlnaHQgKyBcIiBcIiArIHdpbmRvdy5kZXZpY2VQaXhlbFJhdGlvKTtcclxuICAgICAgLy9lQ2FudmFzVGhpcy5zdHlsZS5oZWlnaHQgPSBncmlkLnJvb3Quc3R5bGUuaGVpZ2h0O1xyXG4vKlxyXG4gICAgICBjb25zdCBlQ2FudmFzTmV3ID0gdWkuY2FudmFzKG5XLCBncmlkLnJvb3Qub2Zmc2V0SGVpZ2h0KTtcclxuICAgICAgaWYoaGVhZGVyVGhpcy5tX3Jvb3QucGFyZW50Tm9kZSAhPT0gbnVsbCkge1xyXG4gICAgICAgIGhlYWRlclRoaXMubV9yb290LnBhcmVudE5vZGUucmVwbGFjZUNoaWxkKGVDYW52YXNOZXcsIGhlYWRlclRoaXMubV9yb290KTtcclxuICAgICAgICBoZWFkZXJUaGlzLm1fcm9vdCA9IGVDYW52YXNOZXc7XHJcbiAgICAgIH0qL1xyXG4gICAgICAvL2hlYWRlclRoaXMubV9yb290LmhlaWdodCA9IGdyaWQucm9vdC5vZmZzZXRIZWlnaHQ7XHJcbiAgICAgIGNvbnN0IGcgPSBlQ2FudmFzVGhpcy5nZXRDb250ZXh0KCcyZCcpO1xyXG4gICAgICBmb3IgKGxldCBlbnRyeSBvZiBlbnRyaWVzKSB7XHJcbiAgICAgICAgc2V0VGltZW91dCgoKT0+IHtoZWFkZXJUaGlzLnBhaW50KGcsIGdyaWQpO30sIDEwMCk7XHJcbiAgICAgIH1cclxuICAgIH0pO1xyXG5cclxuICAgIHRoaXMubV9vYnNlcnZlclJlc2l6ZUdyaWQub2JzZXJ2ZShncmlkLmNhbnZhcyk7XHJcblxyXG4gICAgY29uc3Qgc2Nyb2xsVmVydCA9IGdyaWQudmVydFNjcm9sbDtcclxuICAgIHRoaXMubV9oYW5kbGVyVlNjcm9sbCA9IHNjcm9sbFZlcnQub25WYWx1ZXNDaGFuZ2VkLnN1YnNjcmliZSgoKSA9PiB7XHJcbiAgICAgIGNvbnN0IGcgPSBlQ2FudmFzVGhpcy5nZXRDb250ZXh0KCcyZCcpO1xyXG4gICAgICBoZWFkZXJUaGlzLnBhaW50KGcsIGdyaWQpO1xyXG4gICAgfSk7XHJcblxyXG4gICAgdGhpcy5tX2hhbmRsZXJSb3dzRmlsdGVyaW5nID0gZGZyYW1lLm9uUm93c0ZpbHRlcmluZy5zdWJzY3JpYmUoKCkgPT4ge1xyXG4gICAgICBzZXRUaW1lb3V0KCgpID0+IHtcclxuICAgICAgICBjb25zdCBnID0gZUNhbnZhc1RoaXMuZ2V0Q29udGV4dCgnMmQnKTtcclxuICAgICAgICBoZWFkZXJUaGlzLnBhaW50KGcsIGdyaWQpO1xyXG4gICAgICB9LCAxMDApO1xyXG5cclxuICAgIH0pO1xyXG5cclxuICAgIHRoaXMubV9oYW5kbGVyQ3VyclJvdyA9IGRmcmFtZS5vbkN1cnJlbnRSb3dDaGFuZ2VkLnN1YnNjcmliZSgoKSA9PiB7XHJcbiAgICAgICAgY29uc3QgZyA9IGVDYW52YXNUaGlzLmdldENvbnRleHQoJzJkJyk7XHJcbiAgICAgICAgaGVhZGVyVGhpcy5wYWludChnLCBncmlkKTtcclxuICAgICAgfVxyXG4gICAgKTtcclxuXHJcbiAgICB0aGlzLm1faGFuZGxlclNlbCA9IGRmcmFtZS5vblNlbGVjdGlvbkNoYW5nZWQuc3Vic2NyaWJlKChlIDogYW55KSA9PiB7XHJcbiAgICAgICAgY29uc3QgZyA9IGVDYW52YXNUaGlzLmdldENvbnRleHQoJzJkJyk7XHJcbiAgICAgICAgaGVhZGVyVGhpcy5wYWludChnLCBncmlkKTtcclxuICAgICAgfVxyXG4gICAgKTtcclxuXHJcbi8qXHJcbiAgICB0aGlzLm1faGFuZGxlckZpbHRlciA9IGRmcmFtZS5vblJvd3NGaWx0ZXJlZC5zdWJzY3JpYmUoKGUgOiBhbnkpID0+IHtcclxuICAgICAgICBjb25zdCBnID0gZUNhbnZhc1RoaXMuZ2V0Q29udGV4dCgnMmQnKTtcclxuICAgICAgICBoZWFkZXJUaGlzLnBhaW50KGcsIGdyaWQpO1xyXG4gICAgICB9XHJcbiAgICApO1xyXG4qL1xyXG5cclxuICAgIHRoaXMubV9oYW5kbGVyUm93c1Jlc2l6ZWQgPSBncmlkLm9uUm93c1Jlc2l6ZWQuc3Vic2NyaWJlKChlIDogYW55KSA9PiB7XHJcbiAgICAgICAgY29uc3QgZyA9IGVDYW52YXNUaGlzLmdldENvbnRleHQoJzJkJyk7XHJcbiAgICAgICAgaGVhZGVyVGhpcy5wYWludChnLCBncmlkKTtcclxuICAgICAgfVxyXG4gICAgKTtcclxuXHJcbiAgICB0aGlzLm1faGFuZGxlclJvd3NTb3J0ZWQgPSBncmlkLm9uUm93c1NvcnRlZC5zdWJzY3JpYmUoKGUgOiBhbnkpID0+IHtcclxuICAgICAgICBjb25zdCBnID0gZUNhbnZhc1RoaXMuZ2V0Q29udGV4dCgnMmQnKTtcclxuICAgICAgICBoZWFkZXJUaGlzLnBhaW50KGcsIGdyaWQpO1xyXG4gICAgICB9XHJcbiAgICApO1xyXG5cclxuICAgIGxldCBuSFJlc2l6ZVJvd3NCZWZvcmVEcmFnICA9IC0xO1xyXG4gICAgbGV0IG5SZXNpemVSb3dHcmlkRHJhZ2dpbmcgPSAtMTtcclxuICAgIGxldCBuWVJlc2l6ZURyYWdnaW5nQW5jaG9yID0gLTE7XHJcbiAgICBsZXQgblJlc2l6ZVJvd0dyaWRNb3ZpbmcgPSAtMTtcclxuXHJcbiAgICBsZXQgbllEcmFnZ2luZ0FuY2hvciA9IC0xO1xyXG4gICAgbGV0IG5Sb3dHcmlkRHJhZ2dpbmcgPSAtMTtcclxuXHJcbiAgICBsZXQgYXJYWU1vdXNlT25DZWxsRG93biA9IFstMiwgLTJdO1xyXG4gICAgbGV0IGFyWFlNb3VzZU9uQ2VsbFVwID0gWy0xLCAtMV07XHJcblxyXG4gICAgdGhpcy5tX2hhbmRsZXJNb3VzZURvd24gPSByeGpzLmZyb21FdmVudChkb2N1bWVudCwgJ21vdXNlZG93bicpLnN1YnNjcmliZSgoZSA6IEV2ZW50KSA9PiB7XHJcblxyXG4gICAgICBpZihERy50b0RhcnQoZ3Jvay5zaGVsbC52KSAhPT0gREcudG9EYXJ0KHZpZXdUYWJsZSkpXHJcbiAgICAgICAgcmV0dXJuO1xyXG5cclxuICAgICAgY29uc3QgZU1vdXNlID0gZSBhcyBNb3VzZUV2ZW50O1xyXG5cclxuICAgICAgaWYoZU1vdXNlLmJ1dHRvbnMgIT09IDEpXHJcbiAgICAgICAgcmV0dXJuO1xyXG5cclxuICAgICAgblJlc2l6ZVJvd0dyaWRNb3ZpbmcgPSAtMTtcclxuXHJcbiAgICAgIGNvbnN0IGJBZGRUb1NlbCA6IGJvb2xlYW4gPSBlTW91c2UuY3RybEtleSB8fCBlTW91c2Uuc2hpZnRLZXk7XHJcblxyXG4gICAgICBsZXQgblJvd0dyaWQgPSBiQWRkVG9TZWwgPyAtMSA6IFBpbm5lZENvbHVtbi5oaXRUZXN0Um93cyhlQ2FudmFzVGhpcywgZ3JpZCwgZU1vdXNlLCB0cnVlLCB1bmRlZmluZWQpO1xyXG4gICAgICBpZiAoblJvd0dyaWQgPj0gMCkge1xyXG4gICAgICAgIGNvbnN0IG5IUm93cyA9IEdyaWRVdGlscy5nZXRHcmlkUm93SGVpZ2h0KGdyaWQpO1xyXG4gICAgICAgIG5SZXNpemVSb3dHcmlkRHJhZ2dpbmcgPSBuUm93R3JpZDtcclxuICAgICAgICBuWVJlc2l6ZURyYWdnaW5nQW5jaG9yID0gZU1vdXNlLmNsaWVudFk7XHJcbiAgICAgICAgbkhSZXNpemVSb3dzQmVmb3JlRHJhZyA9IG5IUm93cztcclxuICAgICAgfVxyXG4gICAgICBlbHNlXHJcbiAgICAgIHtcclxuICAgICAgICBuUm93R3JpZCA9IFBpbm5lZENvbHVtbi5oaXRUZXN0Um93cyhlQ2FudmFzVGhpcywgZ3JpZCwgZU1vdXNlLCBmYWxzZSwgYXJYWU1vdXNlT25DZWxsRG93bik7XHJcblxyXG4gICAgICAgIG5Sb3dHcmlkRHJhZ2dpbmcgPSBuUm93R3JpZDtcclxuICAgICAgICBuWURyYWdnaW5nQW5jaG9yID0gZU1vdXNlLmNsaWVudFk7XHJcblxyXG4gICAgICAgIGNvbnN0IGNlbGwgPSBncmlkLmNlbGwoY29sR3JpZC5uYW1lLCBuUm93R3JpZCk7XHJcbiAgICAgICAgY29uc3QgcmVuZGVyZXIgPSBnZXRSZW5kZXJlcihjZWxsKTtcclxuICAgICAgICBpZihyZW5kZXJlciBpbnN0YW5jZW9mIEdyaWRDZWxsUmVuZGVyZXJFeCkge1xyXG4gICAgICAgICAgcmVuZGVyZXIub25Nb3VzZURvd25FeChjZWxsLCBlTW91c2UsIGFyWFlNb3VzZU9uQ2VsbERvd25bMF0sIGFyWFlNb3VzZU9uQ2VsbERvd25bMV0pO1xyXG4gICAgICAgIH1cclxuICAgICAgfVxyXG5cclxuICAgICAgZS5wcmV2ZW50RGVmYXVsdCgpO1xyXG4gICAgICBlLnN0b3BQcm9wYWdhdGlvbigpO1xyXG4gICAgfSk7XHJcblxyXG5cclxuICAgIHRoaXMubV9oYW5kbGVyTW91c2VVcCA9IHJ4anMuZnJvbUV2ZW50KGRvY3VtZW50LCAnbW91c2V1cCcpLnN1YnNjcmliZSgoZSkgPT4ge1xyXG5cclxuICAgICAgaWYoREcudG9EYXJ0KGdyb2suc2hlbGwudikgIT09IERHLnRvRGFydCh2aWV3VGFibGUpKSB7XHJcbiAgICAgICAgcmV0dXJuO1xyXG4gICAgICB9XHJcblxyXG4gICAgICBjb25zdCBlTW91c2UgPSBlIGFzIE1vdXNlRXZlbnQ7XHJcblxyXG4gICAgICBpZihuUmVzaXplUm93R3JpZERyYWdnaW5nID49IDApIHtcclxuICAgICAgICBjb25zdCBuSFJvdyA9IEdyaWRVdGlscy5nZXRHcmlkUm93SGVpZ2h0KGdyaWQpO1xyXG4gICAgICAgIG5vdGlmeUFsbFBpbm5lZENvbHNSb3dzUmVzaXplZChoZWFkZXJUaGlzLCBuSFJvdywgZmFsc2UpO1xyXG4gICAgICAgIG5vdGlmeUFsbENvbHNSb3dzUmVzaXplZChncmlkLCBuSFJvdywgZmFsc2UpO1xyXG4gICAgICB9XHJcblxyXG5cclxuICAgICAgbkhSZXNpemVSb3dzQmVmb3JlRHJhZyA9IC0xO1xyXG4gICAgICBuUmVzaXplUm93R3JpZERyYWdnaW5nID0gLTE7XHJcbiAgICAgIG5ZUmVzaXplRHJhZ2dpbmdBbmNob3IgPSAtMTtcclxuICAgICAgblJlc2l6ZVJvd0dyaWRNb3ZpbmcgPSAtMTtcclxuXHJcbiAgICAgIGRvY3VtZW50LmJvZHkuc3R5bGUuY3Vyc29yID0gXCJhdXRvXCI7XHJcblxyXG4gICAgICBpZihuUm93R3JpZERyYWdnaW5nID49IDApIHtcclxuICAgICAgICBjb25zdCBiQWRkVG9TZWwgPSBlTW91c2UuY3RybEtleSB8fCBlTW91c2Uuc2hpZnRLZXk7XHJcblxyXG4gICAgICAgIGNvbnN0IG5Sb3dHcmlkID0gUGlubmVkQ29sdW1uLmhpdFRlc3RSb3dzKGVDYW52YXNUaGlzLCBncmlkLCBlTW91c2UsIGZhbHNlLCBhclhZTW91c2VPbkNlbGxVcCk7XHJcbiAgICAgICAgaWYoIWJBZGRUb1NlbCAmJiBuUm93R3JpZCA9PT0gblJvd0dyaWREcmFnZ2luZykge1xyXG5cclxuICAgICAgICAgIGxldCBjZWxsUkggPSBudWxsO1xyXG4gICAgICAgICAgdHJ5IHtcclxuICAgICAgICAgICAgY2VsbFJIID0gZ3JpZC5jZWxsKFwiXCIsIG5Sb3dHcmlkKTtcclxuICAgICAgICAgIH1cclxuICAgICAgICAgIGNhdGNoKGUpIHtcclxuICAgICAgICAgICAgbGV0IGNvbEcgPSBudWxsO1xyXG4gICAgICAgICAgICBjb25zdCBsc3RDb2xzID0gZ3JpZC5jb2x1bW5zO1xyXG4gICAgICAgICAgICBmb3IobGV0IG5DPTE7IG5DPGxzdENvbHMubGVuZ3RoOyArK25DKSB7XHJcbiAgICAgICAgICAgICAgY29sRyA9IGxzdENvbHMuYnlJbmRleChuQyk7XHJcbiAgICAgICAgICAgICAgY2VsbFJIID0gY29sRyA9PT0gbnVsbCA/IG51bGwgOiBncmlkLmNlbGwoY29sRy5uYW1lLCBuUm93R3JpZCk7XHJcbiAgICAgICAgICAgICAgaWYoY2VsbFJIICE9PSBudWxsKVxyXG4gICAgICAgICAgICAgICAgYnJlYWs7XHJcbiAgICAgICAgICAgIH1cclxuICAgICAgICAgIH1cclxuICAgICAgICAgIGlmKGNlbGxSSCAhPT0gbnVsbCkge1xyXG4gICAgICAgICAgICBjb25zdCBuUm93VGFibGUgOiBhbnkgPSBjZWxsUkgudGFibGVSb3dJbmRleDtcclxuICAgICAgICAgICAgaWYoblJvd1RhYmxlICE9PSBudWxsKVxyXG4gICAgICAgICAgICBkZnJhbWUuY3VycmVudFJvdyA9IG5Sb3dUYWJsZTtcclxuICAgICAgICAgIH1cclxuICAgICAgICB9XHJcbiAgICAgICAgZWxzZVxyXG4gICAgICAgIHtcclxuICAgICAgICAgIGNvbnN0IGJpdHNldFNlbCA9IGRmcmFtZS5zZWxlY3Rpb247XHJcblxyXG4gICAgICAgICAgaWYoIWJBZGRUb1NlbClcclxuICAgICAgICAgICAgYml0c2V0U2VsLnNldEFsbChmYWxzZSwgdHJ1ZSk7XHJcblxyXG4gICAgICAgICAgY29uc3QgblJvd01pbiA9IG5Sb3dHcmlkRHJhZ2dpbmcgPCBuUm93R3JpZCA/IG5Sb3dHcmlkRHJhZ2dpbmcgOiBuUm93R3JpZDtcclxuICAgICAgICAgIGNvbnN0IG5Sb3dNYXggPSBuUm93R3JpZERyYWdnaW5nID4gblJvd0dyaWQgPyBuUm93R3JpZERyYWdnaW5nIDogblJvd0dyaWQ7XHJcbiAgICAgICAgICBsZXQgY2VsbFJIID0gbnVsbDtcclxuICAgICAgICAgIGxldCBuUm93VGFibGUgPSAtMTtcclxuICAgICAgICAgIGZvcihsZXQgblJvdz1uUm93TWluOyBuUm93PD1uUm93TWF4OyArK25Sb3cpIHtcclxuXHJcbiAgICAgICAgICAgIHRyeSB7XHJcbiAgICAgICAgICAgICAgY2VsbFJIID0gZ3JpZC5jZWxsKFwiXCIsIG5Sb3cpO1xyXG4gICAgICAgICAgICB9XHJcbiAgICAgICAgICAgIGNhdGNoKGUpIHtcclxuICAgICAgICAgICAgICBsZXQgY29sRyA9IG51bGw7XHJcbiAgICAgICAgICAgICAgY29uc3QgbHN0Q29scyA9IGdyaWQuY29sdW1ucztcclxuICAgICAgICAgICAgICBmb3IobGV0IG5DPTE7IG5DPGxzdENvbHMubGVuZ3RoOyArK25DKSB7XHJcbiAgICAgICAgICAgICAgICBjb2xHID0gbHN0Q29scy5ieUluZGV4KG5DKTtcclxuICAgICAgICAgICAgICAgIGNlbGxSSCA9IGNvbEcgPT09IG51bGwgPyBudWxsIDogZ3JpZC5jZWxsKGNvbEcubmFtZSwgblJvd0dyaWQpO1xyXG4gICAgICAgICAgICAgICAgaWYoY2VsbFJIICE9PSBudWxsKVxyXG4gICAgICAgICAgICAgICAgICBicmVhaztcclxuICAgICAgICAgICAgICB9XHJcbiAgICAgICAgICAgIH1cclxuXHJcbiAgICAgICAgICAgIGlmKGNlbGxSSCAhPT0gbnVsbCAmJiBjZWxsUkgudGFibGVSb3dJbmRleCAhPT0gbnVsbCkge1xyXG4gICAgICAgICAgICAgIG5Sb3dUYWJsZSA9IGNlbGxSSC50YWJsZVJvd0luZGV4O1xyXG4gICAgICAgICAgICAgIGJpdHNldFNlbC5zZXQoblJvd1RhYmxlLCB0cnVlLCB0cnVlKTtcclxuICAgICAgICAgICAgfVxyXG4gICAgICAgICAgfVxyXG4gICAgICAgIH1cclxuXHJcbiAgICAgICAgY29uc3QgY2VsbCA9IGdyaWQuY2VsbChjb2xHcmlkLm5hbWUsIG5Sb3dHcmlkKTtcclxuICAgICAgICBjb25zdCByZW5kZXJlciA9IGdldFJlbmRlcmVyKGNlbGwpO1xyXG4gICAgICAgIGlmKHJlbmRlcmVyIGluc3RhbmNlb2YgR3JpZENlbGxSZW5kZXJlckV4KSB7XHJcbiAgICAgICAgICByZW5kZXJlci5vbk1vdXNlVXBFeChjZWxsLCBlTW91c2UsIGFyWFlNb3VzZU9uQ2VsbFVwWzBdLCBhclhZTW91c2VPbkNlbGxVcFsxXSk7XHJcbiAgICAgICAgfVxyXG5cclxuICAgICAgICBpZihhclhZTW91c2VPbkNlbGxVcFswXSA9PT0gYXJYWU1vdXNlT25DZWxsRG93blswXSAmJiBhclhZTW91c2VPbkNlbGxEb3duWzFdID09PSBhclhZTW91c2VPbkNlbGxVcFsxXSkge1xyXG4gICAgICAgICAgaWYocmVuZGVyZXIgaW5zdGFuY2VvZiBHcmlkQ2VsbFJlbmRlcmVyRXgpIHtcclxuICAgICAgICAgICAgcmVuZGVyZXIub25DbGlja0V4KGNlbGwsIGVNb3VzZSwgYXJYWU1vdXNlT25DZWxsVXBbMF0sIGFyWFlNb3VzZU9uQ2VsbFVwWzFdKTtcclxuICAgICAgICAgIH1cclxuICAgICAgICB9XHJcblxyXG4gICAgICAgIG5Sb3dHcmlkRHJhZ2dpbmcgPSAtMTtcclxuICAgICAgICBuWURyYWdnaW5nQW5jaG9yID0gLTE7XHJcbiAgICAgICAgYXJYWU1vdXNlT25DZWxsRG93blswXSA9IC0yO1xyXG4gICAgICAgIGFyWFlNb3VzZU9uQ2VsbERvd25bMV0gPSAtMjtcclxuICAgICAgICBhclhZTW91c2VPbkNlbGxVcFswXSA9IC0xO1xyXG4gICAgICAgIGFyWFlNb3VzZU9uQ2VsbFVwWzFdID0gLTE7XHJcbiAgICAgIH1cclxuICAgICAgLy9lLnByZXZlbnREZWZhdWx0KCk7XHJcbiAgICAgIC8vZS5zdG9wUHJvcGFnYXRpb24oKTtcclxuXHJcbiAgICB9KTtcclxuXHJcbiAgICBsZXQgY2VsbEN1cnJlbnQgOiBERy5HcmlkQ2VsbCB8IG51bGwgPSBudWxsO1xyXG4gICAgdGhpcy5tX2hhbmRsZXJNb3VzZUxlYXZlID0gcnhqcy5mcm9tRXZlbnQoZG9jdW1lbnQsICdtb3VzZWxlYXZlJykuc3Vic2NyaWJlKChlKSA9PiB7XHJcblxyXG4gICAgICBpZihjZWxsQ3VycmVudCAhPT0gbnVsbCkge1xyXG4gICAgICAgIGNvbnN0IHJlbmRlcmVyID0gZ2V0UmVuZGVyZXIoY2VsbEN1cnJlbnQpO1xyXG4gICAgICAgIGlmIChyZW5kZXJlciBpbnN0YW5jZW9mIEdyaWRDZWxsUmVuZGVyZXJFeCkge1xyXG4gICAgICAgICAgY29uc3QgZU1vdXNlID0gZSBhcyBNb3VzZUV2ZW50O1xyXG4gICAgICAgICAgcmVuZGVyZXIub25Nb3VzZUxlYXZlRXgoY2VsbEN1cnJlbnQsIGVNb3VzZSwgLTEsIC0xKTtcclxuICAgICAgICB9XHJcbiAgICAgICAgY2VsbEN1cnJlbnQgPSBudWxsO1xyXG4gICAgICB9XHJcbiAgICB9KTtcclxuXHJcblxyXG4gICAgdGhpcy5tX2hhbmRsZXJNb3VzZU1vdmUgPSByeGpzLmZyb21FdmVudChkb2N1bWVudCwgJ21vdXNlbW92ZScpLnN1YnNjcmliZSgoZSkgPT4ge1xyXG5cclxuICAgICAgaWYoREcudG9EYXJ0KGdyb2suc2hlbGwudikgIT09IERHLnRvRGFydCh2aWV3VGFibGUpKSB7XHJcbiAgICAgICAgcmV0dXJuO1xyXG4gICAgICB9XHJcblxyXG4gICAgICBjb25zdCBiUmVzaXppbmcgPSBuUmVzaXplUm93R3JpZERyYWdnaW5nID49IDA7XHJcbiAgICAgIGlmIChiUmVzaXppbmcpIHtcclxuXHJcbiAgICAgICAgLy9jb25zb2xlLmxvZyhcIkRyYWdnaW5nIDogXCIgKyBoZWFkZXJUaGlzLm1fc3RyQ29sTmFtZSk7XHJcbiAgICAgICAgY29uc3QgZU1vdXNlID0gZSBhcyBNb3VzZUV2ZW50O1xyXG4gICAgICAgIGNvbnN0IG5ZRGlmZiA9IGVNb3VzZS5jbGllbnRZIC0gbllSZXNpemVEcmFnZ2luZ0FuY2hvcjtcclxuICAgICAgICBsZXQgbkhSb3dHcmlkID0gbkhSZXNpemVSb3dzQmVmb3JlRHJhZyArIG5ZRGlmZjtcclxuXHJcbiAgICAgICAgaWYgKG5IUm93R3JpZCA8IFBpbm5lZENvbHVtbi5NSU5fUk9XX0hFSUdIVClcclxuICAgICAgICAgIG5IUm93R3JpZCA9IFBpbm5lZENvbHVtbi5NSU5fUk9XX0hFSUdIVDtcclxuICAgICAgICBlbHNlIGlmIChuSFJvd0dyaWQgPiBQaW5uZWRDb2x1bW4uTUFYX1JPV19IRUlHSFQpXHJcbiAgICAgICAgICBuSFJvd0dyaWQgPSBQaW5uZWRDb2x1bW4uTUFYX1JPV19IRUlHSFQ7XHJcblxyXG4gICAgICAgIGxldCBnID0gZUNhbnZhc1RoaXMuZ2V0Q29udGV4dCgnMmQnKTtcclxuICAgICAgICBpZihnID09PSBudWxsKVxyXG4gICAgICAgICAgcmV0dXJuO1xyXG5cclxuICAgICAgICBnLmZpbGxTdHlsZSA9IFwid2hpdGVcIjtcclxuICAgICAgICBjb25zdCBuSEhlYWRlckNvbHMgPSBHcmlkVXRpbHMuZ2V0R3JpZENvbHVtbkhlYWRlckhlaWdodChncmlkKTtcclxuICAgICAgICBnLmZpbGxSZWN0KDAsbkhIZWFkZXJDb2xzLCBlQ2FudmFzVGhpcy5vZmZzZXRXaWR0aCwgZUNhbnZhc1RoaXMub2Zmc2V0SGVpZ2h0KTtcclxuXHJcbiAgICAgICAgZ3JpZC5zZXRPcHRpb25zKHtcclxuICAgICAgICAgIHJvd0hlaWdodDogbkhSb3dHcmlkIC8vdGhpcyB3b24ndCB0cmlnZ2VyIG9uUm93c1Jleml6ZWQgZXZlbnQsIHdoaWNoIGlzIGEgREcgYnVnXHJcbiAgICAgICAgfSk7XHJcblxyXG4gICAgICAgIG5vdGlmeUFsbFBpbm5lZENvbHNSb3dzUmVzaXplZChoZWFkZXJUaGlzLCBuSFJvd0dyaWQsIHRydWUpO1xyXG4gICAgICAgIG5vdGlmeUFsbENvbHNSb3dzUmVzaXplZChncmlkLCBuSFJvd0dyaWQsIHRydWUpO1xyXG5cclxuICAgICAgICBsZXQgaGVhZGVyID0gbnVsbDtcclxuICAgICAgICBjb25zdCBhciA9IGdyaWQuZGFydC5tX2FyUGlubmVkQ29scztcclxuICAgICAgICBmb3IobGV0IG49MDsgbjxhci5sZW5ndGg7ICsrbikge1xyXG4gICAgICAgICAgaGVhZGVyID0gYXJbbl07XHJcbiAgICAgICAgICBnID0gaGVhZGVyLm1fcm9vdC5nZXRDb250ZXh0KCcyZCcpO1xyXG4gICAgICAgICAgaGVhZGVyLnBhaW50KGcsIGdyaWQpO1xyXG4gICAgICAgIH1cclxuXHJcbiAgICAgICAgdHJ5IHtcclxuICAgICAgICAgIGNvbnN0IGNvbEdyaWQwID0gZ3JpZC5jb2x1bW5zLmJ5SW5kZXgoMCk7XHJcbiAgICAgICAgICBpZiAoY29sR3JpZDAgIT09IG51bGwpXHJcbiAgICAgICAgICAgIGNvbEdyaWQwLnZpc2libGUgPSBmYWxzZTsvL3RlbXBvcmFyeSBhZGRyZXNzZWQgdGhlIERHIGJ1Z1xyXG4gICAgICAgIH1cclxuICAgICAgICBjYXRjaChlKSB7XHJcbiAgICAgICAgICAvL0RHIGJ1Z1xyXG4gICAgICAgIH1cclxuICAgICAgICByZXR1cm47XHJcbiAgICAgIH1cclxuXHJcbiAgICAgIGNvbnN0IGFyWFlPbkNlbGwgPSBbLTEsLTFdO1xyXG4gICAgICBsZXQgblJvd0dyaWQgPSBQaW5uZWRDb2x1bW4uaGl0VGVzdFJvd3MoZUNhbnZhc1RoaXMsIGdyaWQsIGUgYXMgTW91c2VFdmVudCwgZmFsc2UsIGFyWFlPbkNlbGwpO1xyXG4gICAgICBpZihuUm93R3JpZCA+PSAwKSB7XHJcbiAgICAgICAgY29uc3QgY2VsbCA9IGdyaWQuY2VsbChjb2xHcmlkLm5hbWUsIG5Sb3dHcmlkKTtcclxuICAgICAgICBjb25zdCByZW5kZXJlciA9IGdldFJlbmRlcmVyKGNlbGwpO1xyXG5cclxuICAgICAgICBpZiAocmVuZGVyZXIgaW5zdGFuY2VvZiBHcmlkQ2VsbFJlbmRlcmVyRXgpIHtcclxuXHJcbiAgICAgICAgICBpZiAoY2VsbEN1cnJlbnQgPT09IG51bGwpIHtcclxuICAgICAgICAgICAgcmVuZGVyZXIub25Nb3VzZUVudGVyRXgoY2VsbCwgZSBhcyBNb3VzZUV2ZW50LCBhclhZT25DZWxsWzBdLCBhclhZT25DZWxsWzFdKTtcclxuICAgICAgICAgIH1cclxuXHJcbiAgICAgICAgICBpZiAoY2VsbEN1cnJlbnQgIT09IG51bGwgJiYgblJvd0dyaWQgIT09IGNlbGxDdXJyZW50LmdyaWRSb3cpIHtcclxuICAgICAgICAgICAgIHJlbmRlcmVyLm9uTW91c2VMZWF2ZUV4KGNlbGxDdXJyZW50LCBlIGFzIE1vdXNlRXZlbnQsIC0xLCAtMSk7XHJcblxyXG4gICAgICAgICAgIHJlbmRlcmVyLm9uTW91c2VFbnRlckV4KGNlbGwsIGUgYXMgTW91c2VFdmVudCwgYXJYWU9uQ2VsbFswXSwgYXJYWU9uQ2VsbFsxXSk7XHJcbiAgICAgICAgICB9XHJcblxyXG4gICAgICAgICAgcmVuZGVyZXIub25Nb3VzZU1vdmVFeChjZWxsLCBlIGFzIE1vdXNlRXZlbnQsIGFyWFlPbkNlbGxbMF0sIGFyWFlPbkNlbGxbMV0pO1xyXG4gICAgICAgICB9XHJcblxyXG4gICAgICAgIGNlbGxDdXJyZW50ID0gY2VsbDtcclxuICAgICAgfVxyXG4gICAgICBlbHNlIGlmIChjZWxsQ3VycmVudCAhPT0gbnVsbCkge1xyXG4gICAgICAgIGNvbnN0IHJlbmRlcmVyID0gZ2V0UmVuZGVyZXIoY2VsbEN1cnJlbnQpO1xyXG4gICAgICAgIGlmIChyZW5kZXJlciBpbnN0YW5jZW9mIEdyaWRDZWxsUmVuZGVyZXJFeCkge1xyXG4gICAgICAgICAgcmVuZGVyZXIub25Nb3VzZUxlYXZlRXgoY2VsbEN1cnJlbnQsIGUgYXMgTW91c2VFdmVudCwgLTEsIC0xKTtcclxuICAgICAgICB9XHJcblxyXG4gICAgICAgIGNlbGxDdXJyZW50ID0gbnVsbDtcclxuICAgICAgfVxyXG5cclxuICAgICAgblJvd0dyaWQgPSBQaW5uZWRDb2x1bW4uaGl0VGVzdFJvd3MoZUNhbnZhc1RoaXMsIGdyaWQsIGUgYXMgTW91c2VFdmVudCwgdHJ1ZSwgdW5kZWZpbmVkKTtcclxuICAgICAgaWYgKG5Sb3dHcmlkID49IDApIHtcclxuICAgICAgICBuUmVzaXplUm93R3JpZE1vdmluZyA9IG5Sb3dHcmlkO1xyXG4gICAgICAgIGRvY3VtZW50LmJvZHkuc3R5bGUuY3Vyc29yID0gXCJyb3ctcmVzaXplXCI7XHJcbiAgICAgICAgcmV0dXJuO1xyXG4gICAgICB9XHJcblxyXG4gICAgICBpZihuUmVzaXplUm93R3JpZE1vdmluZyA+PSAwKSB7XHJcbiAgICAgICAgblJlc2l6ZVJvd0dyaWRNb3ZpbmcgPSAtMTtcclxuICAgICAgICBkb2N1bWVudC5ib2R5LnN0eWxlLmN1cnNvciA9IFwiYXV0b1wiO1xyXG4gICAgICB9XHJcbiAgICB9KTtcclxuXHJcbiAgICBsZXQgbkNvdW50ID0gMDtcclxuICAgIHRoaXMubV9oYW5kbGVyTW91c2VXaGVlbCA9IHJ4anMuZnJvbUV2ZW50KHRoaXMubV9yb290LCAnd2hlZWwnKS5zdWJzY3JpYmUoKGUpID0+IHtcclxuXHJcbiAgICAgIGlmIChERy50b0RhcnQoZ3Jvay5zaGVsbC52KSAhPT0gREcudG9EYXJ0KHZpZXdUYWJsZSkpIHtcclxuICAgICAgICByZXR1cm47XHJcbiAgICAgIH1cclxuXHJcbiAgICAgIGNvbnN0IGVXaGVlbCA9IGUgYXMgV2hlZWxFdmVudDtcclxuICAgICAgaWYoZVdoZWVsLmRlbHRhWCAhPT0gMCB8fCBlV2hlZWwuZGVsdGFaICE9PSAwKSB7XHJcbiAgICAgICAgcmV0dXJuO1xyXG4gICAgICB9XHJcblxyXG4gICAgICBpZihuQ291bnQgPT09IDEpIHtcclxuICAgICAgICAvL3Njcm9sbCArXHJcbiAgICAgICAgY29uc3QgblJvd0NvdW50ID0gR3JpZFV0aWxzLmdldEdyaWRWaXNpYmxlUm93Q291bnQoZ3JpZCk7XHJcbiAgICAgICAgY29uc3Qgc2Nyb2xsWSA9IGdyaWQudmVydFNjcm9sbDtcclxuICAgICAgICBpZihuUm93Q291bnQgLTEgPiBzY3JvbGxZLm1heCkge1xyXG4gICAgICAgICAgc2Nyb2xsWS5zZXRWYWx1ZXMoc2Nyb2xsWS5taW5SYW5nZSwgc2Nyb2xsWS5tYXhSYW5nZSwgc2Nyb2xsWS5taW4gKyAxLCBzY3JvbGxZLm1heCArIDEpO1xyXG4gICAgICAgIH1cclxuICAgICAgICBuQ291bnQgPSAwO1xyXG4gICAgICB9XHJcbiAgICAgIGVsc2UgaWYobkNvdW50ID09PSAtMSlcclxuICAgICAge1xyXG4gICAgICAgIC8vc2Nyb2xsIC1cclxuICAgICAgICBjb25zdCBzY3JvbGxZID0gZ3JpZC52ZXJ0U2Nyb2xsO1xyXG4gICAgICAgIGlmKHNjcm9sbFkubWluID49MSkge1xyXG4gICAgICAgICAgc2Nyb2xsWS5zZXRWYWx1ZXMoc2Nyb2xsWS5taW5SYW5nZSwgc2Nyb2xsWS5tYXhSYW5nZSwgc2Nyb2xsWS5taW4gLSAxLCBzY3JvbGxZLm1heCAtIDEpO1xyXG4gICAgICAgIH1cclxuICAgICAgICBuQ291bnQgPSAwO1xyXG4gICAgICB9XHJcbiAgICAgICBlbHNlIHtcclxuICAgICAgICBuQ291bnQgPSBlV2hlZWwuZGVsdGFZID4gMCA/IDEgOiAtMTtcclxuICAgICAgfVxyXG5cclxuXHJcblxyXG4gICAgICAvL2NvbnNvbGUubG9nKGVXaGVlbC5kZWx0YVggKyBcIiBcIiArIGVXaGVlbC5kZWx0YVkpO1xyXG4gICAgfSk7XHJcbiAgfVxyXG5cclxuICBpc1Bpbm5lZCgpIDogYm9vbGVhbiB7XHJcbiAgICByZXR1cm4gdGhpcy5tX2NvbEdyaWQgIT09IG51bGw7XHJcbiAgfVxyXG5cclxuICBnZXRHcmlkQ29sdW1uKCkgOiBERy5HcmlkQ29sdW1uIHwgbnVsbHtcclxuICAgIHJldHVybiB0aGlzLm1fY29sR3JpZDtcclxuICB9XHJcblxyXG4gIGdldFdpZHRoKCkgOiBudW1iZXIge1xyXG4gICAgcmV0dXJuIHRoaXMubV9yb290ID09PSBudWxsID8gLTEgOiB0aGlzLm1fcm9vdC5vZmZzZXRXaWR0aDtcclxuICB9XHJcblxyXG4gIGdldFJvb3QoKSA6IEhUTUxDYW52YXNFbGVtZW50IHwgbnVsbCB7XHJcbiAgICByZXR1cm4gdGhpcy5tX3Jvb3Q7XHJcbiAgfVxyXG5cclxuICBwdWJsaWMgY2xvc2UoKSA6IHZvaWQge1xyXG5cclxuICAgIGlmKHRoaXMubV9jb2xHcmlkID09PSBudWxsKSB7XHJcbiAgICAgIHRocm93IG5ldyBFcnJvcihcIkNvbHVtbiBoYXMgYWxyZWFkeSBiZWVuIHVucGlubmVkXCIpO1xyXG4gICAgfVxyXG5cclxuICAgIGlmKHRoaXMubV9vYnNlcnZlclJlc2l6ZUdyaWQgIT09IG51bGwpIHtcclxuICAgICAgdGhpcy5tX29ic2VydmVyUmVzaXplR3JpZC5kaXNjb25uZWN0KCk7XHJcbiAgICAgIHRoaXMubV9vYnNlcnZlclJlc2l6ZUdyaWQgPSBudWxsO1xyXG4gICAgfVxyXG4vKm15IGNoYW5nZXNcclxuICAgIGlmKHRoaXMubV9vYnNlcnZlclJlc2l6ZSAhPT0gbnVsbCkge1xyXG4gICAgICB0aGlzLm1fb2JzZXJ2ZXJSZXNpemUuZGlzY29ubmVjdCgpO1xyXG4gICAgICB0aGlzLm1fb2JzZXJ2ZXJSZXNpemUgPSBudWxsO1xyXG4gICAgfVxyXG4gICAgKi9cclxuICAgIHRoaXMubV9oYW5kbGVyVlNjcm9sbC51bnN1YnNjcmliZSgpO1xyXG4gICAgdGhpcy5tX2hhbmRsZXJWU2Nyb2xsID0gbnVsbDtcclxuXHJcbiAgICB0aGlzLm1faGFuZGxlclJvd3NSZXNpemVkLnVuc3Vic2NyaWJlKCk7XHJcbiAgICB0aGlzLm1faGFuZGxlclJvd3NSZXNpemVkID0gbnVsbDtcclxuXHJcbiAgICB0aGlzLm1faGFuZGxlclJvd3NTb3J0ZWQudW5zdWJzY3JpYmUoKTtcclxuICAgIHRoaXMubV9oYW5kbGVyUm93c1NvcnRlZCA9IG51bGw7XHJcblxyXG4gICAgdGhpcy5tX2hhbmRsZXJSb3dzRmlsdGVyaW5nLnVuc3Vic2NyaWJlKCk7XHJcbiAgICB0aGlzLm1faGFuZGxlclJvd3NGaWx0ZXJpbmcgPSBudWxsO1xyXG5cclxuICAgIHRoaXMubV9oYW5kbGVyQ3VyclJvdy51bnN1YnNjcmliZSgpO1xyXG4gICAgdGhpcy5tX2hhbmRsZXJDdXJyUm93ID0gbnVsbDtcclxuXHJcbiAgICB0aGlzLm1faGFuZGxlclNlbC51bnN1YnNjcmliZSgpO1xyXG4gICAgdGhpcy5tX2hhbmRsZXJTZWwgPSBudWxsO1xyXG5cclxuICAgIHRoaXMubV9oYW5kbGVyTW91c2VEb3duLnVuc3Vic2NyaWJlKCk7XHJcbiAgICB0aGlzLm1faGFuZGxlck1vdXNlRG93biA9IG51bGw7XHJcblxyXG4gICAgdGhpcy5tX2hhbmRsZXJNb3VzZVVwLnVuc3Vic2NyaWJlKCk7XHJcbiAgICB0aGlzLm1faGFuZGxlck1vdXNlVXAgPSBudWxsO1xyXG5cclxuICAgIHRoaXMubV9oYW5kbGVyTW91c2VMZWF2ZS51bnN1YnNjcmliZSgpO1xyXG4gICAgdGhpcy5tX2hhbmRsZXJNb3VzZUxlYXZlID0gbnVsbDtcclxuXHJcbiAgICB0aGlzLm1faGFuZGxlck1vdXNlTW92ZS51bnN1YnNjcmliZSgpO1xyXG4gICAgdGhpcy5tX2hhbmRsZXJNb3VzZU1vdmUgPSBudWxsO1xyXG5cclxuICAgIHRoaXMubV9oYW5kbGVyTW91c2VXaGVlbC51bnN1YnNjcmliZSgpO1xyXG4gICAgdGhpcy5tX2hhbmRsZXJNb3VzZVdoZWVsID0gbnVsbDtcclxuXHJcbiAgICBjb25zdCBncmlkID0gZ2V0R3JpZCh0aGlzLm1fY29sR3JpZCk7XHJcbiAgICBpZihncmlkID09PSBudWxsKXtcclxuICAgICAgdGhyb3cgbmV3IEVycm9yKFwiQ29sdW1uICdcIiArIHRoaXMubV9jb2xHcmlkLm5hbWUgKyBcIicgaXMgZGlzY29ubmVjdGVkIGZyb20gZ3JpZC5cIik7XHJcbiAgICB9XHJcblxyXG4gICAgY29uc3QgZGFydCA9IERHLnRvRGFydChncmlkKTtcclxuICAgIGNvbnN0IGFyID0gZGFydC5tX2FyUGlubmVkQ29scztcclxuICAgIGNvbnN0IG5JZHggPSBhci5pbmRleE9mKHRoaXMpO1xyXG4gICAgYXIuc3BsaWNlKG5JZHgsIDEpO1xyXG5cclxuICAgIGlmKHRoaXMubV9yb290ID09PSBudWxsKVxyXG4gICAgICB0aHJvdyBuZXcgRXJyb3IoJ1Jvb3QgY2Fubm90IGJlIG51bGwnKTtcclxuXHJcbiAgICBsZXQgbklkeFBpbm5lZCA9IC0xO1xyXG4gICAgbGV0IGNvbEdyaWRUbXA9IG51bGw7XHJcbiAgICBmb3IobGV0IG49bklkeDsgbjxhci5sZW5ndGg7ICsrbikge1xyXG4gICAgICBjb2xHcmlkVG1wID0gYXJbbl07XHJcbiAgICAgIGNvbEdyaWRUbXAubV9yb290LnN0eWxlLmxlZnQgPSAoY29sR3JpZFRtcC5tX3Jvb3Qub2Zmc2V0TGVmdCAtIHRoaXMubV9yb290Lm9mZnNldFdpZHRoKS50b1N0cmluZygpICsgXCJweFwiO1xyXG5cclxuICAgICAgbklkeFBpbm5lZCA9ICBjb2xHcmlkVG1wLm1fY29sR3JpZC5zZXR0aW5ncy5pZHhQaW5uZWQ7XHJcbiAgICAgIGNvbEdyaWRUbXAubV9jb2xHcmlkLnNldHRpbmdzLmlkeFBpbm5lZCA9IG47XHJcbiAgICB9XHJcblxyXG4gICAgaWYoIUdyaWRVdGlscy5pc1Jvd0hlYWRlcih0aGlzLm1fY29sR3JpZCkpIHtcclxuICAgICAgdGhpcy5tX2NvbEdyaWQuc2V0dGluZ3MuaWR4UGlubmVkID0gLTE7XHJcbiAgICAgIHRoaXMubV9jb2xHcmlkLnNldHRpbmdzLmlzUGlubmVkID0gZmFsc2U7XHJcbiAgICB9XHJcblxyXG5cclxuICAgIGlmKHRoaXMubV9uV2lkdGhCdWcgPj0gMCkge1xyXG4gICAgICB0cnkge1xyXG4gICAgICAgIHRoaXMubV9jb2xHcmlkLndpZHRoID0gdGhpcy5tX25XaWR0aEJ1ZztcclxuICAgICAgfVxyXG4gICAgICBjYXRjaChlKSB7XHJcbiAgICAgICAgLy9ERyBidWdcclxuICAgICAgICBjb25zb2xlLmVycm9yKFwiRVJST1I6IENvdWxkbid0IHNldCB0aGUgd2lkdGggdG8gXCIgKyB0aGlzLm1fbldpZHRoQnVnICsgXCIgZm9yIGNvbHVtbiAnXCIgKyB0aGlzLm1fY29sR3JpZC5uYW1lICsgXCInIGR1ZSB0byBhIERHIGJ1Zy4gVGhpcyBjb3VsZCBiZSBpZ25vcmVkIGlmIHRoZSBjb2x1bW4gdmlzdWFsbHkgbG9va3Mgb2suXCIpO1xyXG4gICAgICB9XHJcbiAgICB9XHJcblxyXG4gICAgdHJ5IHtcclxuICAgICAgdGhpcy5tX2NvbEdyaWQudmlzaWJsZSA9IHRydWU7XHJcbiAgICB9XHJcbiAgICBjYXRjaChlKSB7XHJcbiAgICAgIC8vREcgYnVnXHJcbiAgICAgIGNvbnNvbGUuZXJyb3IoXCJFUlJPUjogQ291bGRuJ3Qgc2hvdyBjb2x1bW4gJ1wiICsgdGhpcy5tX2NvbEdyaWQubmFtZSArIFwiJyBkdWUgdG8gYSBERyBidWcuIFRoaXMgY291bGQgYmUgaWdub3JlZCBpZiB0aGUgY29sdW1uIHZpc3VhbGx5IGxvb2tzIG9rLlwiKTtcclxuICAgIH1cclxuXHJcblxyXG5cclxuICAgIGdyaWQuY2FudmFzLnN0eWxlLmxlZnQgPSAoZ3JpZC5jYW52YXMub2Zmc2V0TGVmdCAtIHRoaXMubV9yb290Lm9mZnNldFdpZHRoKS50b1N0cmluZygpICsgXCJweFwiO1xyXG4gICAgZ3JpZC5vdmVybGF5LnN0eWxlLmxlZnQ9IChncmlkLm92ZXJsYXkub2Zmc2V0TGVmdCAtIHRoaXMubV9yb290Lm9mZnNldFdpZHRoKS50b1N0cmluZygpICsgXCJweFwiO1xyXG4gICAgZ3JpZC5jYW52YXMuc3R5bGUud2lkdGggPSAoZ3JpZC5jYW52YXMub2Zmc2V0V2lkdGggKyB0aGlzLm1fcm9vdC5vZmZzZXRXaWR0aCkudG9TdHJpbmcoKSArIFwicHhcIjtcclxuICAgIGdyaWQub3ZlcmxheS5zdHlsZS53aWR0aD0gKGdyaWQub3ZlcmxheS5vZmZzZXRXaWR0aCArIHRoaXMubV9yb290Lm9mZnNldFdpZHRoKS50b1N0cmluZygpICsgXCJweFwiO1xyXG5cclxuICAgIGlmKHRoaXMubV9yb290LnBhcmVudE5vZGUgIT09IG51bGwpXHJcbiAgICAgdGhpcy5tX3Jvb3QucGFyZW50Tm9kZS5yZW1vdmVDaGlsZCh0aGlzLm1fcm9vdCk7XHJcblxyXG4gICAgdGhpcy5tX3Jvb3QgPSBudWxsO1xyXG5cclxuICAgIGlmIChkYXJ0Lm1fYXJQaW5uZWRDb2xzLmxlbmd0aCA9PT0gMSAmJiBkYXJ0Lm1fYXJQaW5uZWRDb2xzWzBdLm1fY29sR3JpZC5pZHggPT09IDAgJiYgdGhpcy5tX2NvbEdyaWQuaWR4ICE9PSAwKSB7XHJcblxyXG4gICAgICAgIC8vIHRyeXtjb2xHcmlkMC52aXNpYmxlID0gdHJ1ZTt9XHJcbiAgICAgICAgdHJ5IHtcclxuICAgICAgICAgIGRhcnQubV9hclBpbm5lZENvbHNbMF0uY2xvc2UoKTtcclxuICAgICAgICB9IGNhdGNoIChlKSB7XHJcbiAgICAgICAgICBjb25zb2xlLmVycm9yKFwiRVJST1I6IENvdWxkbid0IGNsb3NlIHBpbm5lZCBjb2x1bW4gJ1wiICsgZGFydC5tX2FyUGlubmVkQ29sc1swXS5tX2NvbEdyaWQubmFtZSArIFwiJyBcIik7XHJcbiAgICAgICAgfVxyXG4gICAgfVxyXG4gICAgdGhpcy5tX2NvbEdyaWQgPSBudWxsO1xyXG4gIH1cclxuXHJcblxyXG4gIHByaXZhdGUgcGFpbnQoZyA6IENhbnZhc1JlbmRlcmluZ0NvbnRleHQyRCB8IG51bGwsIGdyaWQgOiBERy5HcmlkKSA6IHZvaWQge1xyXG4gICAgLy9jb25zdCBuV0RpdiA9IGVudHJ5LmNvbnRlbnRCb3hTaXplID8gZW50cnkuY29udGVudEJveFNpemVbMF0uaW5saW5lU2l6ZSA6IGVudHJ5LmNvbnRlbnRSZWN0LndpZHRoO1xyXG5cclxuICAgIGlmKGcgPT09IG51bGwpIHtcclxuICAgICAgcmV0dXJuO1xyXG4gICAgfVxyXG5cclxuICAgIGlmKHRoaXMubV9yb290ID09PSBudWxsKSB7XHJcbiAgICAgIHRocm93IG5ldyBFcnJvcignUm9vdCBjYW5ub3QgYmUgbnVsbC4nKTtcclxuICAgIH1cclxuXHJcbiAgICBpZih0aGlzLm1fY29sR3JpZCA9PT0gbnVsbCkge1xyXG4gICAgICB0aHJvdyBuZXcgRXJyb3IoJ0NvbHVtbiBncmlkIGNhbm5vdCBiZSBudWxsLicpO1xyXG4gICAgfVxyXG4gICAgY29uc3QgZGZyYW1lID0gZ3JpZC5kYXRhRnJhbWU7XHJcbiAgICBjb25zdCBuVyA9IHRoaXMubV9yb290Lm9mZnNldFdpZHRoO1xyXG4gICAgY29uc3QgbkggPSB0aGlzLm1fcm9vdC5vZmZzZXRIZWlnaHQ7XHJcblxyXG4gICAgZy5maWxsU3R5bGUgPSBcIndoaXRlXCI7XHJcbiAgICBnLmZpbGxSZWN0KDAsMCwgblcqd2luZG93LmRldmljZVBpeGVsUmF0aW8sIG5IKndpbmRvdy5kZXZpY2VQaXhlbFJhdGlvKTtcclxuXHJcbiAgICBpZih0aGlzLm1fY29sR3JpZC5uYW1lID09PSBudWxsKVxyXG4gICAgICByZXR1cm47XHJcblxyXG4gICAgY29uc3QgYml0c2V0RmlsdGVyID0gZGZyYW1lLmZpbHRlcjtcclxuICAgIGlmKGJpdHNldEZpbHRlci5mYWxzZUNvdW50ID09PSBkZnJhbWUucm93Q291bnQpXHJcbiAgICAgIHJldHVybjsgLy9ldmVyeXRoaW5nIGlzIGZpbHRlcmVkXHJcblxyXG4gICAgLy9jb2x1bW4gSGVhZGVyXHJcbiAgICBjb25zdCBvcHRpb25zIDogYW55ID0gZ3JpZC5nZXRPcHRpb25zKHRydWUpO1xyXG5cclxuICAgIGNvbnN0IGZvbnRDZWxsRGVmYXVsdCA9IG9wdGlvbnMubG9vay5kZWZhdWx0Q2VsbEZvbnQ7XHJcblxyXG4gICAgbGV0IGZvbnQgPSBvcHRpb25zLmxvb2suY29sSGVhZGVyRm9udCA9PSBudWxsIHx8IG9wdGlvbnMubG9vay5jb2xIZWFkZXJGb250ID09PSB1bmRlZmluZWQgPyBcImJvbGQgMTRweCBWb2x0YSBUZXh0LCBBcmlhbFwiIDogb3B0aW9ucy5sb29rLmNvbEhlYWRlckZvbnQ7XHJcbiAgICBsZXQgZm9udFNjYWxlZCA9IEdyaWRVdGlscy5zY2FsZUZvbnQoZm9udCwgd2luZG93LmRldmljZVBpeGVsUmF0aW8pO1xyXG4gICAgZy5mb250ID0gZm9udFNjYWxlZDtcclxuXHJcbiAgICBsZXQgc3RyID0gVGV4dFV0aWxzLnRyaW1UZXh0KHRoaXMubV9jb2xHcmlkLm5hbWUsIGcsIG5XKTtcclxuXHJcbiAgICBjb25zdCB0bSA9IGcubWVhc3VyZVRleHQoc3RyKTtcclxuICAgIGNvbnN0IG5XTGFiZWwgPSB0bS53aWR0aDtcclxuXHJcbiAgICBjb25zdCBuQXNjZW50ID0gTWF0aC5hYnModG0uYWN0dWFsQm91bmRpbmdCb3hBc2NlbnQpO1xyXG4gICAgY29uc3QgbkRlc2NlbnQgPSB0bS5hY3R1YWxCb3VuZGluZ0JveERlc2NlbnQ7XHJcbiAgICBjb25zdCBuSEZvbnQgPSAgbkFzY2VudCArIG5EZXNjZW50Oy8vICsgMipuWUluc2V0O1xyXG5cclxuICAgIC8vbGV0IGNlbGxDSCA9IGdyaWQuY2VsbCh0aGlzLm1fY29sR3JpZC5uYW1lLCAtMSk7XHJcbiAgICAvL2xldCByZW5kZXJlciA9IGNlbGxDSC5yZW5kZXJlcjtcclxuXHJcbiAgICBsZXQgblggPSAwO1xyXG4gICAgbGV0IG5ZID0gMDtcclxuICAgIGNvbnN0IG5IQ0ggPSBHcmlkVXRpbHMuZ2V0R3JpZENvbHVtbkhlYWRlckhlaWdodChncmlkKSp3aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbztcclxuICAgIGcudGV4dEFsaWduID0gJ3N0YXJ0JztcclxuICAgIGcuZmlsbFN0eWxlID0gXCJCbGFja1wiO1xyXG4gICAgbGV0IG5ZT2Zmc2V0ID0gTWF0aC5mbG9vcigobkhDSCAtIG5IRm9udCkvMik7XHJcbiAgICBjb25zdCBuWFggPSBuWCArICgoblcqd2luZG93LmRldmljZVBpeGVsUmF0aW8gLSBuV0xhYmVsKSA+PiAxKTtcclxuICAgIGxldCBuWVkgPSAoblkgKyBuSENIIC0gTWF0aC5jZWlsKDMqd2luZG93LmRldmljZVBpeGVsUmF0aW8pKTsvLy0yKndpbmRvdy5kZXZpY2VQaXhlbFJhdGlvKTtcclxuICAgIC8vb25zb2xlLmxvZyhcIm5YWCBcIiArIG5YWCArIFwiIG5ZWSA9IFwiICsgbllZICsgXCIgQ0hIIFwiICsgbkhDSCk7XHJcbiAgICBnLmZpbGxUZXh0KHN0ciwgblhYLCBuWVkpO1xyXG5cclxuICAgIC8vaWYob3B0aW9ucy5sb29rLnNob3dSb3dHcmlkbGluZXMpIHtcclxuXHJcblxyXG4gICAgLy99XHJcblxyXG5cclxuXHJcbiAgICAvL1JlZ3VsYXIgY2VsbHNcclxuICAgIGNvbnN0IG5Sb3dDdXJyZW50ID0gIGRmcmFtZS5jdXJyZW50Um93LmlkeDtcclxuICAgIGNvbnN0IGJpdHNldFNlbCA9IGRmcmFtZS5zZWxlY3Rpb247XHJcblxyXG4gICAgY29uc3QgYXJSb3dzTWluTWF4ID0gWy0xLC0xXTtcclxuICAgIEdyaWRVdGlscy5maWxsVmlzaWJsZVZpZXdwb3J0Um93cyhhclJvd3NNaW5NYXgsIGdyaWQpO1xyXG4gICAgY29uc3QgblJvd01pbiA9IGFyUm93c01pbk1heFswXTtcclxuICAgIGNvbnN0IG5Sb3dNYXggPSBhclJvd3NNaW5NYXhbMV07XHJcblxyXG4gICAgLy9jb25zb2xlLmxvZyhuUm93TWluICsgXCIgXCIgKyBuUm93TWF4KTtcclxuICAgIGNvbnN0IG5IUm93ID0gR3JpZFV0aWxzLmdldEdyaWRSb3dIZWlnaHQoZ3JpZCk7XHJcbiAgICBuWU9mZnNldCA9IG5IQ0g7XHJcbiAgICBjb25zdCBuSFJvd0dyaWQgPSBuSFJvdyp3aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbztcclxuICAgIGxldCBjZWxsUkggPSBudWxsO1xyXG5cclxuICAgIGxldCBuV1cgPSBuVyp3aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbztcclxuICAgIC8vY29uc3QgbkhIID0gbkhSb3dHcmlkO1xyXG5cclxuICAgIGNvbnN0IGFyVGFibGVSb3dzID0gbmV3IEFycmF5KG5Sb3dNYXggLSBuUm93TWluICsxKTtcclxuICAgIGxldCBuUm93VGFibGUgPSAtMTtcclxuICAgIGxldCBiU2VsID0gZmFsc2U7XHJcbiAgICBmb3IobGV0IG5SRz1uUm93TWluOyBuUkc8PW5Sb3dNYXg7ICsrblJHKSB7XHJcbiAgICAgIHRyeSB7XHJcbiAgICAgICAgY2VsbFJIID0gZ3JpZC5jZWxsKHRoaXMubV9jb2xHcmlkLm5hbWUsIG5SRyk7XHJcbiAgICAgIH0gY2F0Y2ggKGUpIC8vdG8gYWRkcmVzcyBERyBidWcgd2hlbiBldmVyeXRoaW5nIGlzIGZpbHRlcmVkXHJcbiAgICAgIHtcclxuICAgICAgICBjb250aW51ZTtcclxuICAgICAgfVxyXG5cclxuICAgICAgaWYgKGNlbGxSSC50YWJsZVJvd0luZGV4ID09PSB1bmRlZmluZWQpLy9ERyBidWdcclxuICAgICAgICBjb250aW51ZTtcclxuXHJcbiAgICAgIG5Sb3dUYWJsZSA9IGNlbGxSSC50YWJsZVJvd0luZGV4ID09PSBudWxsID8gLTEgOiBjZWxsUkgudGFibGVSb3dJbmRleDtcclxuICAgICAgYXJUYWJsZVJvd3NbblJHIC0gblJvd01pbl0gPSBuUm93VGFibGU7XHJcblxyXG4gICAgICBuWVkgPSBuWU9mZnNldCArIChuUkcgLSBuUm93TWluKSAqIG5IUm93R3JpZDtcclxuXHJcbiAgICAgIGxldCByZW5kZXJlcjogYW55ID0gR3JpZFV0aWxzLmdldEdyaWRDb2x1bW5SZW5kZXJlcihjZWxsUkguZ3JpZENvbHVtbik7XHJcbiAgICAgIGlmIChyZW5kZXJlciA9PT0gbnVsbCkge1xyXG4gICAgICAgIHRyeSB7XHJcbiAgICAgICAgICByZW5kZXJlciA9IGNlbGxSSC5yZW5kZXJlcjtcclxuICAgICAgICB9IGNhdGNoIChlKSB7XHJcbiAgICAgICAgICBjb25zb2xlLmVycm9yKFwiQ291bGQgbm90IG9idGFpbiByZW5kZXJlciBmb3IgREcgY2VsbC4gREcgYnVnIFwiICsgdGhpcy5tX2NvbEdyaWQubmFtZSArIFwiIHJvdyBcIiArIG5SRyk7XHJcbiAgICAgICAgICBjb250aW51ZTtcclxuICAgICAgICB9XHJcbiAgICAgIH1cclxuXHJcbiAgICAgIGlmIChyZW5kZXJlciA9PT0gbnVsbCB8fCByZW5kZXJlciA9PT0gdW5kZWZpbmVkKSB7XHJcbiAgICAgICAgY29uc29sZS5lcnJvcihcIkNvdWxkbid0IGZpbmQgcmVuZGVyZXIgZm9yIHBpbm5lZCBjb2x1bW4gXCIgKyB0aGlzLm1fY29sR3JpZC5uYW1lICsgXCIgcm93IFwiICsgblJHKTtcclxuICAgICAgICBjb250aW51ZTtcclxuICAgICAgfVxyXG5cclxuICAgICAgLy9sZXQgbllZID0gblk7Ly8qd2luZG93LmRldmljZVBpeGVsUmF0aW87XHJcblxyXG5cclxuICAgICAgZm9udCA9IGNlbGxSSC5zdHlsZS5mb250O1xyXG4gICAgICBmb250U2NhbGVkID0gR3JpZFV0aWxzLnNjYWxlRm9udChmb250LCB3aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbyk7XHJcbiAgICAgIGlmIChmb250U2NhbGVkICE9PSBudWxsKSB7XHJcbiAgICAgICAgY2VsbFJILnN0eWxlLmZvbnQgPSBmb250U2NhbGVkO1xyXG4gICAgICB9XHJcblxyXG4gICAgICBpZiAoblcgPiAwICYmIG5IUm93R3JpZCA+IDApIHsgLy90byBhZGRyZXNzIGEgYnVnIGNhdXNlZCBlaXRoZXIgREcgb3IgY2xpZW50IGFwcFxyXG4gICAgICAgIHRyeSB7XHJcbiAgICAgICAgICBpZiAocmVuZGVyZXIubmFtZSA9PT0gJ01vbGVjdWxlJykge1xyXG4gICAgICAgICAgICByZW5kZXJlci5yZW5kZXIoZywgMCwgbllZL3dpbmRvdy5kZXZpY2VQaXhlbFJhdGlvLCBuV1cvd2luZG93LmRldmljZVBpeGVsUmF0aW8sIG5IUm93R3JpZC93aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbywgY2VsbFJILCBjZWxsUkguc3R5bGUpO1xyXG4gICAgICAgICAgfVxyXG4gICAgICAgICAgZWxzZSByZW5kZXJlci5yZW5kZXIoZywgMCwgbllZLCBuV1csIG5IUm93R3JpZCwgY2VsbFJILCBjZWxsUkguc3R5bGUpO1xyXG5cclxuICAgICAgICB9IGNhdGNoIChlKSB7XHJcbiAgICAgICAgICBjb25zb2xlLmVycm9yKFwiQ291bGQgbm90IHBhaW50IGNlbGwgZm9yIHBpbm5lZCBjb2x1bW4gXCIgKyB0aGlzLm1fY29sR3JpZC5uYW1lICsgXCIgcm93IFwiICsgblJHKTtcclxuICAgICAgICAgIGNvbnRpbnVlO1xyXG4gICAgICAgICAgLy90aHJvdyBlO1xyXG4gICAgICAgIH1cclxuICAgICAgfVxyXG4gICAgfVxyXG5cclxuXHJcbiAgICAvL1BhaW50IEdyaWRcclxuICAgIGcuc3Ryb2tlU3R5bGUgPSBcIkdhaW5zYm9yb1wiO1xyXG4gICAgZy5iZWdpblBhdGgoKTtcclxuICAgIGcubW92ZVRvKDAsIG5ZKndpbmRvdy5kZXZpY2VQaXhlbFJhdGlvKTtcclxuICAgIGcubGluZVRvKDAsIChuWSArIG5IQ0gtMSp3aW5kb3cuZGV2aWNlUGl4ZWxSYXRpbykpO1xyXG4gICAgZy5zdHJva2UoKTtcclxuXHJcbiAgICBnLmJlZ2luUGF0aCgpO1xyXG4gICAgZy5tb3ZlVG8oMCwgbllPZmZzZXQgKyAxKTtcclxuICAgIGcubGluZVRvKG5XVywgbllPZmZzZXQgKyAxKTtcclxuICAgIGcuc3Ryb2tlKCk7XHJcblxyXG4gICAgZm9yKGxldCBuUkc9blJvd01pbjsgblJHPD1uUm93TWF4OyArK25SRylcclxuICAgIHtcclxuICAgICAgbllZID0gbllPZmZzZXQgKyAoblJHIC0gblJvd01pbikgKiBuSFJvd0dyaWQ7XHJcblxyXG4gICAgICAvL2lmKG9wdGlvbnMubG9vay5zaG93Um93R3JpZGxpbmVzKSB7XHJcblxyXG4gICAgICAgIGcuYmVnaW5QYXRoKCk7XHJcbiAgICAgICAgZy5tb3ZlVG8oMCwgbllZICsgbkhSb3dHcmlkKzEpO1xyXG4gICAgICAgIGcubGluZVRvKG5XVywgbllZICsgbkhSb3dHcmlkKzEpO1xyXG4gICAgICAgIGcuc3Ryb2tlKCk7XHJcblxyXG4gICAgICAgIGcuYmVnaW5QYXRoKCk7XHJcbiAgICAgICAgZy5tb3ZlVG8oMCwgbllZKTtcclxuICAgICAgICBnLmxpbmVUbygwLCBuWVkgKyBuSFJvd0dyaWQrMSk7XHJcbiAgICAgICAgZy5zdHJva2UoKTtcclxuICAgICAgLy99XHJcbiAgICAgIG5Sb3dUYWJsZSA9IGFyVGFibGVSb3dzW25SRyAtIG5Sb3dNaW5dO1xyXG4gICAgICBiU2VsID0gblJvd1RhYmxlIDwgMCA/IGZhbHNlIDogYml0c2V0U2VsLmdldChuUm93VGFibGUpO1xyXG4gICAgICBpZihiU2VsKVxyXG4gICAgICB7XHJcbiAgICAgICAgZy5nbG9iYWxBbHBoYSA9IDAuMjtcclxuICAgICAgICBnLmZpbGxTdHlsZSA9IFBpbm5lZENvbHVtbi5TRUxFQ1RJT05fQ09MT1I7XHJcbiAgICAgICAgZy5maWxsUmVjdCgwLCBuWVksIG5XVywgbkhSb3dHcmlkKTtcclxuICAgICAgICBnLmdsb2JhbEFscGhhID0gMTtcclxuICAgICAgfVxyXG5cclxuICAgICAgaWYoblJvd0N1cnJlbnQgPT09IG5Sb3dUYWJsZSlcclxuICAgICAge1xyXG4gICAgICAgIGcuZ2xvYmFsQWxwaGEgPSAwLjI7XHJcbiAgICAgICAgZy5maWxsU3R5bGUgPSBQaW5uZWRDb2x1bW4uQUNUSVZFX0NFTExfQ09MT1I7XHJcbiAgICAgICAgZy5maWxsUmVjdCgwLCBuWVksIG5XVywgbkhSb3dHcmlkKTtcclxuICAgICAgICBnLmdsb2JhbEFscGhhID0gMTtcclxuICAgICAgfVxyXG4gICAgfVxyXG4gIH1cclxuXHJcblxyXG4gIHByaXZhdGUgc3RhdGljIGhpdFRlc3RSb3dzKGVDYW52YXNQaW5uZWQgOiBIVE1MQ2FudmFzRWxlbWVudCwgZ3JpZCA6IERHLkdyaWQsIGUgOiBNb3VzZUV2ZW50LCBiQm9yZGVyIDogYm9vbGVhbiwgYXJYWU9uQ2VsbCA6IEFycmF5PG51bWJlcj4gfCB1bmRlZmluZWQpXHJcbiAge1xyXG4gICAgY29uc3QgcmVjdCA9IGVDYW52YXNQaW5uZWQuZ2V0Qm91bmRpbmdDbGllbnRSZWN0KCk7XHJcbiAgICBjb25zdCBzY3JvbGxMZWZ0PSB3aW5kb3cucGFnZVhPZmZzZXQgfHwgZG9jdW1lbnQuZG9jdW1lbnRFbGVtZW50LnNjcm9sbExlZnQ7XHJcbiAgICBjb25zdCBzY3JvbGxUb3AgPSB3aW5kb3cucGFnZVlPZmZzZXQgfHwgZG9jdW1lbnQuZG9jdW1lbnRFbGVtZW50LnNjcm9sbFRvcDtcclxuICAgIGNvbnN0IG5ZID0gcmVjdC50b3AgICsgc2Nyb2xsVG9wO1xyXG4gICAgY29uc3QgblggPSByZWN0LmxlZnQgKyBzY3JvbGxMZWZ0O1xyXG5cclxuICAgIGlmKG5YIDw9IGUuY2xpZW50WCAmJiBlLmNsaWVudFggPD0gblggKyBlQ2FudmFzUGlubmVkLm9mZnNldFdpZHRoKSAgIC8vb24gdGhlIHJvd3MgaGVhZGVyXHJcbiAgICB7XHJcbiAgICAgIGNvbnN0IG5ISGVhZGVyQ29scyA9IEdyaWRVdGlscy5nZXRHcmlkQ29sdW1uSGVhZGVySGVpZ2h0KGdyaWQpO1xyXG4gICAgICBjb25zdCBuSFJvd0dyaWQgPSBHcmlkVXRpbHMuZ2V0R3JpZFJvd0hlaWdodChncmlkKTtcclxuXHJcbiAgICAgIGNvbnN0IGFyTWluTWF4Um93cyA9IFstMSwtMV07XHJcbiAgICAgIEdyaWRVdGlscy5maWxsVmlzaWJsZVZpZXdwb3J0Um93cyhhck1pbk1heFJvd3MsIGdyaWQpO1xyXG4gICAgICBjb25zdCBuUm93TWluID0gYXJNaW5NYXhSb3dzWzBdO1xyXG4gICAgICBjb25zdCBuUm93TWF4ID0gYXJNaW5NYXhSb3dzWzFdO1xyXG5cclxuICAgICAgY29uc3QgbllNb3VzZU9uSGVhZGVyID0gZS5jbGllbnRZIC0gblk7XHJcblxyXG4gICAgICBsZXQgbllCb3JkZXIgPSAtMTtcclxuICAgICAgbGV0IG5ZRGlmZiA9IC0xO1xyXG5cclxuICAgICAgZm9yKGxldCBuUm93PW5Sb3dNaW47IG5Sb3c8PSBuUm93TWF4OyArK25Sb3cpXHJcbiAgICAgIHtcclxuICAgICAgICBuWUJvcmRlciA9IG5ISGVhZGVyQ29scyArIChuUm93IC0gblJvd01pbisxKSpuSFJvd0dyaWQ7XHJcbiAgICAgICAgbllEaWZmID0gbllNb3VzZU9uSGVhZGVyIC0gbllCb3JkZXI7XHJcblxyXG4gICAgICAgIGlmKGJCb3JkZXIgJiYgTWF0aC5hYnMobllEaWZmKSA8PSBQaW5uZWRDb2x1bW4uWV9SRVNJWkVfU0VOU0lUSVZJVFkpXHJcbiAgICAgICAge1xyXG4gICAgICAgICAgcmV0dXJuIG5Sb3c7XHJcbiAgICAgICAgfVxyXG5cclxuICAgICAgICBpZighYkJvcmRlciAmJiBuWUJvcmRlciAtIG5IUm93R3JpZCA8PSBuWU1vdXNlT25IZWFkZXIgJiYgbllNb3VzZU9uSGVhZGVyIDw9IG5ZQm9yZGVyKSB7XHJcblxyXG4gICAgICAgICAgaWYoYXJYWU9uQ2VsbCAhPT0gdW5kZWZpbmVkKSB7XHJcbiAgICAgICAgICAgIGFyWFlPbkNlbGxbMF0gPSBlLmNsaWVudFggLSBuWDtcclxuICAgICAgICAgICAgYXJYWU9uQ2VsbFsxXSA9IG5ZTW91c2VPbkhlYWRlciAtIG5ZQm9yZGVyICsgbkhSb3dHcmlkO1xyXG4gICAgICAgICAgfVxyXG5cclxuICAgICAgICAgIHJldHVybiBuUm93O1xyXG4gICAgICAgIH1cclxuICAgICAgfVxyXG4gICAgfVxyXG5cclxuICAgIHJldHVybiAtMTtcclxuICB9XHJcbn1cclxuIl19