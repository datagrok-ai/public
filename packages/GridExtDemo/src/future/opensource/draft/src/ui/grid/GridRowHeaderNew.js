import * as ui from "datagrok-api/ui";
import * as DG from "datagrok-api/dg";
import {GridUtils} from "../../utils/GridUtils";
import {ErrorUtils} from "../../utils/ErrorUtils";
import {DGUtils} from "../../utils/DGUtils";
import {toDart} from "datagrok-api/dg";
import {TextUtils} from "../../utils/TextUtils";
import * as grok from "datagrok-api/grok";



export function closeAllPinnedColumns(grid)
{
 ErrorUtils.verifyClass(grid, DG.Grid);

 const dart = toDart(grid);
    let colPinned = null;
    const nPinnedColCount =dart.m_arRowHeaders === undefined ? 0 : dart.m_arRowHeaders.length;
    for(let n=0; n<nPinnedColCount; ++n) {
        colPinned = dart.m_arRowHeaders[0];
        colPinned.close();
    }
}


export function installPinnedColumns(grid)
{
    ErrorUtils.verifyClass(grid, DG.Grid);
    closeAllPinnedColumns(grid);

    let colGrid = null;
    let settings = null;
    const lstCols = grid.columns;
    const arPinnedCols = new Array();

    for(let nCol=0;nCol<lstCols.length; ++nCol) {
        colGrid = lstCols.byIndex(nCol);
        settings = colGrid.settings;
        if(settings !== null && settings !== undefined && settings.isPinned) {
            arPinnedCols.push(colGrid);
        }
    }

    arPinnedCols.sort((colOne, colTwo) => {
        if(colOne.settings.idxPinned === colTwo.settings.idxPinned)
            throw new Error("Pinned indices cannot be equal for different columns");

        return colOne.settings.idxPinned < colTwo.settings.idxPinned ? -1 : 1;
    });

    for(let n=0;n<arPinnedCols.length; ++n) {
        colGrid = arPinnedCols[n];
        let header = new GridRowHeaderNew(colGrid);
    }
}


function getRenderer(cell) {
    const colGrid = cell.gridColumn;
    if (colGrid === null || colGrid === undefined) {
        throw new Error('Grid cell is detached from the Grid column');
    }
    const dart = toDart(colGrid);

    let renderer = dart.m_renderer;
    if(renderer !== undefined) {
        return renderer;
    }

    renderer = cell.renderer;
    return renderer;
}

function getGrid(colGrid) {
    ErrorUtils.verifyClass(colGrid, DG.GridColumn);

    let grid = colGrid.grid;
    if( grid === null) {
        const dart = toDart(colGrid);
        if(dart.m_grid !== undefined)
         return dart.m_grid;
    }

    return grid;
}

export function isPinnedColumn(colGrid)
{
    ErrorUtils.verifyClass(colGrid, DG.GridColumn);
    const grid = getGrid(colGrid);
    const dart = toDart(grid);

    if(dart.m_arRowHeaders === undefined)
     return false;

    let colPinned = null;
    const nPinnedColCount = dart.m_arRowHeaders.length;
    for(let n=0; n<nPinnedColCount; ++n) {
        colPinned = dart.m_arRowHeaders[n];
        if(toDart(colPinned.m_colGrid) === toDart(colGrid))
            return true;
    }

    return false;
}


export class GridRowHeaderNew {
    constructor(colGrid) {
        ErrorUtils.verifyClass(colGrid, DG.GridColumn);

        if(colGrid.cellType === "html") {

            const dartColG = DG.toDart(colGrid);
            const renderer = dartColG.m_renderer;
            if(renderer === undefined)
             throw new Error("HTML columns cannot be pinned.");
        }

        const grid = getGrid(colGrid);
        if(grid === null){
            throw new Error("Column '" + colGrid.name + "' is not attached to the grid.");
        }

        if(isPinnedColumn(colGrid))
            throw new Error("Column '" + colGrid.name + "' is already pinned.");

        const dart = toDart(grid);

        if(dart.m_arRowHeaders === undefined)
            dart.m_arRowHeaders = [];

        dart.m_arRowHeaders.push(this);

        const viewTable = grid.view;
        const dframe = grid.dataFrame;

        const nW = colGrid.width;

        this.m_colGrid = colGrid;
        try{colGrid.visible = false;}
        catch(e) {}

        if(colGrid.settings === null || colGrid.settings === undefined)
         colGrid.settings = {};

        colGrid.settings.isPinned = true; //this will be saved with the layout
        colGrid.settings.idxPinned = dart.m_arRowHeaders.length -1;

        grid.canvas.style.left = (grid.canvas.offsetLeft + nW).toString() + "px";
        grid.overlay.style.left= (grid.overlay.offsetLeft + nW).toString() + "px";

        const eCanvasThis = ui.canvas(nW, 1500);
        grid.canvas.parentNode.insertBefore(eCanvasThis, grid.canvas);
        this.m_root = eCanvasThis;

        const colGrid0 = grid.columns.byIndex(0);
        if(colGrid0 !== null && colGrid0 !== undefined) {//DG Bug from reading layout
            try{colGrid0.visible = false;}
            catch(e) {
                console.error("ERROR: Couldn't set visible property.");
            }
        }
        //const gRowHeader = eCanvasRowHeader.getContext('2d');
        //OnResize Row header
        const headerThis = this;
        this.m_observerResize = new ResizeObserver(entries => {
            const g = eCanvasThis.getContext('2d');
            for (let entry of entries) {
                headerThis.paint(g, grid);
           }
        });

       this.m_observerResize.observe(eCanvasThis);


        //OnResize Row header
        this.m_observerResizeGrid = new ResizeObserver(entries => {

            const bCurrent = DGUtils.isCurrentView(viewTable);
            if(!bCurrent)
                return;

            const g = eCanvasThis.getContext('2d');
            for (let entry of entries) {
                headerThis.paint(g, grid);
            }
        });

        this.m_observerResizeGrid.observe(grid.root);

        const scrollVert = grid.vertScroll;
        this.m_handlerVScroll = scrollVert.onValuesChanged.subscribe((_) => {
            const g = eCanvasThis.getContext('2d');
            headerThis.paint(g, grid);
         });


        this.m_handlerRowsFiltering = dframe.onRowsFiltering.subscribe(function () {
            setTimeout(() => {
                const g = eCanvasThis.getContext('2d');
                headerThis.paint(g, grid);
            }, 100);

        });

        this.m_handlerCurrRow = dframe.onCurrentRowChanged.subscribe(function () {
                //const nRow = dframe.currentRow.idx;
                const g = eCanvasThis.getContext('2d');
                headerThis.paint(g, grid);
            }
        );


        this.m_handlerSel = dframe.onSelectionChanged.subscribe(function (e) {
              const g = eCanvasThis.getContext('2d');
              headerThis.paint(g, grid);
            }
        );

        this.m_handlerRowResized = grid.onRowsResized.subscribe(function (e) {
             const g = eCanvasThis.getContext('2d');
             headerThis.paint(g, grid);
            }
        );


        let nHResizeRowsBeforeDrag = -1;
        let nResizeRowGridDragging = -1;
        let nYResizeDraggingAnchor = -1;

        let nResizeRowGridMoving = -1;


        let nYDraggingAnchor = -1;
        let nRowGridDragging = -1;
        let arXYMouseOnCellDown = [-2, -2];
        let arXYMouseOnCellUp = [-1, -1];

        this.m_handlerMouseDown = rxjs.fromEvent(document, 'mousedown').subscribe((e) => {

            if(!DGUtils.isCurrentView(viewTable))
                return;

            if(e.button !== 0)
                return;

            nResizeRowGridMoving = -1;

            const bAddToSel = e.ctrlKey || e.shiftKey;

            let nRowGrid = bAddToSel ? -1 : GridRowHeaderNew.hitTestRows(eCanvasThis, grid, e, true);
            if (nRowGrid >= 0) {
                nResizeRowGridDragging = nRowGrid;
                nYResizeDraggingAnchor = e.clientY;
                nHResizeRowsBeforeDrag = GridUtils.getRowHeight(grid);
            }
            else
            {
                nRowGrid = GridRowHeaderNew.hitTestRows(eCanvasThis, grid, e, false, arXYMouseOnCellDown);

                nRowGridDragging = nRowGrid;
                nYDraggingAnchor = e.clientY;

                const cell = grid.cell(colGrid.name, nRowGrid);
                const renderer = getRenderer(cell);//toDart(colGrid).m_renderer;
                if(renderer !== undefined && renderer !== null && renderer.onMouseDown !== undefined) {
                    renderer.onMouseDown(cell, e, arXYMouseOnCellDown[0], arXYMouseOnCellDown[1]);
                }

            }

            e.preventDefault();
            e.stopPropagation();
        });


        this.m_handlerMouseUp = rxjs.fromEvent(document, 'mouseup').subscribe((e) => {

            if(!DGUtils.isCurrentView(viewTable))
                return;

            if(e.button !== 0)
                return;

            nResizeRowGridDragging = -1;
            nResizeRowGridDragging = -1;
            nYResizeDraggingAnchor = -1;

            nResizeRowGridMoving = -1;

            document.body.style.cursor = "auto";

            if(nRowGridDragging >= 0) {
                const bAddToSel = e.ctrlKey || e.shiftKey;

                const nRowGrid = GridRowHeaderNew.hitTestRows(eCanvasThis, grid, e, false, arXYMouseOnCellUp);
                if(!bAddToSel && nRowGrid === nRowGridDragging) {

                    let cellRH = null;
                    try{cellRH = grid.cell("", nRowGrid);}
                    catch(e) {
                     let colG = null;
                     const lstCols = grid.columns;
                     for(let nC=1; nC<lstCols.length; ++nC) {
                        colG = lstCols.byIndex(nC);
                        cellRH = grid.cell(colG.name, nRowGrid);
                        if(cellRH !== null)
                          break;
                     }
                    }

                    const nRowTable = cellRH.tableRowIndex;

                    dframe.currentRow = nRowTable;
                }
                else
                {
                    const bitsetSel = dframe.selection;

                    if(!bAddToSel)
                        bitsetSel.setAll(false, true);

                    const nRowMin = nRowGridDragging < nRowGrid ? nRowGridDragging : nRowGrid;
                    const nRowMax = nRowGridDragging > nRowGrid ? nRowGridDragging : nRowGrid;
                    let cellRH = null;
                    let nRowTable = -1;
                    for(let nRow=nRowMin; nRow<=nRowMax; ++nRow) {

                        try{cellRH = grid.cell("", nRow);}
                        catch(e) {
                            let colG = null;
                            const lstCols = grid.columns;
                            for(let nC=1; nC<lstCols.length; ++nC) {
                                colG = lstCols.byIndex(nC);
                                cellRH = grid.cell(colG.name, nRowGrid);
                                if(cellRH !== null)
                                    break;
                            }
                        }

                        nRowTable = cellRH.tableRowIndex;
                        bitsetSel.set(nRowTable, true, true);
                    }
                }


                const cell = grid.cell(colGrid.name, nRowGrid);
                const renderer = getRenderer(cell);
                if(renderer !== undefined && renderer !== null && renderer.onMouseUp !== undefined) {
                    renderer.onMouseUp(cell, e, arXYMouseOnCellUp[0], arXYMouseOnCellUp[1]);
                }

                if(arXYMouseOnCellUp[0] === arXYMouseOnCellDown[0] && arXYMouseOnCellDown[1] === arXYMouseOnCellUp[1]) {
                    if(renderer !== undefined && renderer !== null && renderer.onClick !== undefined) {
                        renderer.onClick(cell, e, arXYMouseOnCellUp[0], arXYMouseOnCellUp[1]);
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

            if(cellCurrent !== null) {
                const renderer = getRenderer(cellCurrent);
                if (renderer !== undefined && renderer !== null) {
                    renderer.onMouseLeave(cellCurrent, e, -1, -1);
                }
                cellCurrent = null;
            }
        });


        this.m_handlerMouseMove = rxjs.fromEvent(document, 'mousemove').subscribe((e) => {

            if(!DGUtils.isCurrentView(viewTable))
                return;

            const bDragging = nResizeRowGridDragging >= 0;
            if (bDragging) {

                //console.log("Dragging : " + headerThis.m_strColName);

                const nYDiff = e.clientY - nYResizeDraggingAnchor;
                let nHRowGrid = nHResizeRowsBeforeDrag + nYDiff;

                if (nHRowGrid < GridRowHeaderNew.MIN_ROW_HEIGHT)
                    nHRowGrid = GridRowHeaderNew.MIN_ROW_HEIGHT;
                else if (nHRowGrid > GridRowHeaderNew.MAX_ROW_HEIGHT)
                    nHRowGrid = GridRowHeaderNew.MAX_ROW_HEIGHT;

                let g = eCanvasThis.getContext('2d');
                g.fillStyle = "white";
                const nHHeaderCols = GridUtils.getColumnHeaderHeight(grid);
                g.fillRect(0,nHHeaderCols, eCanvasThis.offsetWidth, eCanvasThis.offsetHeight);

                GridUtils.setRowHeight(grid, nHRowGrid);

                let header = null;
                const ar = grid.dart.m_arRowHeaders;
                for(let n=0; n<ar.length; ++n) {
                    header = ar[n];
                    g = header.m_root.getContext('2d');
                    header.paint(g, grid);
                }

                try{grid.columns.byIndex(0).visible = false;}
                catch(e) {}//temporary addressed the DG bug
                return;
            }


            const arXYOnCell = [-1,-1];
            let nRowGrid = GridRowHeaderNew.hitTestRows(eCanvasThis, grid, e, false, arXYOnCell);
            if(nRowGrid >= 0) {
                const cell = grid.cell(colGrid.name, nRowGrid);
                const renderer = getRenderer(cell);

                if (renderer !== undefined && renderer !== null) {

                    if (cellCurrent === null && renderer.onMouseEnter !== undefined) {
                        renderer.onMouseEnter(cell, e, arXYOnCell[0], arXYOnCell[1]);
                    }

                    if (cellCurrent !== null && nRowGrid !== cellCurrent.gridRow) {
                        if(renderer.onMouseLeave !== undefined)
                         renderer.onMouseLeave(cellCurrent, e, -1, -1);

                        if(renderer.onMouseEnter !== undefined)
                         renderer.onMouseEnter(cell, e, arXYOnCell[0], arXYOnCell[1]);
                    }

                    if (renderer.onMouseMove !== undefined) {
                        renderer.onMouseMove(cell, e, arXYOnCell[0], arXYOnCell[1]);
                    }
                }

                cellCurrent = cell;
            }
            else if (cellCurrent !== null) {
                const renderer = getRenderer(cellCurrent);
                if (renderer !== undefined && renderer !== null && renderer.onMouseLeave !== undefined) {
                    renderer.onMouseLeave(cellCurrent, e, -1, -1);
                }

              cellCurrent = null;
            }

            nRowGrid = GridRowHeaderNew.hitTestRows(eCanvasThis, grid, e, true);
            if (nRowGrid >= 0) {
                nResizeRowGridMoving = nRowGrid;
                document.body.style.cursor = "row-resize";
                return;
            }

            if(nResizeRowGridMoving >= 0) {
                nResizeRowGridMoving = -1;
                document.body.style.cursor = "auto";
            }
        });

        const thisRowHeader = this;
        this.m_handlerContextMenu = grok.events.onContextMenu.subscribe((args) => {
        //this.m_handlerContextMenu = rxjs.fromEvent(viewTable.root, 'contextmenu').subscribe((e) => {
            const e = args.causedBy;
            const elem = document.elementFromPoint(e.clientX, e.clientY);//e.offsetY);
            const b = elem === eCanvasThis;

            if(b) {
                let menu = args.args.menu;//DG.Menu.popup();

                menu = menu.item("Unpin Column", (str) => {
                    thisRowHeader.close();
                });

                menu = menu.item("Unpin All Columns", (str) => {
                    closeAllPinnedColumns(grid);
                });

                e.preventDefault();
                e.stopPropagation();
            }
        });
    }

    isPinned()
    {
        return this.m_colGrid !== null;
    }

    getGridColumn()
    {
        return this.m_colGrid;
    }

    close() {

        if(this.m_colGrid === null)
            throw new Error("Column has already been unpinned");

        this.m_observerResizeGrid.disconnect();
        this.m_observerResizeGrid = null;

        this.m_observerResize.disconnect();
        this.m_observerResize = null;

        this.m_handlerVScroll.unsubscribe();
        this.m_handlerVScroll = null;

        this.m_handlerRowResized.unsubscribe();
        this.m_handlerRowResized = null;

        this.m_handlerRowsFiltering.unsubscribe();
        this.m_handlerRowsFilterin = null;

        this.m_handlerCurrRow.unsubscribe();
        this.m_handlerCurrRow = null;

        this.m_handlerSel.unsubscribe();
        this.m_handlerSel = null;

        this.m_handlerMouseDown.unsubscribe();
        this.m_handlerMouseDown = null;

        this.m_handlerMouseUp.unsubscribe();
        this.m_handlerMouseUp = null;

        this.m_handlerMouseMove.unsubscribe();
        this.m_handlerMouseMove = null;

        this.m_handlerMouseLeave.unsubscribe();
        this.m_handlerMouseLeave = null;

        this.m_handlerContextMenu.unsubscribe();
        this.m_handlerContextMenu = null;

        const grid = getGrid(this.m_colGrid);
        const dart = toDart(grid);
        const ar = dart.m_arRowHeaders;
        const nIdx = ar.indexOf(this);
        ar.splice(nIdx, 1);

        let nIdxPinned = -1;
        let colGridTmp= null;
        for(let n=nIdx; n<ar.length; ++n) {
            colGridTmp = ar[n];
            nIdxPinned =  colGridTmp.m_colGrid.settings.idxPinned;
            colGridTmp.m_colGrid.settings.idxPinned = n;
        }

        this.m_colGrid.settings.idxPinned = -1;
        this.m_colGrid.settings.isPinned = false;
        try{this.m_colGrid.visible = true;}
        catch(e) {}
        this.m_colGrid = null;


        if (dart.m_arRowHeaders.length === 0) {
            const colGrid0 = grid.columns.byIndex(0);
            if(colGrid0 !== null && colGrid0 !== undefined) {
                try{colGrid0.visible = true;}//DG Bug from reading layout
                catch(e) {
                    console.error("ERROR: Couldn't set visible property.");
                }
            }

         //   grid.columns.byIndex(0).visible = true;

        }
        grid.canvas.style.left = (grid.canvas.offsetLeft - this.m_root.offsetWidth).toString() + "px";
        grid.overlay.style.left= (grid.overlay.offsetLeft - this.m_root.offsetWidth).toString() + "px";

        this.m_root.parentNode.removeChild(this.m_root);
        this.m_root = null;
    }

    getWidth() {
        return this.m_root.offsetWidth;
    }


    paint(g, grid)
    {
        //const nWDiv = entry.contentBoxSize ? entry.contentBoxSize[0].inlineSize : entry.contentRect.width;

        const dframe = grid.dataFrame;
        const nW = this.m_root.offsetWidth;

        g.fillStyle = "white";
        //const nHHeaderCols = GridUtils.getColumnHeaderHeight(grid);
        g.fillRect(0,0, this.m_root.offsetWidth, this.m_root.offsetHeight);

        if(this.m_colGrid.name === null)
            return;

        const bitsetFilter = dframe.filter;
        if(bitsetFilter.falseCount === dframe.rowCount)
            return;

/* column header
        const rendererRowColHeader = this.getRowsHeaderRenderer();

        if (rendererRowColHeader !== null) {

            if (rendererRowColHeader instanceof ButtonGridColumnHeaderRenderer)
                rendererRowColHeader.setFilterEnabled(true);

            const colGrid = grid.columns.byName(this.m_colGrid.name);

            const nHHeaderCols = GridUtils.getColumnHeaderHeight(grid);
            rendererRowColHeader.paint(g, colGrid, 0, 0, nW, nHHeaderCols);
        }
*/
        let str = TextUtils.trimText(this.m_colGrid.name, g, nW);
        g.font = "bold 13px Roboto, Roboto Local";
        const tm = g.measureText(str);
        const nWLabel = tm.width;

        const nAscent = Math.abs(tm.actualBoundingBoxAscent);
        const nDescent = tm.actualBoundingBoxDescent;
        const nHFont =  nAscent + nDescent;// + 2*nYInset;

        let nX = 0;
        let nY = 0;
        const nH = GridUtils.getColumnHeaderHeight(grid);
        g.fillStyle = "Black";
        g.fillText(str, nX + ((nW - nWLabel)>>1), nY + nH-2);


         const nRowCurrent =  dframe.currentRow.idx;
         const bitsetSel = dframe.selection;

         const nGridRowCount = dframe.filter.trueCount;
         const nGridfalseCount = dframe.filter.falseCount;

         const scrollV = grid.vertScroll;
         const nRowMin = Math.floor(scrollV.min);
         let nRowMax = Math.ceil(scrollV.max);

         let nHH = grid.root.offsetHeight - GridUtils.getColumnHeaderHeight(grid);//.style.height;
         let nHRow = GridUtils.getRowHeight(grid);
         let nRCount = Math.round(nHH/nHRow) +1;
         nRowMax = nRowMin + nRCount;

         if(nRowMax >= nGridRowCount)
          nRowMax = nGridRowCount -1;

            //console.log(nRowMin + " " + nRowMax);

         const nYOffset = GridUtils.getColumnHeaderHeight(grid);
         const nHRowGrid = GridUtils.getRowHeight(grid);
         let cellRH = null;


         let nRowTable = -1;
          nY = -1;
         let bSel = false;
         let bFiltered = false;
         let renderer = null;
         for(let nRG=nRowMin; nRG<=nRowMax; ++nRG)
            {
                try{cellRH = grid.cell(this.m_colGrid.name, nRG);}
                catch(e)     //to address DG bug when everything is filtered
                {
                    continue;
                }

                if(cellRH.tableRowIndex === undefined)//DG bug
                    continue;

                nRowTable = cellRH.tableRowIndex;
                nY = nYOffset + (nRG - nRowMin)*nHRowGrid;

                renderer = DG.toDart(cellRH.gridColumn).m_renderer;
                if(renderer === undefined)
                 renderer = cellRH.renderer;

                try{renderer.render(g, 0, nY, nW, nHRowGrid, cellRH, cellRH.style);}
                catch(e) {
                  throw e;
                }


                bSel = bitsetSel.get(nRowTable);
                if(bSel)
                {
                    g.globalAlpha = 0.2;
                    g.fillStyle = GridRowHeaderNew.SELECTION_COLOR;
                    g.fillRect(0, nY, nW, nHRowGrid);
                    g.globalAlpha = 1;
                }

                if(nRowCurrent === nRowTable)
                {
                    g.globalAlpha = 0.2;
                    g.fillStyle = GridRowHeaderNew.ACTIVE_CELL_COLOR;
                    g.fillRect(0, nY, nW, nHRowGrid);
                    g.globalAlpha = 1;
                }
            }
    }

}

GridRowHeaderNew.MIN_ROW_HEIGHT = 20;
GridRowHeaderNew.MAX_ROW_HEIGHT = 500;
GridRowHeaderNew.SELECTION_COLOR = DG.Color.toRgb(DG.Color.colSelection);// "rgba(237, 220, 88, 0.15)";
GridRowHeaderNew.ACTIVE_CELL_COLOR = DG.Color.toRgb(DG.Color.currentRow);// "rgba(153, 237, 82, 0.35)";
GridRowHeaderNew.Y_RESIZE_SENSITIVITY = 2;

GridRowHeaderNew.convertMouseY2Cell = function(e, nRowGrid) {

}

GridRowHeaderNew.hitTestRows = function (eCanvasPinned, grid, e, bBorder, arXYOnCell)
{
    const rect = eCanvasPinned.getBoundingClientRect();
    const scrollLeft= window.pageXOffset || document.documentElement.scrollLeft;
    const scrollTop = window.pageYOffset || document.documentElement.scrollTop;
    const nY = rect.top  + scrollTop;
    const nX = rect.left + scrollLeft;

    if(nX <= e.clientX && e.clientX <= nX + eCanvasPinned.offsetWidth)   //on the rows header
    {
        const nHHeaderCols = GridUtils.getColumnHeaderHeight(grid);
        const nHRowGrid = GridUtils.getRowHeight(grid);

        const scroll = grid.vertScroll;
        const nRowMin = Math.floor(scroll.min);
        const nRecorfCount = grid.dataFrame.rowCount;
        let nRowMax = Math.min(Math.floor(scroll.max) +1, nRecorfCount-1);
        const nYMouseOnHeader = e.clientY - nY;

        let nYBorder = -1;
        let nYDiff = -1;

        for(let nRow=nRowMin+1; nRow<= nRowMax; ++nRow)
        {
            nYBorder = nHHeaderCols + (nRow - nRowMin)*nHRowGrid;
            nYDiff = nYMouseOnHeader - nYBorder;

            if(bBorder && Math.abs(nYDiff) <= GridRowHeaderNew.Y_RESIZE_SENSITIVITY)
            {
                return nRow -1;
            }

            if(!bBorder && nYBorder - nHRowGrid <= nYMouseOnHeader && nYMouseOnHeader <= nYBorder) {
                if(arXYOnCell !== undefined) {
                    arXYOnCell[0] = e.clientX - nX;
                    arXYOnCell[1] = nYMouseOnHeader - nYBorder + nHRowGrid;
                }
                return nRow - 1;
            }
        }
    }

    return -1;
}
