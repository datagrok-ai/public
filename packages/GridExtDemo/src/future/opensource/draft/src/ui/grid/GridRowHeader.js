import * as ui from "datagrok-api/ui";
import * as DG from "datagrok-api/dg";
import {ButtonGridColumnHeaderRenderer} from "../../entity/ui/ButtonGridColumnHeaderRenderer";
import {GridUtils} from "../../utils/GridUtils";
import {SemType} from "../../entity/SemType";
import {ErrorUtils} from "../../utils/ErrorUtils";
import {DGUtils} from "../../utils/DGUtils";

function removePanelHeader(eDivRowHeader)
{
    const eDivRowHeaderP = eDivRowHeader.parentElement.parentElement;
    const arChildren = eDivRowHeaderP.childNodes;
    const nChildCount = arChildren.length;
    for (let nChild = 0; nChild < nChildCount; ++nChild) {
        if (arChildren[nChild] !== eDivRowHeader) {
            eDivRowHeaderP.removeChild(arChildren[nChild]);
            return true;
        }
    }
    return false;
}



export class GridRowHeaderEmpty extends DG.JsViewer {
    constructor() {
        super();

        const eCanvasRowHeader = ui.canvas(350, 500);
        this.root.appendChild(eCanvasRowHeader);
    }

}

function registerEmptyViewer()
{
    grok.shell.registerViewer("emptyviewer", "Empty Viewer", function() {
       return new GridRowHeaderEmpty();
    });
}

registerEmptyViewer();

export class GridRowHeader {//extends DG.JsViewer {
    constructor(viewTable, strColName) {
    //super();
        ErrorUtils.verifyClass(viewTable, DG.TableView);
        ErrorUtils.verifyType(strColName, String);

        this.m_strColName = strColName;

        const manager = viewTable.dockManager;
        const grid = viewTable.grid;
        const dframe = grid.dataFrame;


        const col = dframe.columns.byName(this.m_strColName);
        const typeSem = GridUtils.getSemType(col);
        const config = GridUtils.getRenderersConfig(grid);

        this.m_rendererCell = typeSem === null || config === null ? null : config.getCellRenderer(typeSem.constructor);


        this.m_bitsetSel = dframe.selection;

        this.m_viewTable = viewTable;

        const eDivRowHeader = document.createElement("div");
        this.m_root = eDivRowHeader;
        //this.root = eDivRowHeader;

        const eCanvasRowHeader = ui.canvas(2350, 1500);
        eDivRowHeader.appendChild(eCanvasRowHeader);


        const nodeRoot = manager.rootNode;
        const nodeGrid = manager.findNode(grid.root);

        const node = manager.dock(eDivRowHeader, DG.DOCK_TYPE.LEFT, nodeGrid, "rows header", 0.11);
        const gRowHeader = eCanvasRowHeader.getContext('2d');
        this.m_graphics = gRowHeader;

        GridUtils.setRowsHeader(grid, this, false);

        removePanelHeader(eDivRowHeader);

        //OnResize Row header
        const headerThis = this;
        this.m_observerResize = new ResizeObserver(entries => {
            for (let entry of entries) {
                headerThis.paint(gRowHeader, grid);
           }
        });

       this.m_observerResize.observe(eDivRowHeader);
       grid.columns.byIndex(0).visible = false;

        //OnResize Row header

        this.m_eGridParent = grid.root;
        this.m_observerResizeGrid = new ResizeObserver(entries => {

            const bCurrent = DGUtils.isCurrentView(viewTable);
            if(!bCurrent)
                return;

            for (let entry of entries) {

                if( grid.root.offsetParent !== null)
                {
                  //let aaa = grid.root.offsetParent;

                  if(eDivRowHeader.offsetParent === null)
                  {
                   manager.dock(eDivRowHeader, DG.DOCK_TYPE.LEFT, nodeGrid, "rows header", 0.11);
                   removePanelHeader(eDivRowHeader);
                   headerThis.m_observerResize.observe(eDivRowHeader);
                  }
                }
                else
                {
                   // eDivRowHeader.style.display = "none";

                   headerThis.m_observerResize.disconnect();
                   manager.close(eDivRowHeader);
                }


                //headerThis.paint(gRowHeader, grid);
            }
        });

        this.m_observerResizeGrid.observe(grid.root);


        const scrollVert = grid.vertScroll;
        this.m_handlerVScroll = scrollVert.onValuesChanged.subscribe((_) => {
            headerThis.paint(gRowHeader, grid);
         });


        this.m_handlerRowsFiltering = dframe.onRowsFiltering.subscribe(function () {
            setTimeout(() => {
                headerThis.paint(gRowHeader, grid);
            }, 100);

        });

        this.m_handlerCurrRow = dframe.onCurrentRowChanged.subscribe(function () {
                const nRow = dframe.currentRow.idx;
                headerThis.paint(gRowHeader, grid);
            }
        );


        this.m_handlerSel = dframe.onSelectionChanged.subscribe(function (e) {
            //headerThis.m_bitsetSel = dframe.selection;
            headerThis.paint(gRowHeader, grid);
            }
        );

        this.m_handlerRowResized = grid.onRowsResized.subscribe(function (e) {
                //headerThis.m_bitsetSel = dframe.selection;
                headerThis.paint(gRowHeader, grid);
            }
        );


        let nHResizeRowsBeforeDrag = -1;
        let nResizeRowGridDragging = -1;
        let nYResizeDraggingAnchor = -1;

        let nResizeRowGridMoving = -1;


        let nYDraggingAnchor = -1;
        let nRowGridDragging = -1;


        this.m_handlerMouseDown = rxjs.fromEvent(document, 'mousedown').subscribe((e) => {

            if(!DGUtils.isCurrentView(viewTable))
                return;

            if(e.button !== 0)
                return;

            nResizeRowGridMoving = -1;

            const bAddToSel = e.ctrlKey || e.shiftKey;

            let nRowGrid = bAddToSel ? -1 : GridRowHeader.hitTestRows(eDivRowHeader, grid, e, true);
            if (nRowGrid >= 0) {
                nResizeRowGridDragging = nRowGrid;
                nYResizeDraggingAnchor = e.clientY;
                nHResizeRowsBeforeDrag = GridUtils.getRowHeight(grid);
            }
            else
            {
                nRowGrid = GridRowHeader.hitTestRows(eDivRowHeader, grid, e, false);

                nRowGridDragging = nRowGrid;
                nYDraggingAnchor = e.clientY;
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

                const nRowGrid = GridRowHeader.hitTestRows(eDivRowHeader, grid, e, false);
                if(!bAddToSel && nRowGrid === nRowGridDragging) {

                    const cellRH = grid.cell("", nRowGrid);
                    const nRowTable = cellRH.tableRowIndex;

                    dframe.currentRow = nRowTable;
                }
                else
                {
                    if(!bAddToSel)
                        headerThis.m_bitsetSel.setAll(false, true);

                    const nRowMin = nRowGridDragging < nRowGrid ? nRowGridDragging : nRowGrid;
                    const nRowMax = nRowGridDragging > nRowGrid ? nRowGridDragging : nRowGrid;
                    let cellRH = null;
                    let nRowTable = -1;
                    for(let nRow=nRowMin; nRow<=nRowMax; ++nRow) {
                        cellRH = grid.cell("", nRow);
                        nRowTable = cellRH.tableRowIndex;

                        headerThis.m_bitsetSel.set(nRowTable, true, true);
                    }
                }


                nRowGridDragging = -1;
                nYDraggingAnchor = -1;
            }
            //e.preventDefault();
            //e.stopPropagation();

        });


        this.m_handlerMouseMove = rxjs.fromEvent(document, 'mousemove').subscribe((e) => {

            if(!DGUtils.isCurrentView(viewTable))
                return;

            const bDragging = nResizeRowGridDragging >= 0;
            if (bDragging) {

                //console.log("Dragging : " + headerThis.m_strColName);

                const nYDiff = e.clientY - nYResizeDraggingAnchor;
                let nHRowGrid = nHResizeRowsBeforeDrag + nYDiff;

                if (nHRowGrid < GridRowHeader.MIN_ROW_HEIGHT)
                    nHRowGrid = GridRowHeader.MIN_ROW_HEIGHT;
                else if (nHRowGrid > GridRowHeader.MAX_ROW_HEIGHT)
                    nHRowGrid = GridRowHeader.MAX_ROW_HEIGHT;

                gRowHeader.fillStyle = "white";
                const nHHeaderCols = GridUtils.getColumnHeaderHeight(grid);
                gRowHeader.fillRect(0,nHHeaderCols, eDivRowHeader.offsetWidth, eDivRowHeader.offsetHeight);

                GridUtils.setRowHeight(grid, nHRowGrid);
                //headerThis.paint(gRowHeader, grid);
                return;
            }


            if (nResizeRowGridDragging >= 0) {
                document.body.style.cursor = "row-resize";
                return;
            }

            const nRowGrid = GridRowHeader.hitTestRows(eDivRowHeader, grid, e, true);
            if (nRowGrid >= 0) {
                nResizeRowGridMoving = nRowGrid;
                document.body.style.cursor = "row-resize";
                return;
            }

            if(nResizeRowGridMoving >= 0) {
                nResizeRowGridMoving = -1;
                document.body.style.cursor = "auto";
            }


            //TOOLTIP
            const nRowGridTT = GridRowHeader.hitTestRows(eDivRowHeader, grid, e, false);

            if(nRowGridTT >= 0)
            {
                const bTooltipOn = GridUtils.isTooltipEnabled(grid);
                if(e.ctrlKey || bTooltipOn)
                {
                    const rect = eDivRowHeader.getBoundingClientRect();
                    const scrollLeft= window.pageXOffset || document.documentElement.scrollLeft;
                    const scrollTop = window.pageYOffset || document.documentElement.scrollTop;
                    const nY = rect.top  + scrollTop;
                    const nX = rect.left + scrollLeft;

                    const cell = grid.cell(headerThis.m_strColName, nRowGridTT);
                    const nXOnCell = headerThis.getWidth() + nX;//cell.bounds.x + cell.gridColumn.width + nX;;
                    const nYOnCell = cell.bounds.y + nY;
                    const col = dframe.columns.byName(headerThis.m_strColName);

                    const typeSem = GridUtils.getSemType(col);
                    if(typeSem instanceof SemType)
                    {
                        const renderer = config.getCellRenderer(typeSem.constructor);
                        const eTooltip = renderer.createTootipContent(cell);
                        ui.tooltip.show(ui.divV([
                            eTooltip
                        ]), nXOnCell, nYOnCell);
                    }

                }
            }


           // const nRowGridTmp = GridRowHeader.hitTestRows(eDivRowHeader, grid, e, false);
            //const renderer = this.getRowsRenderer();
            //   console.log("Moved : " + nRowGrid);

        });

        this.m_handlerContextMenu = rxjs.fromEvent(eCanvasRowHeader, 'contextmenu').subscribe((e) => {

            //const nHColHrader = GridUtils.getColumnHeaderHeight(grid);
            let cell = grid.hitTest(eDivRowHeader.offsetWidth + 1, e.offsetY);
            if (cell === undefined || cell.cellType === null)//bufg in DG , top left cell
                return;

            cell = grid.cell(headerThis.m_strColName, cell.gridRow);

            //const config = GridUtils.getRenderersConfig(grid);
            //const nRowGrid = cell.gridRow;
            const rendererRowHeader = headerThis.getRowsRenderer();

            const nGroupCount = rendererRowHeader.getPopupGroupCount(cell);
            let menu = DG.Menu.popup();
            let bHasContent = false;
            if (nGroupCount > 0) {
                bHasContent = true;
                let strGroupName = "";
                let arItems = [];
                for (let nGroup = 0; nGroup < nGroupCount; ++nGroup) {
                    menu = menu.separator();
                    strGroupName = rendererRowHeader.getPopupGroupName(cell, nGroup);
                    arItems = rendererRowHeader.toPopupGrouptemsArray(cell, nGroup);
                    let nGroupTmp = nGroup;
                    let arItemsTmp = arItems;
                    menu.group(strGroupName).items(arItems, function (strItem) {

                        const nItemClicked = arItemsTmp.indexOf(strItem);
                        rendererRowHeader.onPopupGroupAction(cell, nGroupTmp, nItemClicked);
                    });
               }
              menu = menu.separator();
            }


            const arItems = rendererRowHeader.toPopupItemsArray(cell);
            if (arItems.length > 0) {
                bHasContent = true;

                function onItem(strItem) {
                    const nItemClicked = arItems.indexOf(strItem);
                    rendererRowHeader.onPopupAction(cell, nItemClicked);
                }

                menu.items(arItems, onItem);
            }

            if(bHasContent)
             menu.show();

            e.preventDefault();
        });
    }


    getTableView()
    {
        return this.m_viewTable;
    }

    getColumnName()
    {
        return this.m_strColName;
    }

    setColumnName(strColName)
    {
        ErrorUtils.verifyType(strColName, String);

        this.m_strColName = strColName;
        const g = this.getGraphics();
        this.paint(g, this.m_viewTable.grid);
    }


    close()
    {

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

        this.m_handlerContextMenu.unsubscribe();
        this.m_handlerContextMenu = null;

        const manager = this.m_viewTable.dockManager;
        manager.close(this.m_root);
        const grid = this.m_viewTable.grid;
        GridUtils.setRowsHeader(grid, this, true);

        if(GridUtils.getRowsHeaderCount(grid) === 0)
         grid.columns.byIndex(0).visible = true;
    }

    getWidth() {
        return this.m_root.offsetWidth;
    }


    getGraphics() {
        return this.m_graphics;
    }

    getRowsRenderer() {

        if(this.m_strColName === null)
            return null;


        return this.m_rendererCell;
    }

    getRowsHeaderRenderer() {

        if(this.m_strColName === null)
            return null;

        const grid = this.m_viewTable.grid;
        const dframe = grid.dataFrame;
        const col = dframe.columns.byName(this.m_strColName);
        const typeSem = GridUtils.getSemType(col);
        const config = GridUtils.getRenderersConfig(grid);

        let renderer = config === null ? null : config.getDefaultColHeaderRenderer();//this.getRowsHeaderRenderer();
        if (typeSem instanceof SemType)
          renderer = config.getColHeaderRenderer(typeSem.constructor);

        return renderer;
    }

    /*
    paintRow(g, nRowGrid, grid)
    {

        if(this.m_strColName === null)
            return;

        const dframe = grid.dataFrame;
        const bitsetSel = this.m_bitsetSel;
        const nRowCurrent =  dframe.currentRow.idx;

        const nW = this.m_root.offsetWidth;
        const nHRowGrid = GridUtils.getRowHeight(grid);

        const cellRH = grid.cell(this.m_strColName, nRowGrid);
        const lstCells = cellRH.gridColumn.getVisibleCells();
        for(let c of lstCells)
        {
            let fgfh = 0;
        }

        const nRowTable = cellRH.tableRowIndex;
        const nYOffset = GridUtils.getColumnHeaderHeight(grid);
        const scrollV = grid.vertScroll;
        const nRowMin = scrollV.min;
        const nY = nYOffset + (nRowGrid - nRowMin)*nHRowGrid;


        g.fillStyle = "white";
        g.fillRect(0, nY, nW, nHRowGrid);

        const renderer = this.getRowsRenderer();
        renderer.paint(g, cellRH, 0, nY, nW, nHRowGrid);

        const bSel = bitsetSel.get(nRowTable);
        if(bSel)
        {
            g.fillStyle = GridRowHeader.SELECTION_COLOR;
            g.fillRect(0, nY, nW, nHRowGrid);
        }

        if(nRowCurrent === nRowTable)
        {
            g.fillStyle = GridRowHeader.ACTIVE_CELL_COLOR;
            g.fillRect(0, nY, nW, nHRowGrid);
        }
   } */

    invalidate()
    {
        this.paint(this.m_graphics, this.m_viewTable.grid);
    }

    paint(g, grid)
    {
        //const nWDiv = entry.contentBoxSize ? entry.contentBoxSize[0].inlineSize : entry.contentRect.width;

        const dframe = grid.dataFrame;
        const nW = this.m_root.offsetWidth;

        g.fillStyle = "white";
        //const nHHeaderCols = GridUtils.getColumnHeaderHeight(grid);
        g.fillRect(0,0, this.m_root.offsetWidth, this.m_root.offsetHeight);

        if(this.m_strColName === null)
            return;

        const bitsetFilter = dframe.filter;
        if(bitsetFilter.falseCount === dframe.rowCount)
            return;


        const rendererRowColHeader = this.getRowsHeaderRenderer();

        if (rendererRowColHeader !== null) {

            if (rendererRowColHeader instanceof ButtonGridColumnHeaderRenderer)
                rendererRowColHeader.setFilterEnabled(true);

            const colGrid = grid.columns.byName(this.m_strColName);

            const nHHeaderCols = GridUtils.getColumnHeaderHeight(grid);
            rendererRowColHeader.paint(g, colGrid, 0, 0, nW, nHHeaderCols);
        }


        const rendererRowHeader = this.getRowsRenderer();

      //  if (rendererRowHeader !== null) {

            const nRowCurrent =  dframe.currentRow.idx;
            const bitsetSel = this.m_bitsetSel;

            const nGridRowCount = dframe.filter.trueCount;
            const nGridfalseCount = dframe.filter.falseCount;

            const scrollV = grid.vertScroll;
            const nRowMin = Math.floor(scrollV.min);
            let nRowMax = Math.ceil(scrollV.max);
            if(nRowMax >= nGridRowCount)
                nRowMax = nGridRowCount -1;

            //console.log(nRowMin + " " + nRowMax);

            const nYOffset = GridUtils.getColumnHeaderHeight(grid);
            const nHRowGrid = GridUtils.getRowHeight(grid);
            let cellRH = null;
            let nRowTable = -1;
            let nY = -1;
            let bSel = false;
            let bFiltered = false;
            for(let nRG=nRowMin; nRG<=nRowMax; ++nRG)
            {
                try{cellRH = grid.cell(this.m_strColName, nRG);}
                catch(e)     //to address DG bug when everything is filtered
                {
                    continue;
                }

                if(cellRH.tableRowIndex === undefined)//DG bug
                    continue;

                nRowTable = cellRH.tableRowIndex;
                nY = nYOffset + (nRG - nRowMin)*nHRowGrid;

                if(rendererRowHeader === null) {
                    try{cellRH.renderer.render(g, 0, nY, nW, nHRowGrid, cellRH, cellRH.style);}
                    catch(e) {
                        throw e;
                    }
                }
                else {
                    try {
                        rendererRowHeader.paint(g, cellRH, 0, nY, nW, nHRowGrid);
                    } catch (e) {
                        let fghfg = 0;
                        //throw e;
                    }
                }

                bSel = bitsetSel.get(nRowTable);
                if(bSel)
                {
                    g.fillStyle = GridRowHeader.SELECTION_COLOR;
                    g.fillRect(0, nY, nW, nHRowGrid);
                }

                if(nRowCurrent === nRowTable)
                {
                    g.fillStyle = GridRowHeader.ACTIVE_CELL_COLOR;
                    g.fillRect(0, nY, nW, nHRowGrid);
                }
            }
        //}
    }

}

GridRowHeader.MIN_ROW_HEIGHT = 20;
GridRowHeader.MAX_ROW_HEIGHT = 500;
GridRowHeader.SELECTION_COLOR = "rgba(237, 220, 88, 0.15)";
GridRowHeader.ACTIVE_CELL_COLOR = "rgba(153, 237, 82, 0.35)";
GridRowHeader.Y_RESIZE_SENSITIVITY = 2;

GridRowHeader.hitTestRows = function (eDivRowHeader, grid, e, bBorder)
{
    const rect = eDivRowHeader.getBoundingClientRect();
    const scrollLeft= window.pageXOffset || document.documentElement.scrollLeft;
    const scrollTop = window.pageYOffset || document.documentElement.scrollTop;
    const nY = rect.top  + scrollTop;
    const nX = rect.left + scrollLeft;

    if(nX <= e.clientX && e.clientX <= nX + eDivRowHeader.offsetWidth)   //on the rows header
    {
        const nHHeaderCols = GridUtils.getColumnHeaderHeight(grid);
        const nHRowGrid = GridUtils.getRowHeight(grid);

        const scroll = grid.vertScroll;
        const nRowMin = scroll.min;
        let nRowMax = Math.min(Math.floor(scroll.max) +1,  grid.dataFrame.rowCount -1);
        const nYMouseOnHeader = e.clientY - nY;

        let nYBorder = -1;
        let nYDiff = -1;

        for(let nRow=nRowMin+1; nRow<= nRowMax; ++nRow)
        {
            nYBorder = nHHeaderCols + (nRow - nRowMin)*nHRowGrid;
            nYDiff = nYMouseOnHeader - nYBorder;

            if(bBorder && Math.abs(nYDiff) <= GridRowHeader.Y_RESIZE_SENSITIVITY)
            {
                return nRow -1;
            }

            if(!bBorder && nYBorder - nHRowGrid <= nYMouseOnHeader && nYMouseOnHeader <= nYBorder)
                return nRow -1;
        }
    }

    return -1;
}
