import * as DG from "datagrok-api/dg";
import * as ui from 'datagrok-api/ui';
import {GridUtils} from "../../utils/GridUtils";
import {MCBUtils} from "./MCBUtils";
import {SemType} from "../SemType";
import ImageAdjustRowsHeights from "../../../images/adjust_rows_height.png";
import ImageFitZoom from "../../../images/fit_zoom.png";
import ImageDfaultZoom from "../../../images/default_zoom.png";
import {GridCellRenderer} from "../ui/GridCellRenderer";
import {GridColumnHeaderRenderer} from "../ui/GridColumnHeaderRenderer";
import {GridRowHeaderRenderer} from "../ui/GridRowHeaderRenderer";
import {ButtonGridColumnHeaderRenderer} from "../ui/ButtonGridColumnHeaderRenderer";
import {FiltersUtils} from "../../ui/filters/FiltersUtils";
import {FilterPanel} from "../../ui/filters/FilterPanel";
import {DGUtils} from "../../utils/DGUtils";
import {ErrorUtils} from "../../utils/ErrorUtils";

let m_handlerViewerAddedDock = null;
let m_bFiltersInProgress = false;

export class DGApp {
}

DGApp.USE_ALT_COL_NAMES = false;

DGApp.PRIMITIVE_COLS_TAG_NAME = "PRIMITIVE_COLS_TAG_NAME";
DGApp.PRIMITIVE_COL_LABEL_TAG_NAME = "PRIMITIVE_COL_LABEL_TAG_NAME";
DGApp.VIRTUAL_COL_ALT_NAME_TAG_NAME = "VIRTUAL_COL_ALT_NAME_TAG_NAME";
DGApp.VIRTUAL_PARENT_TAG_NAME = "VIRTUAL_PARENT_TAG_NAME";
DGApp.isVirtual = function(col) {

    ErrorUtils.verifyClass(col, DG.Column);
    const b = col.isVirtual;
    if(b !== true && b !== false)
     throw new Error("IsVirtual DG BUG on column " + col.name);
     /* Uncomment the following if DG continue to produce the BUG
    const arPrimiCols = col.getTag(DGApp.PRIMITIVE_COLS_TAG_NAME);
    const bV = arPrimiCols !== null && arPrimiCols !== undefined;

    if((b && !bV) || (!b && bV))
    {
      throw new Error("IsVirtual fails 1 on " + col.name);
    } */
    return b;
}

DGApp.hasPrimitiveColumns = function(col)
{
  ErrorUtils.verifyClass(col, DG.Column);
  const arColNames = col.getTag(DGApp.PRIMITIVE_COLS_TAG_NAME);
  const b = arColNames !== null && arColNames !== undefined;
  return b;
}


DGApp.hasVirtualColumns = function(dframe)
{
    let col = null;
    const nColCount = dframe.columns.length;
    for(var nCol=0; nCol<nColCount; ++nCol)
    {
        col = dframe.columns.byIndex(nCol);
        if(DGApp.isVirtual(col))
            return true;
    }

    return false;
}


DGApp.setVistualColumnsVisible = function(grid, bVisible)
{
    let colG = null;
    let col = null;
    let b = false;
    const lstCols = grid.columns;
    const nColCount = lstCols.length;
    for(let nColG= 1; nColG<nColCount; ++nColG)
      {
            colG = lstCols.byIndex(nColG);
            col = colG.column;
            b = DGApp.isVirtual(col);
            if(b)
             colG.visible = bVisible;
     }
}



DGApp.open = async function(obPackage, viewSS, configGridRenderers, arFilterClassess)
{
    const grid = viewSS.grid;
    const dframe = grid.dataFrame;

    GridUtils.setRenderersConfig(grid, configGridRenderers);

    if(arFilterClassess === undefined)
        arFilterClassess = [];

    const strURLWR = obPackage.webRoot;
               /*
    grok.events.onViewerAdded.subscribe((args) => {
        const viewer = args.args.viewer;
        const props = viewer.getProperties();
        const info = viewer.getInfo();
        const typeViewer = viewer.type;
        if(typeViewer === DG.VIEWER.GRID)
        {
            let eCanvas = viewer.canvas;
            let rrrrrr = 0;
        }
    });          */

       /*
    grok.events.onTooltipRequest.subscribe((args) => {
        let context = args.args.context;
        //if (context instanceof DG.User && context.name.includes('a'))
            args.preventDefault();
    });


    grok.events.onTooltipShown.subscribe((args) => {
        let context = args.args.context;
        let element = args.args.element;
        args.preventDefault();
        while (element.firstChild) {
            element.removeChild(element.firstChild);
        }
       // if (context instanceof DG.User)
            element.appendChild(ui.h3("Customized tooltip for"));
    });  */



    grid.onCellTooltip(function (cell, x, y) {


      return true;
    });


    GridUtils.assignSemTypes(dframe);

    //my changes const viewSS = grok.shell.addTableView(dframe);
    grid.dart.m_ctxSort = null;

    console.log("url: " + strURLWR);

    const parent = grid.root.parentElement.parentElement;//moving to the grandparent  fixes jumping of the horz splitter

    let img = new Image();
    img.src = ImageAdjustRowsHeights;

    const btnH = ui.button("H", function() {
        GridUtils.adjustRowSizeToColumnWidth(grid);
    }, "Adjust Rows Height");

    btnH.innerHTML = img.outerHTML;


    img = new Image();
    img.src = ImageFitZoom;

    const btnF = ui.button("F", function() {
        GridUtils.autoFitColumns(grid, configGridRenderers);
    }, "Zoom to Fit");

    btnF.innerHTML = img.outerHTML;

    img = new Image();
    img.src = ImageDfaultZoom;

    const btnD = ui.button("D", function() {
        GridUtils.setDefaultColumnSize(grid, configGridRenderers);
    }, "Default Zoom");

    btnD.innerHTML = img.outerHTML;

    const split = ui.splitH([btnH, btnF, btnD]);
    split.style.background = "white";
    parent.appendChild(split);

    //const headerRows = new GridRowHeader(viewSS);
    //my changes grid.columns.byIndex(0).visible = false;

    const nInsets = 2;
    const nDefaultColWidth = 75;
    const nMaxHeaderHeight = 100;

    let eCanvas =  grid.canvas;

    const g = eCanvas.getContext("2d");
    //var arLines = new Array();
    g.font = "13px Arial";

    let nHPrefHeader = GridUtils.calcPrefColHeaderHeightNew(g, dframe, configGridRenderers, nDefaultColWidth - (nInsets << 1), nInsets);
    if(nHPrefHeader > nMaxHeaderHeight)
        nHPrefHeader = nMaxHeaderHeight;

    grid.setOptions({
        colHeaderHeight: nHPrefHeader,
        //rowHeight: 150
    });


    //my changes grid.columns.rowHeader.width = 150;
    //grid.columns.byIndex(0).width = 250;


    //const nColCpd = GridUtils.findGridColumnBySemType(grid, CpdSemType);
    //const colGcpd = grid.columns.byIndex(nColCpd);  //my changes 1
    //colGcpd.visible = false;
    //GridUtils.initFilterIcons(dframe, arFilterClassess);


    GridUtils.initZoom(grid);
    GridUtils.initDefaultCellsSizes(grid, false, configGridRenderers);


    let nTopLeftCount = 0;
    let nColHeaderCount = 0;
    let nRowHeaderCount = 0;


    grid.onCellRender.subscribe(function (args) {

        try{

        const nWCell = args.bounds.width;
        const nHCell = args.bounds.height;

        //let nCol = args.cell.gridColumn.idx;
        const nRow = args.cell.gridRow;
        const nCol = args.cell.gridColumn.idx;

        const rr = args.cell.renderer;
        if (nRow < 0 && nCol === 0)//Top Left my changes !args.cell.isTableCell && !args.cell.isRowHeader && !args.cell.isRowHeader)
        {

            const renderer = configGridRenderers.getColRowHeaderRenderer();

            if (renderer !== null && !(renderer instanceof GridColumnHeaderRenderer))
                throw new Error("Column header renderer must be an instance of the " + GridColumnHeaderRenderer);

            if (renderer instanceof ButtonGridColumnHeaderRenderer) {
              renderer.setFilterEnabled(true);//FiltersUtils.hasCompatibleFilter(colAdjusted));//args.cell.tableColumn));
            }

            if (renderer !== null) {
                renderer.paint(args.g, null, args.bounds.x, args.bounds.y, nWCell, nHCell);
                args.preventDefault();
            }
            return;

            //console.log("nTopLeftCount= " + nRow + " " + nCol + " " + nTopLeftCount);
            ++nTopLeftCount;
            //args.preventDefault();
        } else if (nRow < 0 && nCol > 0)//my changes args.cell.isColHeader)
        {
            let tp = args.cell.cellType;
            let nCol = args.cell.gridColumn.idx;
            let nRow = args.cell.gridRow;
            if (nRow !== -1) {
                throw new Error("Row must be -1.")
            }


            //let ctxRenderer = new SemEntityRendererContext(args.cell.cell, args.cell.column, dframe);
            let renderer = null;

            if (GridUtils.getSemType(args.cell.tableColumn) instanceof SemType) {
                const typeSem = GridUtils.getSemType(args.cell.tableColumn);
                renderer = configGridRenderers.getColHeaderRenderer(typeSem.constructor);
            } else renderer = configGridRenderers.getDefaultColHeaderRenderer();


            if (renderer === undefined || renderer === null) {
                let typeSem = GridUtils.getSemType(args.cell.tableColumn);
                renderer = configGridRenderers.getColHeaderRenderer(typeSem.constructor);

                throw new Error("Renderer cannot be null.")
            }

            if (!(renderer instanceof GridColumnHeaderRenderer))
                throw new Error("Column header renderer must be an instance of the " + GridColumnHeaderRenderer);

            if (renderer instanceof ButtonGridColumnHeaderRenderer) {
                const colAdjusted = renderer.adjustColumn(args.cell);
                renderer.setFilterEnabled(true);//FiltersUtils.hasCompatibleFilter(colAdjusted));//args.cell.tableColumn));
            }

            renderer.setFont(args.cell.style.font);
            const b = renderer.paint(args.g, args.cell.gridColumn, args.bounds.x, args.bounds.y, nWCell, nHCell);
            if (b) {
                args.preventDefault();
                ++nColHeaderCount;
                return;
            }


        }
        else if (args.cell.isRowHeader)//args.cell.cellType === "row header")
        {
              if (args.cell.cell !== null) {
                const rendererRowHeader = configGridRenderers.getRowHeaderRenderer();
                if (rendererRowHeader !== null) {

                    if(!(rendererRowHeader instanceof GridRowHeaderRenderer))
                      throw new Error("Column header renderer must be an instance of the " + GridRowHeaderRenderer);

                    const nRowGrid = args.cell.gridRow;
                    const nRowTable = args.cell.tableRowIndex;

                    rendererRowHeader.paint(args.g, nRowGrid, nRowTable, args.cell.grid, args.bounds.x, args.bounds.y, nWCell, nHCell);
                    //rendererRowHeader.paint(gRowHeader, nRowGrid, nRowTable, args.cell.grid, 0, args.bounds.y, eDivRowHeader.offsetWidth, nHCell);

                    args.preventDefault();
                }


                //console.log("nRowHeaderCount= " + nRowHeaderCount);
                ++nRowHeaderCount;
            } else {
                let rty = 0;
            }

        } else if (args.cell.isTableCell) {
            const nCol = args.cell.gridColumn.idx;
            const nRow = args.cell.gridRow;
            const nRecord = args.cell.tableRowIndex;
            const val = args.cell.cell.value;

            if (GridUtils.getSemType(args.cell.tableColumn) instanceof SemType) {
                const dframe = args.cell.cell.dataFrame;

                args.g.fillStyle = "white";
                args.g.fillRect(args.bounds.x, args.bounds.y, nWCell, nHCell);

                let b = false;
                if (val !== null) {
                    const typeSem = GridUtils.getSemType(args.cell.tableColumn);
                    const renderer = configGridRenderers.getCellRenderer(typeSem.constructor);

                    if (!(renderer instanceof GridCellRenderer))
                        throw new Error("Column header renderer must be an instance of the " + GridCellRenderer);


                    try {
                        b = renderer.paint(args.g, args.cell, args.bounds.x, args.bounds.y, nWCell, nHCell);
                    } catch (e) {
                        throw new Error(e);
                    }

                }
                if (!b)
                    args.preventDefault();
            }


         //External Row Header
            /*
                const headerRowsTmp = GridUtils.getRowsHeader(grid);
                if(headerRowsTmp !== null && GridUtils.isFirstViewportColumn(args.cell.gridColumn, args.cell.grid))
                {
                    const gRowHeaderTmp = headerRowsTmp.getGraphics();
                    const nRowGrid = args.cell.gridRow;
                    headerRowsTmp.paintRow(gRowHeaderTmp, nRowGrid, grid);
                 }
              */                                                                               ``
        } else {
            let ddd = 0;
        }

    }
    catch(e)
    {
        throw e;
    }
    });


    grok.events.onContextMenu.subscribe((args) => {

        try{console.log("type: " + args.args.context.dart.m_fZoomValue !== undefined);}
        catch(e)
        {
            let rrr = 0;
        }

        if (args.args.context instanceof DG.Viewer && args.args.context.dart.m_fZoomValue !== undefined)
        {
            //if (view.type === DG.VIEW_TYPE.TABLE_VIEW && view.name === 'Main')
            //this.layoutMain(view);

            const e = args.causedBy;
            const gridTmp = args.args.context;
            const cell = grid.hitTest(e.offsetX, e.offsetY);
            if(cell === undefined || cell.cellType === null)//bufg in DG , top left cell
                return;



            const gridd = cell.grid;
            const colGrid = cell.gridColumn;
            let bRH = cell.isRowHeader;
            let bCH = cell.isColHeader;

            let menu = args.args.menu;
            let popup = menu;

            menu.clear();

            function zoom(strItem)
            {
                if(strItem === "Default")
                    GridUtils.setDefaultColumnSize(cell.grid, configGridRenderers);
                else if(strItem === "To Fit")
                    GridUtils.autoFitColumns(cell.grid, configGridRenderers);
                else if(strItem === "Compound Struct")
                    GridUtils.adjustRowSizeToColumnWidth(cell.grid);

            }
            menu.group("Zoom").items(["Default", "To Fit", "Compound Struct"], zoom);

            let renderer = null;
            if(colGrid.idx === 0)
              renderer = configGridRenderers.getRowHeaderRenderer();
            else if(bCH)
            {
             const colT = cell.tableColumn;
             if(GridUtils.getSemType(colT) instanceof SemType)
              renderer = configGridRenderers.getColHeaderRenderer(GridUtils.getSemType(colT).constructor);
            }
            else if(cell.isTableCell)
            {
                const colT = cell.tableColumn;
                const typeSem = GridUtils.getSemType(colT);
                if(typeSem instanceof SemType)
                 renderer = configGridRenderers.getCellRenderer(typeSem.constructor);
            }

            if(renderer !== null)
            {
                const nGroupCount = renderer.getPopupGroupCount(cell);
                if(nGroupCount > 0)
                {

                    let strGroupName = "";
                    let arItems = [];
                    for(var nGroup = 0; nGroup<nGroupCount; ++nGroup)
                    {
                        menu = menu.separator();
                        strGroupName = renderer.getPopupGroupName(cell, nGroup);
                        arItems = renderer.toPopupGrouptemsArray(cell, nGroup);
                        let nGroupTmp = nGroup;
                        let arItemsTmp = arItems;
                        menu.group(strGroupName).items(arItems, function(strItem){

                            const nItemClicked = arItemsTmp.indexOf(strItem);
                            renderer.onPopupGroupAction(cell, nGroupTmp, nItemClicked);
                        });
                    }
                }

             const arItems = renderer.toPopupItemsArray(cell);
             if(arItems.length > 0) {

                menu = menu.separator();

                function onItem(strItem) {
                    const nItemClicked = arItems.indexOf(strItem);
                    renderer.onPopupAction(cell, nItemClicked);
                }

                menu.items(arItems, onItem);
              }
            }
         }
    });


    document.addEventListener("wheel", function (e) {

        if(!DGUtils.isCurrentView(viewSS))
            return;


        let bCtrl =e.shiftKey;
        if(!bCtrl)
            return;


        let opt = grid.getOptions();

        let nXScroll1 = grid.horzScroll.min;
        let nYScroll2 = grid.vertScroll.min;
        let nXScroll = grid.horzScroll.min;
        let nYScroll = grid.vertScroll.min*GridUtils.getRowHeight(grid);

        let cellGrid = grid.hitTest(e.offsetX, e.offsetY);

        let nWColOld = GridUtils.getColumnWidth(cellGrid.gridColumn);
        let nHColOld=  GridUtils.getRowHeight(grid);
        let nXOld =  GridUtils.getColumnX(cellGrid.gridColumn);
        let nYOld =  cellGrid.gridRow*GridUtils.getRowHeight(grid);
        let nXMouse = nXScroll + e.offsetX;
        let nYMouse = nYScroll + e.offsetY;

        if(bCtrl)
        {
            let nWheelLengthTillFull = 19;

            let fMinimum = grid.dart.getMinimumZoomValue();
            let fMaximum = grid.dart.getMaximumZoomValue();
            let fK = (fMaximum - fMinimum)/nWheelLengthTillFull;
            fK = 0.2;

            let nCount = -1;
            let fZoomValue = 1.0;

            if(e.deltaY > 0)//zoom in
            {
                if(bCtrl)
                    fK = GridUtils.getZoomValue(grid)*0.05;

                nCount = Math.floor((GridUtils.getZoomValue(grid) - grid.dart.getMinimumZoomValue())/fK);
                fZoomValue = nCount === 0 ? grid.dart.getMinimumZoomValue() : GridUtils.getZoomValue(grid)- 1*fK;

                if(fZoomValue < grid.dart.getMinimumZoomValue())
                    fZoomValue = grid.dart.getMinimumZoomValue();

            }
            else
            {
                if(bCtrl)
                    fK = GridUtils.getZoomValue(grid)*0.05;

                nCount = Math.floor((grid.dart.getMaximumZoomValue() - grid.dart.getZoomValue())/fK);
                fZoomValue = nCount === 0 ? grid.dart.getMaximumZoomValue() : grid.dart.getZoomValue()+ 1*fK;

                if(fZoomValue > grid.dart.getMaximumZoomValue())
                    fZoomValue = grid.dart.getMaximumZoomValue();

            }
            //console.log("Zoom: " + fZoomValue);
            GridUtils.setZoomValue(grid, fZoomValue);

            let nWColNew =  GridUtils.getColumnWidth(cellGrid.gridColumn);
            let nHColNew=   GridUtils.getRowHeight(grid);
            let nXNew = GridUtils.getColumnX(cellGrid.gridColumn);
            let nYNew = cellGrid.gridRow*GridUtils.getRowHeight(grid);
            let fKX = nWColNew/nWColOld;
            let fKY = nHColNew/nHColOld;

            let nX = nXNew + Math.floor(fKX*(nXMouse - nXOld));//nXNew + Math.floor(fKX*(e.offsetX - nXOld));
            let nY = nYNew + Math.floor(fKY*(nYMouse - nYOld));

            grid.scrollToPixels(nXScroll + nX - nXMouse, nYScroll + nY - nYMouse);
            // e.preventDefault();


            e.stopPropagation();
        }


    }, true);



    rxjs.fromEvent(grid.overlay, 'mousedown').subscribe((e) => {

        const nButton = e.button;
        const cell = grid.hitTest(e.offsetX, e.offsetY);
        if(cell === undefined || cell.dart === undefined)
            return;//DDG Bug


        const colG = cell.gridColumn;
        const bCH = cell.isColHeader;

        if(bCH)
        {
            let renderer = null;
            if(colG.idx === 0)
             renderer = configGridRenderers.getColRowHeaderRenderer();
             else
             {
                const colT = cell.tableColumn;
                if(GridUtils.getSemType(colT) instanceof SemType)
                 renderer = configGridRenderers.getColHeaderRenderer(GridUtils.getSemType(colT).constructor);
                else renderer = configGridRenderers.getDefaultColHeaderRenderer();
             }

             if(renderer !== null)
             {
                 const strName = colG.name;
                 const sbHorz = cell.grid.horzScroll;
                 const nXOnCell = e.offsetX - GridUtils.getColumnX(colG) + sbHorz.min;
                 const nYOnCell = e.offsetY - 0;

                 renderer.onMousePressed(cell, nXOnCell, nYOnCell, nButton);
             }
        }
    });


    rxjs.fromEvent(grid.overlay, 'mouseup').subscribe((e) => {

        const nButton = e.button;
        const cell = grid.hitTest(e.offsetX, e.offsetY);

        if(cell === undefined || cell.dart === undefined)
            return;//DDG Bug

        const colG = cell.gridColumn;
        const bTC = cell.isTableCell;
        const bRH = cell.isRowHeader;
        const bCH = cell.isColHeader;

        if(bCH)
        {
            let renderer = null;
            if(colG.idx === 0)
                renderer = configGridRenderers.getColRowHeaderRenderer();
            else
            {
                const colT = cell.tableColumn;
                if(GridUtils.getSemType(colT) instanceof SemType)
                   renderer = configGridRenderers.getColHeaderRenderer(GridUtils.getSemType(colT).constructor);
                else renderer = configGridRenderers.getDefaultColHeaderRenderer();
            }

            if(renderer !== null)
            {
                const strName = colG.name;
                const sbHorz = cell.grid.horzScroll;
                const nXOnCell = e.offsetX - GridUtils.getColumnX(colG) + sbHorz.min;
                const nYOnCell = e.offsetY - 0;

                renderer.onMouseReleased(cell, nXOnCell, nYOnCell, nButton);
            }
        }
    });




    rxjs.fromEvent( grid.overlay, 'click').subscribe(async (e) =>  {

        const nButton = e.button;
        const cell = grid.hitTest(e.offsetX, e.offsetY);
        const colG = cell.gridColumn;
        const bCH = cell.isColHeader;
        if(bCH)
        {
            let renderer = null;
            if(colG.idx === 0)
                renderer = configGridRenderers.getColRowHeaderRenderer();
            else
            {
                const colT = cell.tableColumn;
                if(GridUtils.getSemType(colT) instanceof SemType)
                    renderer = configGridRenderers.getColHeaderRenderer(GridUtils.getSemType(colT).constructor);
                else renderer = configGridRenderers.getDefaultColHeaderRenderer();
            }

            if(renderer !== null)
            {
                const strName = colG.name;
                const sbHorz = cell.grid.horzScroll;

                const nXOnCell = e.offsetX - GridUtils.getColumnX(colG) + sbHorz.min;
                const nYOnCell = e.offsetY - 0;

                renderer.onMouseClicked(cell, nXOnCell, nYOnCell, nButton);
            }
       }

        e.preventDefault();
        e.stopPropagation();
    });

    let nXMouseOnGrid = NaN;
    let nYMouseOnGrid = NaN;
    rxjs.fromEvent( document, 'mousemove').subscribe(async (e) => {

        const rect = grid.overlay.getBoundingClientRect();
        nXMouseOnGrid = e.clientX - rect.left; //x position within the element.
        nYMouseOnGrid = e.clientY - rect.top;  //y position within the element.

        if(!document.hasFocus())
            grid.overlay.focus();

         //console.log(x + " " + y);
    });


    rxjs.fromEvent( grid.overlay, 'mousemove').subscribe(async (e) =>  {

        const nButton = e.button;
        let cell = null;

        try{cell = grid.hitTest(e.offsetX, e.offsetY);}
        catch(e)
        {
            //Huge DG bug
            return;
        }

        if(cell === undefined || cell.dart === undefined)
            return; // to address bug in DG when the all the cells are filtered out cell should be null, not cell.d === undefined

        const colG = cell.gridColumn;
        const bCH = cell.isColHeader;
        if(bCH)
        {
            const bTooltipOn = GridUtils.isTooltipEnabled(grid);
            if(e.ctrlKey || bTooltipOn)  //Tooltip
            {
                const rect = grid.overlay.getBoundingClientRect();
                const scrollLeft= window.pageXOffset || document.documentElement.scrollLeft;
                const scrollTop = window.pageYOffset || document.documentElement.scrollTop;
                const nY = rect.top  + scrollTop;
                const nX = rect.left + scrollLeft;

                const nXTT = cell.bounds.x + nX;
                const nYTT = cell.bounds.y + nY;

                const typeSem = cell.tableColumn === null ? null : GridUtils.getSemType(cell.tableColumn);
                if(typeSem instanceof SemType)
                {
                    const renderer = configGridRenderers.getColHeaderRenderer(typeSem.constructor);
                    const eTooltip = renderer.createTootipContent(cell);
                    ui.tooltip.show(ui.divV([
                        eTooltip
                    ]), nXTT, nYTT);
                }
                else ui.tooltip.hide();
                return;
            }

            //Regular Move
            let renderer = null;
            if(colG.idx === 0)
                renderer = configGridRenderers.getColRowHeaderRenderer();
            else
            {
                const colT = cell.tableColumn;
                if(GridUtils.getSemType(colT) instanceof SemType)
                    renderer = configGridRenderers.getColHeaderRenderer(GridUtils.getSemType(colT).constructor);
                else renderer = configGridRenderers.getDefaultColHeaderRenderer();
            }

            if(renderer !== null)
            {
                if(renderer instanceof ButtonGridColumnHeaderRenderer) {
                    const colAdjusted = renderer.adjustColumn(cell);
                    renderer.setFilterEnabled(true); //FiltersUtils.hasCompatibleFilter(colAdjusted));
                }


                const strName = colG.name;
                const sbHorz = cell.grid.horzScroll;
                const nXOnCell = e.offsetX - GridUtils.getColumnX(colG) + sbHorz.min;
                const nYOnCell = e.offsetY - 0;

                renderer.onMouseMoved(cell, nXOnCell, nYOnCell, nButton);
            }
        }

        else if(cell.isTableCell)
        {
          const bTooltipOn = GridUtils.isTooltipEnabled(grid);
          if(e.ctrlKey || bTooltipOn)
          {
              const rect = grid.overlay.getBoundingClientRect();
              const scrollLeft= window.pageXOffset || document.documentElement.scrollLeft;
              const scrollTop = window.pageYOffset || document.documentElement.scrollTop;
              const nY = rect.top  + scrollTop;
              const nX = rect.left + scrollLeft;

              const nXOnCell = cell.bounds.x + cell.gridColumn.width + nX;;
              const nYOnCell = cell.bounds.y + nY;

              const typeSem = GridUtils.getSemType(cell.tableColumn);
              if(typeSem instanceof SemType)
              {
                  const renderer = configGridRenderers.getCellRenderer(typeSem.constructor);
                  const eTooltip = renderer.createTootipContent(cell);
                  ui.tooltip.show(ui.divV([
                      eTooltip
                  ]), nXOnCell, nYOnCell);
              }
              else ui.tooltip.hide();
          }
        }
        else  ui.tooltip.hide();
    });


    rxjs.fromEvent(document, 'keydown').subscribe(async (e) =>  {

        if(!DGUtils.isCurrentView(viewSS))
            return;

        if(0 <= nXMouseOnGrid && nYMouseOnGrid < grid.overlay.width)
        {

            let cell = null;
            try{cell = grid.hitTest(nXMouseOnGrid, nYMouseOnGrid);}
            catch(e)
            {
                //Huge DG bug
                return;
            }

            if(cell === undefined || cell.dart === undefined || (!cell.isTableCell && !cell.isColHeader)) {
                ui.tooltip.hide();
                return;
            }

            const rect = grid.overlay.getBoundingClientRect();
            const scrollLeft= window.pageXOffset || document.documentElement.scrollLeft;
            const scrollTop = window.pageYOffset || document.documentElement.scrollTop;
            const nY = rect.top  + scrollTop;   //y-coord of the grid overlay in the document's coordinate space
            const nX = rect.left + scrollLeft;  //x-coord of the grid overlay in the document's coordinate space

            const nXTT = cell.bounds.x + nX + (cell.isTableCell ? cell.gridColumn.width : 0);
            const nYTT = cell.bounds.y + nY;

            const typeSem = cell.tableColumn === null ? null : GridUtils.getSemType(cell.tableColumn);
            if(typeSem instanceof SemType)
            {
                const renderer = cell.isTableCell ? configGridRenderers.getCellRenderer(typeSem.constructor) :
                                                    configGridRenderers.getColHeaderRenderer(typeSem.constructor);
                const eTooltip = renderer.createTootipContent(cell);
                ui.tooltip.show(ui.divV([
                    eTooltip
                ]), nXTT, nYTT);
            }
       }
    });


    rxjs.fromEvent( document, 'keyup').subscribe(async (e) =>  {

        if(!DGUtils.isCurrentView(viewSS))
            return;

        if(e.code === "ControlLeft")
         ui.tooltip.hide();
    });


    grid.onColumnResized.subscribe((ev) => {

        const strColName = ev.args.column.name;
        const bAdjusting = ev.args.dragging;
        const colGrid = grid.columns.byName(strColName);
        const col = colGrid.column;
        if(col === null)
            return;

        const typeSem = GridUtils.getSemType(col);
        if(typeSem instanceof SemType)
        {
            const renderer = configGridRenderers.getCellRenderer(typeSem.constructor);
            const nW = colGrid.width;
            renderer.onResizeWidth(colGrid, grid, nW, bAdjusting);
        }

        //grok.shell.info(`Resizing ${ev.args.column.name}: ` + (ev.args.dragging ? "in progress" : "done"));
    });


    grid.onRowsResized.subscribe((ev) => {

        const bAdjusting = ev.args.dragging;

        const arColRowIdxs = new Array(4);
        GridUtils.fillVisibleGridCells(arColRowIdxs, grid);

        const nHRow = GridUtils.getRowHeight(grid);

        const nColMin = arColRowIdxs[0];
        const nColMax = arColRowIdxs[1];

        //const nRowMin = arColRowIdxs[2];
        //const nRowMax = arColRowIdxs[3];
        let colGrid = null;
        let col = null;
        let typeSem = null;
        let renderer= null;

        //console.log("Requesting cols " + nColMin + " " + nColMax + " for " + nHRow);
        let nColVisited = 0;
        for(var nCol=nColMin; nCol<=nColMax; ++nCol)
        {
            if(nCol === 0)
             continue;

            colGrid = grid.columns.byIndex(nCol);
            if(!colGrid.visible)
                continue;

            col = colGrid.column;
            typeSem = GridUtils.getSemType(col);

            if(!(typeSem instanceof SemType))
                continue;

            renderer = configGridRenderers.getCellRenderer(typeSem.constructor);
            if(renderer === null)
                continue;

            renderer.onResizeHeight(colGrid, grid, nHRow, bAdjusting);
            ++nColVisited;
        }


       //grok.shell.info("Resizing row height: " + (ev.args.dragging ? "in progress" : "done"));
    });


       /*
    rxjs.fromEvent( grid.overlay, 'dblclick').subscribe(async (e) =>  {

        e.preventDefault();
        e.stopPropagation();
    });  */
    grok.events.onViewerClosed.subscribe((args) => {
        const viewer = args.args.viewer;

        const type = viewer.type;
        const root = viewer.root;

        if(type === DG.VIEWER.FILTERS && viewer.dart.m_panelFilters !== undefined && viewer.dart.m_panelFilters !== null)
        {
            viewer.dart.m_panelFilters.detach();
           //root.removeChild(viewer.d.m_panelFilters.root); root no longer available

            viewer.dart.m_panelFilters = null;
        }
    });


    //ViewerLayoutManager.install(new FlowTileLayoutManager(2));
    //ViewerLayoutManager.uninstall();
    //ViewerLayoutManager.install(new FlowTileLayoutManager(2));


    grok.events.onViewerAdded.subscribe((args) => {

        const viewer = args.args.viewer;

        const b = DGUtils.containsViewer(viewSS, viewer);
        if(!b)
          return;

        //DGApp.setVistualColumnsVisible(viewSS.grid, true);

        const type = viewer.type;
        const root = viewer.root;

        if(type === DG.VIEWER.FILTERS)
        {
         FiltersUtils.showFilters(false);

           // const filter = compoundFilter();
            //filter.attach(viewer.dataFrame);
            //FiltersUtils.addFilter(filter, viewer);

            viewer.dart.m_panelFilters = null;
         if(arFilterClassess !== undefined)
         {
             if(viewer.dart.m_panelFilters === null) {
                 const panelFilters = new FilterPanel(arFilterClassess);
                 root.appendChild(panelFilters.root);
                 panelFilters.attach(viewer.dataFrame);
                 viewer.dart.m_panelFilters = panelFilters;
             }
         }
        }


        viewer.onEvent('d4-column-combo-box-popup-show').subscribe((args) => {

            const viewer = args.args.viewer;
            const eDivTree = ui.div();
            eDivTree.id = "datajtree";
            viewer.root.appendChild(eDivTree);

            //alert("Column Chooser: " + args.args.selectorName);

            const strselectormName = args.args.selectorProperty.name;
            const strPropType = args.args.selectorProperty.propertyType;
            const strDescr = args.args.selectorProperty.description;

            let strTitle = "";
            if(strselectormName === "xColumnName")
              strTitle = "X-Axis Column Selector";

            else if(strselectormName === "yColumnName")
                strTitle = "Y-Axis Column Selector";

            else if(strselectormName === "colorColumnName") {
                strTitle = "Color Column Selector";
            }

            else if(strselectormName === "sizeColumnName") {
                strTitle = "Size Column Selector";
            }

                //}
                //else if(viewer.type === DG.VIEWER.HISTOGRAM)
            //{
            else if(strselectormName === "valueColumnName") {
                strTitle = "Value Column Selector";
            }

            else if(astrselectormName === "splitColumnName") {
               strTitle = "Split Column Selector";
            }

            else if(strselectormName === "categoryColumnName") {
               strTitle = "Category Column Selector";
            }


            if(DGApp.m_dlgColSelector != null)
                DGApp.m_dlgColSelector.close();

            const dial = ui.dialog(strTitle);
            DGApp.m_dlgColSelector = dial;

            const cbb = args.args.comboBox;



            let treedata = null;
            try {
                treedata = MCBUtils.createColumnSelector(eDivTree, viewer.dataFrame,function(strColName, bRollover){

                    //if(viewer.type.toLowerCase() === DG.VIEWER.SCATTER_PLOT)
                    //{

                    if(strselectormName === "xColumnName") {
                        viewer.setOptions({x: strColName});

                    }
                    else if(strselectormName === "yColumnName") {
                        viewer.setOptions({y: strColName});
                    }

                    else if(strselectormName === "colorColumnName") {
                        viewer.setOptions({colorColumnName: strColName});
                    }

                    else if(strselectormName === "sizeColumnName") {
                        viewer.setOptions({sizeColumnName: strColName});
                    }

                        //}
                        //else if(viewer.type === DG.VIEWER.HISTOGRAM)
                    //{
                    else if(strselectormName === "valueColumnName") {
                        viewer.setOptions({valueColumnName: strColName});
                    }

                    else if(astrselectormName === "splitColumnName") {
                        viewer.setOptions({splitColumnName: strColName});
                    }

                    else if(strselectormName === "categoryColumnName") {
                        viewer.setOptions({splitColumnName: strColName});
                    }

                    //}

                    if(!bRollover)
                        dial.close();
                });
            }
            catch(e)
            {
                console.log("ERROR: Failed column selector")
               return;
            }

            //ui.showPopup(eDiv, args.args.comboBox.root, args.args.comboBox.vertical);
            const eLabel = ui.divText("");
            eLabel.innerHTML = "<b>Total Columns: "+ treedata.length + "</b>";

            const eEdit = ui.stringInput("search", "");
            eEdit.root.addEventListener('keyup', function(e){

                const source = e.currentTarget;
                const strText = eEdit.value;

                const arNodes = MCBUtils.createColumnSelectorData(dframe, strText);
                eLabel.innerHTML = "<b>Total Columns: "+ arNodes.length + "</b>";

                // console.log("text: " + strText + " node_count: "  + arNodes.length);

                MCBUtils.updateColumnSelector(eDivTree, arNodes);
             });



            //let eDivPanel = ui.splitV([eEdit.root, eLabel, eDivTree]);

            dial.add(eEdit.root);
            dial.add(eLabel);
            dial.add(eDivTree);

            const bodyRect = document.body.getBoundingClientRect();
            const viewerRect = viewer.root.getBoundingClientRect();
            const nX = viewerRect.left - bodyRect.left;
            const nY = viewerRect.top - bodyRect.top;

            dial.show({modal: false, x: nX + 0.1*viewerRect.width, y: nY + 0.1*viewerRect.height, width: 0.8*viewerRect.width, height: 0.8*viewerRect.height});
            args.preventDefault();
        });

    });

    return viewSS;
}

DGApp.m_dlgColSelector = null;
