import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as GridUtils from '../utils/GridUtils';
import * as TextUtils from '../utils/TextUtils';
import {ColorUtils} from '../utils/ColorUtils';
import * as rxjs from 'rxjs';
import { GridCellRendererEx} from "../renderer/GridCellRendererEx";
import * as PinnedUtils from "./PinnedUtils";
import {TableView} from "datagrok-api/dg";
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

function getRenderer(cell : DG.GridCell) : GridCellRendererEx | DG.GridCellRenderer {
  const colGrid = cell.gridColumn;
  if (colGrid === null || colGrid === undefined) {
    throw new Error('Grid cell is detached from the Grid column');
  }

  let renderer = GridUtils.getGridColumnRenderer(colGrid);
  if(renderer instanceof GridCellRendererEx) {
    return renderer;
  }

  return cell.renderer;
}


function getGrid(colGrid : DG.GridColumn) : DG.Grid | null {
  let grid : DG.Grid | null = colGrid.grid;
  if( grid === null) {
    grid = GridUtils.getInstalledGridForColumn(colGrid);
    if(grid instanceof DG.Grid)
      return grid;
  }

  return grid;
}


function notifyAllColsRowsResized(grid : DG.Grid, nHRows : number, bAdjusting : boolean) : void {

  let renderer : GridCellRendererEx | null = null
  let colGrid = null;
  const lstColsGrid = grid.columns;
  const nColCount = lstColsGrid.length;
  for(let nCol=0; nCol<nColCount; ++nCol) {
    colGrid = lstColsGrid.byIndex(nCol);
    if(colGrid === null || !colGrid.visible){
      continue
    }

    renderer = GridUtils.getGridColumnRenderer(colGrid);
    if (renderer instanceof GridCellRendererEx) {
      renderer.onResizeHeight(colGrid, grid, nHRows, bAdjusting);
    }
  }
}


function notifyAllPinnedColsRowsResized(colPinnedSource : PinnedColumn, nHRows : number, bAdjusting : boolean) : void {

  const colGridSource  = colPinnedSource.getGridColumn();
  if(colGridSource === null){
    return;
  }

  const grid = colGridSource.grid;
  const dart = DG.toDart(grid);
  if(dart.m_arRowHeaders === undefined) {
    throw new Error('Pinned Columns are not installed.');
  }

  let renderer : GridCellRendererEx | null = null
  let colPinned = null;
  let colGrid = null;
  const nPinnedColCount = dart.m_arRowHeaders.length;
  for(let nColPin=0; nColPin<nPinnedColCount; ++nColPin) {
    colPinned = dart.m_arRowHeaders[nColPin];
    colGrid = colPinned.m_colGrid;
    if(colGrid === null) {
      throw new Error('Pinned Column is detached.');
    }

    renderer = GridUtils.getGridColumnRenderer(colGrid);
    if (renderer instanceof GridCellRendererEx  && colPinned.m_root !== null) {
      renderer.onResizeHeight(colPinned, grid, nHRows, bAdjusting);
    }
  }
}




export class PinnedColumn {

  private static MIN_ROW_HEIGHT = 20;
  private static MAX_ROW_HEIGHT = 500;
  private static SELECTION_COLOR = ColorUtils.toRgb(ColorUtils.colSelection); //"rgba(237, 220, 88, 0.15)";
  private static ACTIVE_CELL_COLOR = ColorUtils.toRgb(ColorUtils.currentRow); //"rgba(153, 237, 82, 0.25)";
  private static Y_RESIZE_SENSITIVITY = 2;

  private m_colGrid : DG.GridColumn | null;
  private m_root : HTMLCanvasElement | null;
  //private m_observerResize : ResizeObserver | null;
  private m_observerResizeGrid : ResizeObserver | null;
  private m_handlerVScroll : any;
  private m_handlerRowsFiltering : any;
  private m_handlerCurrRow : any;
  private m_handlerSel : any;
  //private m_handlerFilter : any;
  private m_handlerRowsResized : any;
  private m_handlerRowsSorted : any;
  private m_handlerMouseDown : any;
  private m_handlerMouseUp : any;
  private m_handlerMouseLeave : any;
  private m_handlerMouseMove : any;

  constructor(colGrid : DG.GridColumn) {

    const grid = getGrid(colGrid);
    if(grid === null) {
      throw new Error("Column '" + colGrid.name + "' is not attached to the grid.");
    }

    if(!PinnedUtils.isPinnableColumn(colGrid)) {
      throw new Error("Column '" + colGrid.name + "' cannot be pinned. It either pinned or HTML.");
    }

    const dart = DG.toDart(grid);

    if(dart.m_arRowHeaders === undefined)
      dart.m_arRowHeaders = [];

    dart.m_arRowHeaders.push(this);

    const viewTable = grid.view;
    const dframe = grid.dataFrame;

    const nW = colGrid.width;
    this.m_colGrid = colGrid;

    try {
      colGrid.visible = false;
    }
    catch(e) {
      //DG bug
      console.error("ERROR: Couldn't hide column.");
    }

    if(colGrid.settings === null || colGrid.settings === undefined)
      colGrid.settings = {};

    colGrid.settings.isPinned = true; //this will be saved with the layout
    colGrid.settings.idxPinned = dart.m_arRowHeaders.length -1;

    grid.canvas.style.left = (grid.canvas.offsetLeft + nW).toString() + "px";
    grid.overlay.style.left= (grid.overlay.offsetLeft + nW).toString() + "px";

    grid.canvas.style.width = (grid.canvas.offsetWidth - nW).toString() + "px";
    grid.overlay.style.width= (grid.overlay.offsetWidth - nW).toString() + "px";

    const eCanvasThis = ui.canvas(nW, 2500);

    if(grid.canvas.parentNode === null)
      throw new Error("Parent node for canvas cannot bre null.");

    grid.canvas.parentNode.insertBefore(eCanvasThis, grid.canvas);
    this.m_root = eCanvasThis;

    const colGrid0 = grid.columns.byIndex(0);
    if(colGrid0 !== null && colGrid0 !== undefined) {//DG Bug from reading layout
    try{
        colGrid0.visible = false;
      }
      catch(e) {
        console.error("ERROR: Couldn't hide row header.");
      }
    }



    //OnResize Row header
    const headerThis = this;/*
    this.m_observerResize = new ResizeObserver(entries => {
      const g = headerThis.m_root.getContext('2d');
      for (let entry of entries) {
        headerThis.paint(g, grid);
      }
    });



    this.m_observerResize.observe(headerThis.m_root);*/

    //OnResize Row header
    this.m_observerResizeGrid = new ResizeObserver(entries => {

      const bCurrent =  DG.toDart(grok.shell.v) === DG.toDart(viewTable);
      if(!bCurrent)
        return;

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
        headerThis.paint(g, grid);
      }
    });

    this.m_observerResizeGrid.observe(grid.root);

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
      }
    );

    this.m_handlerSel = dframe.onSelectionChanged.subscribe((e : any) => {
        const g = eCanvasThis.getContext('2d');
        headerThis.paint(g, grid);
      }
    );

/*
    this.m_handlerFilter = dframe.onRowsFiltered.subscribe((e : any) => {
        const g = eCanvasThis.getContext('2d');
        headerThis.paint(g, grid);
      }
    );
*/

    this.m_handlerRowsResized = grid.onRowsResized.subscribe((e : any) => {
        const g = eCanvasThis.getContext('2d');
        headerThis.paint(g, grid);
      }
    );

    this.m_handlerRowsSorted = grid.onRowsSorted.subscribe((e : any) => {
        const g = eCanvasThis.getContext('2d');
        headerThis.paint(g, grid);
      }
    );

    let nHResizeRowsBeforeDrag  = -1;
    let nResizeRowGridDragging = -1;
    let nYResizeDraggingAnchor = -1;
    let nResizeRowGridMoving = -1;

    let nYDraggingAnchor = -1;
    let nRowGridDragging = -1;

    let arXYMouseOnCellDown = [-2, -2];
    let arXYMouseOnCellUp = [-1, -1];

    this.m_handlerMouseDown = rxjs.fromEvent(document, 'mousedown').subscribe((e : Event) => {

      if(DG.toDart(grok.shell.v) !== DG.toDart(viewTable))
        return;

      const eMouse = e as MouseEvent;

      if(eMouse.buttons !== 1)
        return;

      nResizeRowGridMoving = -1;

      const bAddToSel : boolean = eMouse.ctrlKey || eMouse.shiftKey;

      let nRowGrid = bAddToSel ? -1 : PinnedColumn.hitTestRows(eCanvasThis, grid, eMouse, true, undefined);
      if (nRowGrid >= 0) {
        const nHRows = GridUtils.getGridRowHeight(grid);
        nResizeRowGridDragging = nRowGrid;
        nYResizeDraggingAnchor = eMouse.clientY;
        nHResizeRowsBeforeDrag = nHRows;
      }
      else
      {
        nRowGrid = PinnedColumn.hitTestRows(eCanvasThis, grid, eMouse, false, arXYMouseOnCellDown);

        nRowGridDragging = nRowGrid;
        nYDraggingAnchor = eMouse.clientY;

        const cell = grid.cell(colGrid.name, nRowGrid);
        const renderer = getRenderer(cell);
        if(renderer instanceof GridCellRendererEx) {
          renderer.onMouseDownEx(cell, eMouse, arXYMouseOnCellDown[0], arXYMouseOnCellDown[1]);
        }
      }

      e.preventDefault();
      e.stopPropagation();
    });


    this.m_handlerMouseUp = rxjs.fromEvent(document, 'mouseup').subscribe((e) => {

      if(DG.toDart(grok.shell.v) !== DG.toDart(viewTable)) {
        return;
      }

      const eMouse = e as MouseEvent;

      if(nResizeRowGridDragging >= 0) {
        const nHRow = GridUtils.getGridRowHeight(grid);
        notifyAllPinnedColsRowsResized(headerThis, nHRow, false);
        notifyAllColsRowsResized(grid, nHRow, false);
      }


      nHResizeRowsBeforeDrag = -1;
      nResizeRowGridDragging = -1;
      nYResizeDraggingAnchor = -1;
      nResizeRowGridMoving = -1;

      document.body.style.cursor = "auto";

      if(nRowGridDragging >= 0) {
        const bAddToSel = eMouse.ctrlKey || eMouse.shiftKey;

        const nRowGrid = PinnedColumn.hitTestRows(eCanvasThis, grid, eMouse, false, arXYMouseOnCellUp);
        if(!bAddToSel && nRowGrid === nRowGridDragging) {

          let cellRH = null;
          try {
            cellRH = grid.cell("", nRowGrid);
          }
          catch(e) {
            let colG = null;
            const lstCols = grid.columns;
            for(let nC=1; nC<lstCols.length; ++nC) {
              colG = lstCols.byIndex(nC);
              cellRH = colG === null ? null : grid.cell(colG.name, nRowGrid);
              if(cellRH !== null)
                break;
            }
          }
          if(cellRH !== null) {
            const nRowTable : any = cellRH.tableRowIndex;
            if(nRowTable !== null)
            dframe.currentRow = nRowTable;
          }
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

            try {
              cellRH = grid.cell("", nRow);
            }
            catch(e) {
              let colG = null;
              const lstCols = grid.columns;
              for(let nC=1; nC<lstCols.length; ++nC) {
                colG = lstCols.byIndex(nC);
                cellRH = colG === null ? null : grid.cell(colG.name, nRowGrid);
                if(cellRH !== null)
                  break;
              }
            }

            if(cellRH !== null && cellRH.tableRowIndex !== null) {
              nRowTable = cellRH.tableRowIndex;
              bitsetSel.set(nRowTable, true, true);
            }
          }
        }

        const cell = grid.cell(colGrid.name, nRowGrid);
        const renderer = getRenderer(cell);
        if(renderer instanceof GridCellRendererEx) {
          renderer.onMouseUpEx(cell, eMouse, arXYMouseOnCellUp[0], arXYMouseOnCellUp[1]);
        }

        if(arXYMouseOnCellUp[0] === arXYMouseOnCellDown[0] && arXYMouseOnCellDown[1] === arXYMouseOnCellUp[1]) {
          if(renderer instanceof GridCellRendererEx) {
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

    let cellCurrent : DG.GridCell | null = null;
    this.m_handlerMouseLeave = rxjs.fromEvent(document, 'mouseleave').subscribe((e) => {

      if(cellCurrent !== null) {
        const renderer = getRenderer(cellCurrent);
        if (renderer instanceof GridCellRendererEx) {
          const eMouse = e as MouseEvent;
          renderer.onMouseLeaveEx(cellCurrent, eMouse, -1, -1);
        }
        cellCurrent = null;
      }
    });


    this.m_handlerMouseMove = rxjs.fromEvent(document, 'mousemove').subscribe((e) => {

      if(DG.toDart(grok.shell.v) !== DG.toDart(viewTable)) {
        return;
      }

      const bResizing = nResizeRowGridDragging >= 0;
      if (bResizing) {

        //console.log("Dragging : " + headerThis.m_strColName);
        const eMouse = e as MouseEvent;
        const nYDiff = eMouse.clientY - nYResizeDraggingAnchor;
        let nHRowGrid = nHResizeRowsBeforeDrag + nYDiff;

        if (nHRowGrid < PinnedColumn.MIN_ROW_HEIGHT)
          nHRowGrid = PinnedColumn.MIN_ROW_HEIGHT;
        else if (nHRowGrid > PinnedColumn.MAX_ROW_HEIGHT)
          nHRowGrid = PinnedColumn.MAX_ROW_HEIGHT;

        let g = eCanvasThis.getContext('2d');
        if(g === null)
          return;

        g.fillStyle = "white";
        const nHHeaderCols = GridUtils.getGridColumnHeaderHeight(grid);
        g.fillRect(0,nHHeaderCols, eCanvasThis.offsetWidth, eCanvasThis.offsetHeight);

        grid.setOptions({
          rowHeight: nHRowGrid //this won't trigger onRowsRezized event, which is a DG bug
        });

        notifyAllPinnedColsRowsResized(headerThis, nHRowGrid, true);
        notifyAllColsRowsResized(grid, nHRowGrid, true);

        let header = null;
        const ar = grid.dart.m_arRowHeaders;
        for(let n=0; n<ar.length; ++n) {
          header = ar[n];
          g = header.m_root.getContext('2d');
          header.paint(g, grid);
        }

        try {
          const colGrid0 = grid.columns.byIndex(0);
          if (colGrid0 !== null)
            colGrid0.visible = false;//temporary addressed the DG bug
        }
        catch(e) {
          //DG bug
        }
        return;
      }

      const arXYOnCell = [-1,-1];
      let nRowGrid = PinnedColumn.hitTestRows(eCanvasThis, grid, e as MouseEvent, false, arXYOnCell);
      if(nRowGrid >= 0) {
        const cell = grid.cell(colGrid.name, nRowGrid);
        const renderer = getRenderer(cell);

        if (renderer instanceof GridCellRendererEx) {

          if (cellCurrent === null) {
            renderer.onMouseEnterEx(cell, e as MouseEvent, arXYOnCell[0], arXYOnCell[1]);
          }

          if (cellCurrent !== null && nRowGrid !== cellCurrent.gridRow) {
             renderer.onMouseLeaveEx(cellCurrent, e as MouseEvent, -1, -1);

           renderer.onMouseEnterEx(cell, e as MouseEvent, arXYOnCell[0], arXYOnCell[1]);
          }

          renderer.onMouseMoveEx(cell, e as MouseEvent, arXYOnCell[0], arXYOnCell[1]);
         }

        cellCurrent = cell;
      }
      else if (cellCurrent !== null) {
        const renderer = getRenderer(cellCurrent);
        if (renderer instanceof GridCellRendererEx) {
          renderer.onMouseLeaveEx(cellCurrent, e as MouseEvent, -1, -1);
        }

        cellCurrent = null;
      }

      nRowGrid = PinnedColumn.hitTestRows(eCanvasThis, grid, e as MouseEvent, true, undefined);
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
  }

  isPinned() : boolean {
    return this.m_colGrid !== null;
  }

  getGridColumn() : DG.GridColumn | null{
    return this.m_colGrid;
  }

  getWidth() : number {
    return this.m_root === null ? -1 : this.m_root.offsetWidth;
  }

  getRoot() {
    return this.m_root;
  }

  public close() : void {

    if(this.m_colGrid === null) {
      throw new Error("Column has already been unpinned");
    }

    if(this.m_observerResizeGrid !== null) {
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

    const grid = this.m_colGrid.grid;
    const dart = DG.toDart(grid);
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
    try {
      this.m_colGrid.visible = true;
    }
    catch(e) {
      //DG bug
    }
    this.m_colGrid = null;


    if (dart.m_arRowHeaders.length === 0) {
      const colGrid0 = grid.columns.byIndex(0);
      if(colGrid0 !== null) {
        try{colGrid0.visible = true;}
        catch(e) {
          console.error("ERROR: Couldn't set visible property to true");
        }
      }
    }
    if(this.m_root === null)
      throw new Error('Root cannot be null');

    grid.canvas.style.left = (grid.canvas.offsetLeft - this.m_root.offsetWidth).toString() + "px";
    grid.overlay.style.left= (grid.overlay.offsetLeft - this.m_root.offsetWidth).toString() + "px";
    grid.canvas.style.width = (grid.canvas.offsetWidth + this.m_root.offsetWidth).toString() + "px";
    grid.overlay.style.width= (grid.overlay.offsetWidth + this.m_root.offsetWidth).toString() + "px";

    if(this.m_root.parentNode !== null)
     this.m_root.parentNode.removeChild(this.m_root);

    this.m_root = null;
  }


  private paint(g : CanvasRenderingContext2D | null, grid : DG.Grid) : void {
    //const nWDiv = entry.contentBoxSize ? entry.contentBoxSize[0].inlineSize : entry.contentRect.width;

    if(g === null) {
      return;
    }

    if(this.m_root === null) {
      throw new Error('Root cannot be null.');
    }

    if(this.m_colGrid === null) {
      throw new Error('Column grid cannot be null.');
    }
    const dframe = grid.dataFrame;
    const nW = this.m_root.offsetWidth;
    const nH = this.m_root.offsetHeight;

    g.fillStyle = "white";
    g.fillRect(0,0, nW, nH);

    if(this.m_colGrid.name === null)
      return;

    const bitsetFilter = dframe.filter;
    if(bitsetFilter.falseCount === dframe.rowCount)
      return;

    const options : any = grid.getOptions(true);
    g.font = options.colHeaderFont == null || options.colHeaderFont === undefined ? "bold 13px Roboto, Roboto Local" : options.colHeaderFont;
    let str = TextUtils.trimText(this.m_colGrid.name, g, nW);

    const tm = g.measureText(str);
    const nWLabel = tm.width;

    const nAscent = Math.abs(tm.actualBoundingBoxAscent);
    const nDescent = tm.actualBoundingBoxDescent;
    const nHFont =  nAscent + nDescent;// + 2*nYInset;

    let nX = 0;
    let nY = 0;
    const nHCH = GridUtils.getGridColumnHeaderHeight(grid);
    g.translate(0,0);
    g.textAlign = 'start';
    g.fillStyle = "Black";
    const nXX = nX + ((nW - nWLabel) >> 1);
    const nYY = nY + nHCH-2;
    //onsole.log("nXX " + nXX + " nYY = " + nYY + " CHH " + nHCH);
    g.fillText(str, nXX, nYY);


    const nRowCurrent =  dframe.currentRow.idx;
    const bitsetSel = dframe.selection;

    const nGridRowCount = dframe.filter.trueCount;
    const nGridfalseCount = dframe.filter.falseCount;

    const scrollV = grid.vertScroll;
    const nRowMin = Math.floor(scrollV.min);
    let nRowMax = Math.ceil(scrollV.max);

    let nHH = grid.root.offsetHeight - nHCH;//GridUtils.getColumnHeaderHeight(grid);//.style.height;
    const nHRow = GridUtils.getGridRowHeight(grid);
    const nRCount = Math.round(nHH/nHRow) +1;
    nRowMax = nRowMin + nRCount;

    if(nRowMax >= nGridRowCount)
      nRowMax = nGridRowCount -1;

    //console.log(nRowMin + " " + nRowMax);

    const nYOffset = nHCH;//GridUtils.getColumnHeaderHeight(grid);
    const nHRowGrid = nHRow;//GridUtils.getRowHeight(grid);
    let cellRH = null;


    let nRowTable = -1;
    nY = -1;
    let bSel = false;
    let bFiltered = false;
    for(let nRG=nRowMin; nRG<=nRowMax; ++nRG)
    {
      try{cellRH = grid.cell(this.m_colGrid.name, nRG);}
      catch(e)     //to address DG bug when everything is filtered
      {
        continue;
      }

      if(cellRH.tableRowIndex === undefined)//DG bug
        continue;

      nRowTable = cellRH.tableRowIndex === null ? -1 : cellRH.tableRowIndex;
      nY = nYOffset + (nRG - nRowMin)*nHRowGrid;

      let renderer : any = GridUtils.getGridColumnRenderer(cellRH.gridColumn);
      if(renderer === null) {
        renderer = cellRH.renderer;
      }


      if(nW > 0 && nHRowGrid > 0) { //to address a bug caused either DG or client app
        try {
         renderer.render(g, 0, nY, nW, nHRowGrid, cellRH, cellRH.style);
        } catch(e) {
           throw e;
        }
      }

      if(options.look.showRowGridlines) {
        g.strokeStyle = "Gainsboro";
        g.beginPath();
        g.moveTo(0, nY + nHRowGrid);
        g.lineTo(nW - 1, nY + nHRowGrid);
        g.stroke();
      }

      bSel = nRowTable < 0 ? false : bitsetSel.get(nRowTable);
      if(bSel)
      {
        g.globalAlpha = 0.2;
        g.fillStyle = PinnedColumn.SELECTION_COLOR;
        g.fillRect(0, nY, nW, nHRowGrid);
        g.globalAlpha = 1;
      }

      if(nRowCurrent === nRowTable)
      {
        g.globalAlpha = 0.2;
        g.fillStyle = PinnedColumn.ACTIVE_CELL_COLOR;
        g.fillRect(0, nY, nW, nHRowGrid);
        g.globalAlpha = 1;
      }
    }
  }


  private static hitTestRows(eCanvasPinned : HTMLCanvasElement, grid : DG.Grid, e : MouseEvent, bBorder : boolean, arXYOnCell : Array<number> | undefined)
  {
    const rect = eCanvasPinned.getBoundingClientRect();
    const scrollLeft= window.pageXOffset || document.documentElement.scrollLeft;
    const scrollTop = window.pageYOffset || document.documentElement.scrollTop;
    const nY = rect.top  + scrollTop;
    const nX = rect.left + scrollLeft;

    if(nX <= e.clientX && e.clientX <= nX + eCanvasPinned.offsetWidth)   //on the rows header
    {
      const nHHeaderCols = GridUtils.getGridColumnHeaderHeight(grid);
      const nHRowGrid = GridUtils.getGridRowHeight(grid);

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

        if(bBorder && Math.abs(nYDiff) <= PinnedColumn.Y_RESIZE_SENSITIVITY)
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
}
