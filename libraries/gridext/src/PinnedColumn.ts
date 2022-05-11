import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as rxjs from 'rxjs';

function trimText(str : string, g : any, nWidth : number) {
  let tm = g.measureText(str);
  let nW  = tm.width;
  if(nW <= nWidth)
    return str;

  let nHFont = tm.actualBoundingBoxAscent + tm.actualBoundingBoxDescent;
  let strDots = "...";
  tm = g.measureText(strDots);
  let nWDots  = tm.width;
  if(nWDots > nWidth)
  {
    strDots = "..";
    tm = g.measureText(strDots);
    nWDots = tm.width;
    if(nWDots <= nWidth)
      return  strDots;

    strDots = ".";
    tm = g.measureText(strDots);
    nWDots = tm.width;
    return nWDots <= nWidth ? strDots : "";
  }

  let nWAvail = nWidth - nWDots;
  let strW = "W";
  tm = g.measureText(strW);
  let nWW = tm.width;

  let nCharCount = Math.floor(nWAvail / nWW);
  let strAdj = str.substring(0, nCharCount);
  tm = g.measureText(strAdj);
  if(tm.width > nWAvail)
  {
    for(var n=nCharCount -1; n>=0; --n)
    {
      strAdj = str.substring(0, n);
      tm = g.measureText(strAdj);
      if(tm.width <= nWAvail)
        return strAdj + strDots;
    }
  }
  else
  {
    let strAdjOld = strAdj;
    for(var n=nCharCount +1; n<str.length; ++n)
    {
      strAdj = str.substring(0, n);
      tm = g.measureText(strAdj);
      if(tm.width > nWAvail)
        return strAdjOld + strDots;

      strAdjOld = strAdj;
    }
  }
}




export function closeAllPinnedColumns(grid : DG.Grid)
{
  const dart = DG.toDart(grid);
  let colPinned = null;
  const nPinnedColCount =dart.m_arRowHeaders === undefined ? 0 : dart.m_arRowHeaders.length;
  for(let n=0; n<nPinnedColCount; ++n) {
    colPinned = dart.m_arRowHeaders[0];
    colPinned.close();
  }
}

export function installPinnedColumns(grid : DG.Grid)
{
  closeAllPinnedColumns(grid);

  let colGrid = null;
  let settings = null;
  const lstCols = grid.columns;
  const arPinnedCols = new Array();

  for(let nCol=0;nCol<lstCols.length; ++nCol) {
    colGrid = lstCols.byIndex(nCol);
    if(colGrid === null)
      continue;

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
    new PinnedColumn(colGrid);
  }
}

export function isPinnedColumn(colGrid : DG.GridColumn)
{
  const grid = colGrid.grid;
  const dart = DG.toDart(grid);

  if(dart.m_arRowHeaders === undefined)
    return false;

  let colPinned = null;
  const nPinnedColCount = dart.m_arRowHeaders.length;
  for(let n=0; n<nPinnedColCount; ++n) {
    colPinned = dart.m_arRowHeaders[n];
    if(DG.toDart(colPinned.m_colGrid) === DG.toDart(colGrid))
      return true;
  }

  return false;
}

export class PinnedColumn {

  private static MIN_ROW_HEIGHT = 20;
  private static MAX_ROW_HEIGHT = 500;
  private static SELECTION_COLOR = DG.Color.toRgb(DG.Color.colSelection);// "rgba(237, 220, 88, 0.15)";
  private static ACTIVE_CELL_COLOR = DG.Color.toRgb(DG.Color.currentRow);// "rgba(153, 237, 82, 0.35)";
  private static Y_RESIZE_SENSITIVITY = 2;

  private m_colGrid : DG.GridColumn | null;
  private m_root : HTMLCanvasElement | null;
  private m_observerResize : ResizeObserver | null;
  private m_observerResizeGrid : ResizeObserver | null;
  private m_handlerVScroll : any;
  private m_handlerRowsFiltering : any;
  private m_handlerCurrRow : any;
  private m_handlerSel : any;
  private m_handlerRowResized : any;
  private m_handlerMouseDown : any;
  private m_handlerMouseUp : any;
  private m_handlerMouseMove : any;
  private m_handlerContextMenu : any;

  constructor(colGrid : DG.GridColumn) {
    if (colGrid.cellType === 'html') {
      throw new Error("HTML columns cannot be pinned.");
    }

    if(isPinnedColumn(colGrid))
      throw new Error("Column '" + colGrid.name + "' is already pinned.");

    const grid = colGrid.grid;
    const dart = DG.toDart(grid);

    if(dart.m_arRowHeaders === undefined)
      dart.m_arRowHeaders = [];

    dart.m_arRowHeaders.push(this);

    const viewTable = grid.view;
    const dframe = grid.dataFrame;

    const nW = colGrid.width;
    this.m_colGrid = colGrid;

    colGrid.visible = false;

    if(colGrid.settings === null || colGrid.settings === undefined)
      colGrid.settings = {};

    colGrid.settings.isPinned = true; //this will be saved with the layout
    colGrid.settings.idxPinned = dart.m_arRowHeaders.length -1;

    grid.canvas.style.left = (grid.canvas.offsetLeft + nW).toString() + "px";
    grid.overlay.style.left= (grid.overlay.offsetLeft + nW).toString() + "px";

    const eCanvasThis = ui.canvas(nW, 1500);

    if(grid.canvas.parentNode === null)
      throw new Error("Parent node for canvas cannot bre null.");

    grid.canvas.parentNode.insertBefore(eCanvasThis, grid.canvas);
    this.m_root = eCanvasThis;


    //OnResize Row header
    const headerThis = this;
    this.m_observerResize = new ResizeObserver(entries => {
      const g = eCanvasThis.getContext('2d');
      for (let entry of entries) {
        headerThis.paint(g, grid);
      }
    });

    this.m_observerResize.observe(eCanvasThis);
    const colGtrid0 = grid.columns.byIndex(0);
    if(colGtrid0 !== null) {
      try{ colGtrid0.visible = false; }
      catch(e) {
        console.error("ERROR: Couldn't set visible property to false");
      }
    }


    //OnResize Row header
    this.m_observerResizeGrid = new ResizeObserver(entries => {

      const bCurrent =  DG.toDart(grok.shell.v) === DG.toDart(viewTable);
      if(!bCurrent)
        return;

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

    this.m_handlerRowResized = grid.onRowsResized.subscribe((e : any) => {
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


    this.m_handlerMouseDown = rxjs.fromEvent(document, 'mousedown').subscribe((e : Event) => {

      if(DG.toDart(grok.shell.v) !== DG.toDart(viewTable))
        return;

      const eMouse = e as MouseEvent;

      if(eMouse.button !== 0)
        return;

      nResizeRowGridMoving = -1;

      const bAddToSel : boolean = eMouse.ctrlKey || eMouse.shiftKey;

      let nRowGrid = bAddToSel ? -1 : PinnedColumn.hitTestRows(eCanvasThis, grid, eMouse, true);
      if (nRowGrid >= 0) {
        const options : any = grid.getOptions(true);
        const nHRows = options.look.rowHeight;

        nResizeRowGridDragging = nRowGrid;
        nYResizeDraggingAnchor = eMouse.clientY;
        nHResizeRowsBeforeDrag = nHRows;
      }
      else
      {
        nRowGrid = PinnedColumn.hitTestRows(eCanvasThis, grid, eMouse, false);

        nRowGridDragging = nRowGrid;
        nYDraggingAnchor = eMouse.clientY;
      }

      e.preventDefault();
      e.stopPropagation();
    });


    this.m_handlerMouseUp = rxjs.fromEvent(document, 'mouseup').subscribe((e) => {

      if(DG.toDart(grok.shell.v) !== DG.toDart(viewTable))
        return;

      const eMouse = e as MouseEvent;
      if(eMouse.button !== 0)
        return;

      nResizeRowGridDragging = -1;
      nResizeRowGridDragging = -1;
      nYResizeDraggingAnchor = -1;

      nResizeRowGridMoving = -1;

      document.body.style.cursor = "auto";

      if(nRowGridDragging >= 0) {
        const bAddToSel = eMouse.ctrlKey || eMouse.shiftKey;

        const nRowGrid = PinnedColumn.hitTestRows(eCanvasThis, grid, eMouse, false);
        if(!bAddToSel && nRowGrid === nRowGridDragging) {

          const cellRH = grid.cell("", nRowGrid);
          const nRowTable = cellRH.tableRowIndex;

          if(nRowTable !== null)
           dframe.currentRowIdx = nRowTable;
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
            cellRH = grid.cell("", nRow);
            nRowTable = cellRH.tableRowIndex === null ? -1 : cellRH.tableRowIndex;

            if(nRowTable >= 0)
             bitsetSel.set(nRowTable, true, true);
          }
        }


        nRowGridDragging = -1;
        nYDraggingAnchor = -1;
      }
      //e.preventDefault();
      //e.stopPropagation();

    });


    this.m_handlerMouseMove = rxjs.fromEvent(document, 'mousemove').subscribe((e) => {

      if(DG.toDart(grok.shell.v) !== DG.toDart(viewTable))
        return;

      const bDragging = nResizeRowGridDragging >= 0;
      if (bDragging) {

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
        const options : any = grid.getOptions(true);
        const nHHeaderCols = options.look.colHeaderHeight;
        g.fillRect(0,nHHeaderCols, eCanvasThis.offsetWidth, eCanvasThis.offsetHeight);

        grid.setOptions({
          rowHeight: nHRowGrid
        });

        let header = null;
        const ar = grid.dart.m_arRowHeaders;
        for(let n=0; n<ar.length; ++n) {
          header = ar[n];
          g = header.m_root.getContext('2d');
          header.paint(g, grid);
        }

        const colGrid0 = grid.columns.byIndex(0);
        if(colGrid0 !== null)
          colGrid0.visible = false;//temporary addressed the DG bug
        return;
      }


      if (nResizeRowGridDragging >= 0) {
        document.body.style.cursor = "row-resize";
        return;
      }

      const nRowGrid = PinnedColumn.hitTestRows(eCanvasThis, grid, e as MouseEvent, true);
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
    this.m_handlerContextMenu = grok.events.onContextMenu.subscribe((args : any) => {
      //this.m_handlerContextMenu = rxjs.fromEvent(viewTable.root, 'contextmenu').subscribe((e) => {
      const e = args.causedBy;
      const elem = document.elementFromPoint(e.clientX, e.clientY);//e.offsetY);
      const b = elem === eCanvasThis;

      if(b) {
        let menu = args.args.menu;//DG.Menu.popup();

        menu = menu.item("Unpin Column", (str : string) => {
          thisRowHeader.close();
        });

        menu = menu.item("Unpin All Columns", (str : string) => {
          closeAllPinnedColumns(grid);
        });

        e.preventDefault();
        e.stopPropagation();
      }
    });
  }

  public close() {

    if(this.m_colGrid === null)
      throw new Error("Column has already been unpinned");

    if(this.m_observerResizeGrid !== null) {
      this.m_observerResizeGrid.disconnect();
      this.m_observerResizeGrid = null;
    }

    if(this.m_observerResize !== null) {
      this.m_observerResize.disconnect();
      this.m_observerResize = null;
    }

    this.m_handlerVScroll.unsubscribe();
    this.m_handlerVScroll = null;

    this.m_handlerRowResized.unsubscribe();
    this.m_handlerRowResized = null;

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

    this.m_handlerMouseMove.unsubscribe();
    this.m_handlerMouseMove = null;

    this.m_handlerContextMenu.unsubscribe();
    this.m_handlerContextMenu = null;

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
    this.m_colGrid.visible = true;
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

    if(this.m_root.parentNode !== null)
     this.m_root.parentNode.removeChild(this.m_root);

    this.m_root = null;
  }

  getWidth() {
    return this.m_root === null ? -1 : this.m_root.offsetWidth;
  }


  private paint(g : any, grid : DG.Grid)
  {
    //const nWDiv = entry.contentBoxSize ? entry.contentBoxSize[0].inlineSize : entry.contentRect.width;

    if(this.m_root === null)
      throw new Error('Root cannot be null.');

    if(this.m_colGrid === null)
      throw new Error('Column grid cannot be null.');

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


    let str = trimText(this.m_colGrid.name, g, nW);
    g.font = "bold 13px Roboto, Roboto Local";
    const tm = g.measureText(str);
    const nWLabel = tm.width;

    const nAscent = Math.abs(tm.actualBoundingBoxAscent);
    const nDescent = tm.actualBoundingBoxDescent;
    const nHFont =  nAscent + nDescent;// + 2*nYInset;

    let nX = 0;
    let nY = 0;
    const options : any = grid.getOptions(true);
    //const nHHeaderCols = options.look.colHeaderHeight;
    //const nHRowGrid = options.rowHeight;

    const nH = options.look.colHeaderHeight;//GridUtils.getColumnHeaderHeight(grid);
    g.fillStyle = "Black";
    g.fillText(str, nX + ((nW - nWLabel)>>1), nY + nH-2);


    const nRowCurrent =  dframe.currentRow.idx;
    const bitsetSel = dframe.selection;

    const nGridRowCount = dframe.filter.trueCount;
    const nGridfalseCount = dframe.filter.falseCount;

    const scrollV = grid.vertScroll;
    const nRowMin = Math.floor(scrollV.min);
    let nRowMax = Math.ceil(scrollV.max);

    let nHH = grid.root.offsetHeight - nH;//GridUtils.getColumnHeaderHeight(grid);//.style.height;
    const nHRow = options.look.rowHeight;//GridUtils.getRowHeight(grid);
    let nRCount = Math.round(nHH/nHRow) +1;
    nRowMax = nRowMin + nRCount;

    if(nRowMax >= nGridRowCount)
      nRowMax = nGridRowCount -1;

    //console.log(nRowMin + " " + nRowMax);

    const nYOffset = nH;//GridUtils.getColumnHeaderHeight(grid);
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

      try{cellRH.renderer.render(g, 0, nY, nW, nHRowGrid, cellRH, cellRH.style);}
      catch(e) {
        throw e;
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


  private static hitTestRows(eCanvasPinned : HTMLCanvasElement, grid : DG.Grid, e : MouseEvent, bBorder : boolean)
  {
    const rect = eCanvasPinned.getBoundingClientRect();
    const scrollLeft= window.pageXOffset || document.documentElement.scrollLeft;
    const scrollTop = window.pageYOffset || document.documentElement.scrollTop;
    const nY = rect.top  + scrollTop;
    const nX = rect.left + scrollLeft;

    if(nX <= e.clientX && e.clientX <= nX + eCanvasPinned.offsetWidth)   //on the rows header
    {
      const options : any = grid.getOptions(true);
      const nHHeaderCols = options.look.colHeaderHeight;
      const nHRowGrid = options.look.rowHeight;

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

        if(!bBorder && nYBorder - nHRowGrid <= nYMouseOnHeader && nYMouseOnHeader <= nYBorder)
          return nRow -1;
      }
    }

    return -1;
  }
}
