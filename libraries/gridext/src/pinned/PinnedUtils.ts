import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as GridUtils from '../utils/GridUtils';
import * as PinnedUtils from "./PinnedUtils";
import {PinnedColumn} from "./PinnedColumn";
import {GridCellRendererEx} from "../renderer/GridCellRendererEx";
import {TableView} from "datagrok-api/dg";

function getGrid(colGrid : DG.GridColumn) : DG.Grid | null {
  let grid : DG.Grid | null = colGrid.grid;
  if( grid === null) {
    grid = GridUtils.getInstalledGridForColumn(colGrid);
    if(grid instanceof DG.Grid)
      return grid;
  }

  return grid;
}


export function getTotalPinnedRowsHeight(grid : DG.Grid) : number {

  const dart = DG.toDart(grid);

  const nPinnedRowsCounnt = dart.m_arPinnedRows === undefined ? 0 : dart.m_arPinnedRows.length;
  let rowPinned = null;
  let nHeight = 0;
  for(let nR=0; nR<nPinnedRowsCounnt; ++nR) {
    rowPinned = dart.m_arPinnedRows[nR];
    nHeight += rowPinned.getHeight();
  }

  return nHeight;
}


export function getPinnedColumnLeft(colPinned : PinnedColumn) : number {
  const grid = colPinned.getGridColumn()?.grid;
  const dart = DG.toDart(grid);
  let colPinnedTmp = null;

  let nW = 0;
  const nPinnedColCount =dart.m_arPinnedCols === undefined ? 0 : dart.m_arPinnedCols.length;
  for(let n=0; n<nPinnedColCount; ++n) {
    colPinnedTmp = dart.m_arPinnedCols[n];
    if(colPinned === colPinnedTmp) {
      break;
    }
    nW += colPinnedTmp.getWidth();
  }

  return nW;
}


export function setPinnedColumnWidth(colPinned : PinnedColumn, nW : number) : void {
  const grid = colPinned.getGridColumn()?.grid;
  if(grid === null || grid === undefined)
    return;

  const dart = DG.toDart(grid);
  let colPinnedTmp = null;
  const nPinnedColCount =dart.m_arPinnedCols === undefined ? 0 : dart.m_arPinnedCols.length;
  let n=0;
  for(; n<nPinnedColCount; ++n) {
    colPinnedTmp = dart.m_arPinnedCols[n];
    if(colPinned === colPinnedTmp) {
      break;
    }
  }

  if(n === nPinnedColCount)
    return;

  let eCanvas = colPinnedTmp.getRoot();
  if(eCanvas === null)
    return;

  //This Pinned Column
  const nWCol = colPinned.getWidth();
  let nXDiff = nW - nWCol;
  if(grid.canvas.offsetWidth - nXDiff < 20) {
    nXDiff = grid.canvas.offsetWidth - 20;
    nW = nWCol + nXDiff;
  }
  eCanvas.width = nW*window.devicePixelRatio;
  eCanvas.style.width = nW.toString() + "px";
  //console.log('nWWWCol: ' + nW + " nXDiff= " + nXDiff + " new W= " + nXDiff);

  const g = eCanvas.getContext('2d');
  colPinned.paint(g, grid);

  ++n;
  for(; n<nPinnedColCount; ++n) {
    colPinnedTmp = dart.m_arPinnedCols[n];
    eCanvas = colPinnedTmp.getRoot();
    if(eCanvas === null)
      continue;

    eCanvas.style.left = (eCanvas.offsetLeft + nXDiff).toString() + "px";
   }

  //Grid
  grid.canvas.style.left = (grid.canvas.offsetLeft + nXDiff).toString() + "px";
  grid.overlay.style.left= (grid.overlay.offsetLeft + nXDiff).toString() + "px";

  grid.canvas.style.width = (grid.canvas.offsetWidth - nXDiff).toString() + "px";
  grid.overlay.style.width= (grid.overlay.offsetWidth - nXDiff).toString() + "px";
  grid.invalidate();
}


export function getTotalPinnedColsWidth(grid : DG.Grid) : number {
  const dart = DG.toDart(grid);
  const nPinnedColsCounnt = dart.m_arPinnedCols === undefined ? 0 : dart.m_arPinnedCols.length;
  let colPinned = null;
  let nWidth = 0;
  for(let nR=0; nR<nPinnedColsCounnt; ++nR) {
    colPinned = dart.m_arPinnedCols[nR];
    nWidth += colPinned.getWidth();
  }

  return nWidth;
}


export function findPinnedColumnByRoot(eCanvas : HTMLCanvasElement, grid : DG.Grid) : PinnedColumn | null {
  const dart = DG.toDart(grid);
  let colPinned = null;
  const nPinnedColCount =dart.m_arPinnedCols === undefined ? 0 : dart.m_arPinnedCols.length;
  for(let n=0; n<nPinnedColCount; ++n) {
    colPinned = dart.m_arPinnedCols[n];
    if(colPinned.getRoot() === eCanvas) {
      return colPinned;
    }
  }
  return null;
}

export function getPinnedColumnCount(grid : DG.Grid) : number {
  const dart = DG.toDart(grid);
  const nPinnedColCount =dart.m_arPinnedCols === undefined ? 0 : dart.m_arPinnedCols.length;
  return nPinnedColCount;
}


export function getPinnedColumn(nIdx : number, grid : DG.Grid) : PinnedColumn {
  if(nIdx < 0) {
    throw new Error("Pinned column index cannot be negative: " + nIdx);
  }

  const dart = DG.toDart(grid);
  const nPinnedColCount =dart.m_arPinnedCols === undefined ? 0 : dart.m_arPinnedCols.length;
  if(nIdx >= nPinnedColCount) {
    throw new Error("Pinned column index cis out of bounds [0,: " + (nPinnedColCount -1) + "]");
  }

  return dart.m_arPinnedCols[nIdx];
}

export function addPinnedColumn(colGrid : DG.GridColumn) : PinnedColumn {
  const colPinned = new PinnedColumn(colGrid);
  return colPinned;
}

export function closeAllPinnedColumns(grid : DG.Grid) : void {
  const dart = DG.toDart(grid);
  let colPinned = null;
  const nPinnedColCount =dart.m_arPinnedCols === undefined ? 0 : dart.m_arPinnedCols.length;
  for(let n=1; n<nPinnedColCount; ++n) {
    colPinned = dart.m_arPinnedCols[1]; //0 is not a bug
    colPinned.close();
  }
}

export function installPinnedColumns(grid : DG.Grid) : void {
  closeAllPinnedColumns(grid);

  let colGrid = null;
  let settings : any = null;
  const lstCols = grid.columns;
  const arColsToPin = new Array<DG.GridColumn>();

  for(let nCol=0;nCol<lstCols.length; ++nCol) {
    colGrid = lstCols.byIndex(nCol);
    if(colGrid === null || GridUtils.isRowHeader(colGrid))
      continue;

    settings = colGrid.settings;
    if(settings !== null && settings !== undefined && settings.isPinned) {
      arColsToPin.push(colGrid);
    }
  }

  arColsToPin.sort((colOne, colTwo) => {
    if(colOne.settings.idxPinned === colTwo.settings.idxPinned) {
      return 0;//throw new Error("Pinned indices cannot be equal for different columns");
    }

    return colOne.settings.idxPinned < colTwo.settings.idxPinned ? -1 : 1;
  });

  let bPinned = false;
  for(let n=0;n<arColsToPin.length; ++n) {
    colGrid = arColsToPin[n];
    if(isPinnableColumn(colGrid)) {
      new PinnedColumn(colGrid);
      bPinned = true;
    }
  }


  const colGrid0 = !bPinned ? null : grid.columns.byIndex(0);
  if(colGrid0 !== null && colGrid0 !== undefined) {//DG Bug from reading layout
    try{
      setTimeout(() => {colGrid0.visible = false}, 100);
    }
    catch(e) {
      console.error("ERROR: Couldn't hide row header.");
    }
  }

}




export function isPinnedColumn(colGrid : DG.GridColumn) : boolean {
  const grid = getGrid(colGrid);
  const dart = DG.toDart(grid);

  if(dart.m_arPinnedCols === undefined)
    return false;

  let colPinned = null;
  const nPinnedColCount = dart.m_arPinnedCols.length;
  for(let nColPin=0; nColPin<nPinnedColCount; ++nColPin) {
    colPinned = dart.m_arPinnedCols[nColPin];
    if(DG.toDart(colPinned.m_colGrid) === DG.toDart(colGrid))
      return true;
  }

  return false;
}

export function isPinnableColumn(colGrid : DG.GridColumn) : boolean {
  const b = isPinnedColumn(colGrid);
  if(b) {
    return false;
  }


  let grid : DG.Grid | null = getGrid(colGrid);
  if(!(grid instanceof DG.Grid)) {
    grid = GridUtils.getInstalledGridForColumn(colGrid);

    if(!(grid instanceof DG.Grid)) {
      return false;
    }
  }

  //temporary disable to allow de-serialization from layout if(grid.canvas.offsetWidth <= colGrid.width) {
  //return false;
  //}

  if(colGrid.cellType === "html") {
    const renderer = GridUtils.getGridColumnRenderer(colGrid);
    if(!(renderer instanceof GridCellRendererEx)) {
      return false;
    }
  }

  return true;
}


let PINNED_COLUMNS_REGISTERED = false;

export function registerPinnedColumns() {
  if(PINNED_COLUMNS_REGISTERED)
    return;

  grok.events.onContextMenu.subscribe((args : any) => {
    PinnedUtils.handleContextMenu(args, (menu : DG.Menu, colGridOrPinned : DG.GridColumn | PinnedColumn,
                                         grid : DG.Grid) => {
      if (colGridOrPinned instanceof PinnedColumn) {
        const colGrid = colGridOrPinned.getGridColumn();
        if (colGrid !== null && !GridUtils.isRowHeader(colGrid)) {
          menu.item('Unpin Column', () => {
            colGridOrPinned.close();
          });
        }
        menu.item('Unpin All Columns', () => {
          PinnedUtils.closeAllPinnedColumns(grid);
        });
      } else {
        menu.item('Pin Column', async () => {
          PinnedUtils.addPinnedColumn(colGridOrPinned);
        });
      }
    });
  });

  grok.events.onViewLayoutApplied.subscribe((layout : DG.ViewLayout) => {
    const view : DG.TableView = layout.view as TableView;
    if(view === null) {
      console.error("View cannot be null; layout.view = null; grok.events.onViewLayoutApplied");
      return;
    }

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


  PINNED_COLUMNS_REGISTERED = true;
}


export function handleContextMenu(args : any, fnMenuCallback : Function) : void {

  const grid = args.args.context;
  if (!(grid instanceof DG.Grid)) {
    return;
  }
  const e = args.causedBy;
  //Check if we are on a pinned column
  const elem = document.elementFromPoint(e.clientX, e.clientY);
  if (elem instanceof HTMLCanvasElement) {
    const colPinned = findPinnedColumnByRoot(elem, grid);
    if (colPinned !== null) {
      let menu = args.args.menu;
      fnMenuCallback(menu, colPinned, grid);
      e.preventDefault();
      e.stopPropagation();
      return;
    }
  }
  //Regular Columns
  const cell = grid.hitTest(e.offsetX, e.offsetY);
  if (cell === undefined || cell === null || cell.cellType === null) //bug in DG , top left cell
    return;

  const colGrid = cell.gridColumn;
  if(!isPinnableColumn(colGrid)) {
    return;
  }

  if ((cell.isTableCell || cell.isColHeader)  && (elem === grid.canvas || elem === grid.overlay)) {
    const menu = args.args.menu;
    fnMenuCallback(menu, colGrid, grid);
    e.preventDefault();
    e.stopPropagation();
    return;
  }
}
