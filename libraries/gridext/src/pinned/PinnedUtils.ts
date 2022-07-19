import * as DG from 'datagrok-api/dg';
import * as GridUtils from '../utils/GridUtils';
import {PinnedColumn} from "./PinnedColumn";
import {GridCellRendererEx} from "../renderer/GridCellRendererEx";


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

  for(let n=0;n<arColsToPin.length; ++n) {
    colGrid = arColsToPin[n];
    if(isPinnableColumn(colGrid)) {
      new PinnedColumn(colGrid);
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

  if(grid.canvas.offsetWidth <= colGrid.width) {
    return false;
  }

  if(colGrid.cellType === "html") {
    const renderer = GridUtils.getGridColumnRenderer(colGrid);
    if(!(renderer instanceof GridCellRendererEx)) {
      return false;
    }
  }

  return true;
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
