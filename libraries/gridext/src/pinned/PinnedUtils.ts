import * as DG from 'datagrok-api/dg';
import * as GridUtils from '../utils/GridUtils';
import {PinnedColumn} from "./PinnedColumn";
import {GridCellRendererEx} from "../renderer/GridCellRendererEx";

export function findPinnedColumnByRoot(eCanvas : HTMLCanvasElement, grid : DG.Grid) : PinnedColumn | null {
  const dart = DG.toDart(grid);
  let colPinned = null;
  const nPinnedColCount =dart.m_arRowHeaders === undefined ? 0 : dart.m_arRowHeaders.length;
  for(let n=0; n<nPinnedColCount; ++n) {
    colPinned = dart.m_arRowHeaders[n];
    if(colPinned.getRoot() === eCanvas) {
      return colPinned;
    }
  }
  return null;
}

export function addPinnedColumn(colGrid : DG.GridColumn) : PinnedColumn {
  const colPinned = new PinnedColumn(colGrid);
  return colPinned;
}

export function closeAllPinnedColumns(grid : DG.Grid) : void {
  const dart = DG.toDart(grid);
  let colPinned = null;
  const nPinnedColCount =dart.m_arRowHeaders === undefined ? 0 : dart.m_arRowHeaders.length;
  for(let n=0; n<nPinnedColCount; ++n) {
    colPinned = dart.m_arRowHeaders[0]; //0 is not a bug
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
    if(colGrid === null)
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
    new PinnedColumn(colGrid);
  }
}

export function isPinnedColumn(colGrid : DG.GridColumn) : boolean {
  const grid = colGrid.grid;
  const dart = DG.toDart(grid);

  if(dart.m_arRowHeaders === undefined)
    return false;

  let colPinned = null;
  const nPinnedColCount = dart.m_arRowHeaders.length;
  for(let nColPin=0; nColPin<nPinnedColCount; ++nColPin) {
    colPinned = dart.m_arRowHeaders[nColPin];
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

  let grid : DG.Grid | null = colGrid.grid;
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
  const cell = grid.hitTest(e.offsetX, e.offsetY);
  if (cell === undefined || cell === null || cell.cellType === null) //bug in DG , top left cell
    return;

  const colGrid = cell.gridColumn;
  if(!isPinnableColumn(colGrid)) {
    return;
  }

  if ((cell.isTableCell || cell.isColHeader) && (elem === grid.canvas || elem === grid.overlay)) {
    const menu = args.args.menu;
    fnMenuCallback(menu, colGrid, grid);
    e.preventDefault();
    e.stopPropagation();
    return;
  }
}
