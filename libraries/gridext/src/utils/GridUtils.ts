import * as DG from 'datagrok-api/dg';
import {GridCellRendererEx} from "../renderer/GridCellRendererEx";


export function getInstalledGridForColumn(colGrid : DG.GridColumn) : DG.Grid | null {
  const dart : any = DG.toDart(colGrid);
  if(!(dart.m_grid instanceof DG.Grid))
    return null;

  return dart.m_grid;
}

export function installGridForColumn(grid : DG.Grid, colGrid : DG.GridColumn) : boolean {
  if(colGrid.grid instanceof DG.Grid)
    return false;

  const dart : any = DG.toDart(colGrid);
  if(dart.m_grid instanceof DG.Grid)
    return false;

  dart.m_grid = grid;
  return true;
}


export function setGridColumnRenderer(colGrid : DG.GridColumn, renderer : GridCellRendererEx) : void {
  const dart : any = DG.toDart(colGrid);
  dart.m_renderer = renderer;
}

export function getGridColumnRenderer(colGrid : DG.GridColumn) : GridCellRendererEx | null {
  const dart : any = DG.toDart(colGrid);
  const renderer = dart.m_renderer;
  if(renderer === undefined)
    return null;

  return renderer;
}

export function getGridColumnHeaderHeight(grid : DG.Grid) : number {
  const options : any = grid.getOptions(true);
  let nHColHeader = options.look.colHeaderHeight;
  if(nHColHeader === null || nHColHeader === undefined) {//DG bug

    const cellGrid = grid.hitTest(2,2);//.cell(col.name, 0);
    if(cellGrid !== null) {
      const rc = cellGrid.bounds;
      return rc.y;
      //console.log('rc.y ' + rc.y + " rc.h= " + rc.height + " row " + cellGrid.gridRow + " name " +  cellGrid.gridColumn.name);
    }
  }
  return nHColHeader;
}

export function getGridRowHeight(grid : DG.Grid) : number {
  const options : any = grid.getOptions(true);
  const nHRow =  options.look.rowHeight;
  if(nHRow === null || nHRow === undefined) {//DG bug
    let col = null;
    const nColCount = grid.columns.length;
    for(let nCol=0; nCol<nColCount; ++nCol) {
      col = grid.columns.byIndex(nCol);
      if(col === null || !col.visible)
        continue;

      const cellGrid = grid.cell(col.name, 0);
      const rc = cellGrid.bounds;
      return rc.height;
    }
    return -1;
  }

  return nHRow;
}

export function getGridVisibleRowCount(grid : DG.Grid) : number {
  const dframe = grid.dataFrame;
  const bitsetFilter = dframe.filter;
  const nRowCount = bitsetFilter.trueCount;
  return nRowCount;
}

export function fillVisibleViewportRows(arMinMaxRowIdxs : Array<number>, grid : DG.Grid) : void {
  if(arMinMaxRowIdxs.length !== 2)
    throw new Error("Array to cobtain indices must have the length 2.");

  const scroll = grid.vertScroll;
  const nRowMin = Math.floor(scroll.min);
  let nRowMax = Math.ceil(scroll.max);
  const nRowCount = getGridVisibleRowCount(grid);
  if(nRowMax >= nRowCount) {
    nRowMax = nRowCount -1;
  }

  arMinMaxRowIdxs[0] = nRowMin;
  arMinMaxRowIdxs[1] = nRowMax;
  //console.log('min: ' + scroll.min + " max: " + scroll.max + ' nRowMax ' + nRowMax + " nVisRowCount: " + nRowCount);
}

export function fillVisibleViewportGridCells(arColRowIdxs : Array<number>, grid : DG.Grid)
{
  if(arColRowIdxs.length !== 4)
    throw new Error("Array to cobtain bound row column indices must have the length 4.");

  const arRowIdxs : Array<number> = [];
  const arColIdxs : Array<number> = [];
  const lstVisibleCells = grid.getVisibleCells();
  for(let cellGTmp of lstVisibleCells)
  {
    if(!arRowIdxs.includes(cellGTmp.gridRow))
      arRowIdxs.push(cellGTmp.gridRow);

    if(!arColIdxs.includes(cellGTmp.gridColumn.idx))
      arColIdxs.push(cellGTmp.gridColumn.idx);
  }

  const nRowMin = arRowIdxs.length === 0 ? -1 : arRowIdxs[0];
  const nRowMax = arRowIdxs.length === 0 ? -2 : arRowIdxs[arRowIdxs.length-1];

  arColRowIdxs[0] = arColIdxs.length === 0 ? -1 : arColIdxs[0];
  arColRowIdxs[1] = arColIdxs.length === 0 ? -2 : arColIdxs[arColIdxs.length -1];
  arColRowIdxs[2] = nRowMin;
  arColRowIdxs[3] = nRowMax;
}
