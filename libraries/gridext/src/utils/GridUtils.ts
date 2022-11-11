import * as DG from 'datagrok-api/dg';
import {GridCellRendererEx} from "../renderer/GridCellRendererEx";
import * as TextUtils from "./TextUtils";
/*
const canvas = ui.canvas(200*r, 300*r);

cellGrid.renderer.render(10, 10, 200, 300);

const r = window.devicePixelRatio;
x = r * x; y = r * y;
w = r * w; h = r * h;
*/

export function getGridDartPopupMenu() : HTMLElement | null {

  let eDiv = null;
  const cDiv = document.getElementsByClassName('d4-menu-item-container d4-vert-menu d4-menu-popup');
  for(let n=0; n<cDiv.length; ++n) {
    eDiv = (cDiv.item(n) as HTMLElement);
    return eDiv;
  }
  return null;
}


export function getToolIconDiv(grid : DG.Grid) : HTMLElement | null {

  const cImg = document.getElementsByClassName('grok-icon grok-font-icon-menu');

  let eDivHamb : HTMLElement | null = null;
  let eParent = null;
  for(let n=0; n<cImg.length; ++n) {
    eDivHamb = (cImg.item(n) as HTMLElement).parentElement;
    if(eDivHamb == null)
      return null;

    if(eDivHamb?.getAttribute('column_name') !== null) {//'data') === 'ColHamb') {
      eParent = eDivHamb.parentElement;
      while(eParent !== null) {
        if(eParent === grid.root)
          return eDivHamb;

        eParent = eParent.parentElement;
      }
    }
  }
  return null;
}



export function isHitTestOnElement(eElem: HTMLElement, eMouse : MouseEvent) : boolean {
  const eElemHit = document.elementFromPoint(eMouse.clientX, eMouse.clientY);
  const b = eElemHit == eElem;
  return b;
}

export function isRowHeader(colGrid : DG.GridColumn) : boolean {
  return colGrid.idx === 0 || colGrid.name === 'row header';
}

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
  const nRowCount = bitsetFilter.trueCount + Array.from(grid.pinnedRows).length;
  //my changes pinned rows const nRowCount = bitsetFilter.trueCount + Array.from(grid.pinnedRows).length;
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

export function getActiveGridRow(grid: DG.Grid) {
  const dframe = grid.dataFrame;
  const nRowTableActive = dframe.currentRow.idx;
  const nGridColCount = grid.columns.length;
  let colGrid = null;
  let cellGrid = null;
  for(let nCol=0; nCol<nGridColCount; ++nCol) {
    colGrid = grid.columns.byIndex(nCol);
    if(colGrid?.visible) {

      const nGridRowCount = dframe.filter.trueCount;
      for(let nR=0; nR<nGridRowCount; ++nR) {
        cellGrid = grid.cell(colGrid.name, nR);
        if(cellGrid.tableRowIndex === null || nRowTableActive === null)
          continue;

        if(cellGrid.tableRowIndex === nRowTableActive)
          return nR;
      }
      return -1;
    }
  }
  return -1;
}


const m_mapScaledFonts = new Map();

export function scaleFont(font : string, fFactor : number) : string {
  if(fFactor === 1.0) {
    return font;
  }

  const strKey = font + fFactor.toString();
  let fontScaled = m_mapScaledFonts.get(strKey);
  if(fontScaled !== undefined)
    return fontScaled;

  const nFontSize : number = TextUtils.getFontSize(font);
  fontScaled = TextUtils.setFontSize(font, Math.ceil(nFontSize * fFactor));
  m_mapScaledFonts.set(strKey, fontScaled);

  return fontScaled;
}

export function paintColHeaderCell(g : CanvasRenderingContext2D | null, nX : number, nY : number, nW: number, nH: number, colGrid : DG.GridColumn) {

  if(g === null)
    return;

  g.fillStyle = "white";
  g.fillRect(nX*window.devicePixelRatio, nY*window.devicePixelRatio, nW*window.devicePixelRatio, nH*window.devicePixelRatio);

  const grid = colGrid.grid;
  const options : any = grid.getOptions(true);

  const font = options.look.colHeaderFont == null || options.look.colHeaderFont === undefined ? "bold 14px Volta Text, Arial" : options.look.colHeaderFont;
  const fontNew = scaleFont(font, window.devicePixelRatio);
  g.font = fontNew;

  let str = TextUtils.trimText(colGrid.name, g, nW);

  const tm = g.measureText(str);
  const nWLabel = tm.width;

  const nAscent = Math.abs(tm.actualBoundingBoxAscent);
  const nDescent = tm.actualBoundingBoxDescent;
  const nHFont =  nAscent + nDescent;// + 2*nYInset;

  const nHH = nH*window.devicePixelRatio;

  g.textAlign = 'start';
  g.fillStyle = "black";
  const nXX = nX*window.devicePixelRatio + Math.ceil(3*window.devicePixelRatio);//((nW*window.devicePixelRatio - nWLabel) >> 1);
  const nYY = (nY*window.devicePixelRatio + nHH - Math.ceil(3*window.devicePixelRatio));//-2*window.devicePixelRatio);
  g.fillText(str, nXX, nYY);
}


export const LeftArrow = '←';
export const RightArrow = '→';
export const UpArrow = '↑';
export const DownArrow = '↓';
