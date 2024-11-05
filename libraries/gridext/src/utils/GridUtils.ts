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

function parseStyleValue(str: string) : number {
  const idx = str.indexOf('px');
  if (idx < 0)
    return NaN;

  const strVal = str.substring(0, idx);
  return parseInt(strVal);
}

export function getLeft(canvas: HTMLCanvasElement) : number {
  const str = canvas.style.left;
  const left = parseStyleValue(str);
  return left;
}

export function getWidth(canvas: HTMLCanvasElement) : number {
  const str = canvas.style.width;
  const width = parseStyleValue(str);
  return width;
}

export function getGridDartPopupMenu() : HTMLElement | null {
  let eDiv = null;
  const cDiv = document.getElementsByClassName('d4-menu-item-container d4-vert-menu d4-menu-popup');
  for (let n=0; n<cDiv.length; ++n) {
    eDiv = (cDiv.item(n) as HTMLElement);
    return eDiv;
  }
  return null;
}


export function getToolIconDiv(grid : DG.Grid) : HTMLElement | null {
  const cImg = document.getElementsByClassName('grok-icon grok-font-icon-menu');
  let eDivHamb : HTMLElement | null = null;
  let eParent = null;
  for (let n=0; n<cImg.length; ++n) {
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
  if (dart === null)
    return;

  dart.m_renderer = renderer;
}

export function getGridColumnRenderer(colGrid : DG.GridColumn) : GridCellRendererEx | null {
  const dart : any = DG.toDart(colGrid);
  if (dart === null)
    return null;

  const renderer = dart.m_renderer;
  if (renderer === undefined)
    return null;

  return renderer;
}

export function getGridColumnHeaderHeight(grid : DG.Grid, colHeaderHeight: number = -1) : number {
  const nHColHeader = colHeaderHeight !== -1 ? colHeaderHeight : grid.getOptions(true).look.colHeaderHeight;
  if(nHColHeader === null || nHColHeader === undefined) {//DG bug

    const cellGrid = grid.hitTest(2,2);//.cell(col.name, 0);
    if(cellGrid !== null && cellGrid !== undefined) {
      const rc = cellGrid.bounds;
      return rc.y;
      //console.log('rc.y ' + rc.y + " rc.h= " + rc.height + " row " + cellGrid.gridRow + " name " +  cellGrid.gridColumn.name);
    }
    else return 30;
  }
  return nHColHeader;
}

export function getGridRowHeight(grid : DG.Grid, rowHeight: number = -1) : number {
  const nHRow = rowHeight !== -1 ? rowHeight : grid.getOptions(true).look.rowHeight;
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


//Calculates individual words layout from an arbitrary string within a specified bounded range.
//The text is split only by whitespace as other characters are considered to a part of individual words.
//The priority is given to keep the words as a whole. Truncation occurs only used when an individual word cannot fit withing the range's bounds
//ctx the canvas' graphics context. Should be initialized with font and alignment flags before calling this function.
//arLines an array that after the call will contain the words layout. Each element of the array is another array representing the layout of an individual line.
//strText a string containing text to be layouted.
//nWidth the width of bounded range.
export function calcWordsCellLayout(ctx: CanvasRenderingContext2D, arLines: string[][], strText: string, nWidth: number) {
  while(arLines.length > 0) {
    arLines.pop();
  }

  if(strText === "")
    return;

  arLines.push([]);

  const arWords : Array<string> = [];
  splitWithWhitespaces(strText, arWords);

  let tm = null;
  let strWord = "";
  let nLine = 0;
  let nWLine = 0;
  let nWWord = -1;
  let nCharCount = -1;
  for (let nWord = 0; nWord < arWords.length; ++nWord) {
    strWord = arWords[nWord];
    tm = ctx.measureText(strWord);
    nWWord = tm.width;
    if (nWLine + nWWord > nWidth) { //the word cannot fit within the cell's bounds
      if (arLines[nLine].length === 0) { // the word gets truncated if the line starts with it. Otherwise, it goes to the next line (eventually it will be truncated)
        nCharCount = fitWordPart(ctx, strWord, nWidth);
        if (nCharCount === 0)
          return;

        strWord = strWord.substring(0, nCharCount);
        arLines[nLine].push(strWord);   //replace the word with the remaining part to repeat on the next iteration
        arWords[nWord] = arWords[nWord].substring(nCharCount);
      }

      --nWord;
      nWLine = 0;
      arLines.push([]);
      ++nLine;
      continue;
    }

    arLines[nLine].push(strWord);
    nWLine += nWWord;
  }
}

//Split an arbitrary text into words by whitespace character. Whitespaces as considered as individual words, and thus are included into the output.
//strText a string containing some text to be split.
//ar the specified array that after the call will contain individual words including whitespaces as legitimate words.
function splitWithWhitespaces(strText: string, ar: Array<string>) {

  while (ar.length > 0) {
    ar.pop();
  }

  let nIdxStart=0;
  let nIdxSpace=-1;

  while (nIdxStart < strText.length) {
    nIdxSpace = strText.indexOf(" ", nIdxStart);
    if (nIdxSpace >= 0) {
      if(nIdxSpace> nIdxStart)
        ar.push(strText.substring(nIdxStart, nIdxSpace))

      ar.push(" ");
      nIdxStart = nIdxSpace + 1;
    } else {
      ar.push(strText.substring(nIdxStart, strText.length));
      break;
    }
  }
}


//Calculates the number of characters from a word that fits within a specified bounded range.
//strWord the specified word to be split
//nWidth the width of the specified bounded range.
function fitWordPart(g: CanvasRenderingContext2D, strWord: string, nWidth: number) {
  let str = "";
  let nWWWord = -1;
  let nLength = strWord.length;
  for (var n = nLength; n >= 0; --n) {
    str = strWord.substring(0, n);
    nWWWord = g.measureText(str).width;
    if(nWWWord <= nWidth)
      return n;
  }

  return 0;
}

export const LeftArrow = '←';
export const RightArrow = '→';
export const UpArrow = '↑';
export const DownArrow = '↓';
