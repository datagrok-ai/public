import * as DG from 'datagrok-api/dg';
import * as GridUtils from '../utils/GridUtils';
import * as TextUtils from '../utils/TextUtils';
import {GridCellRendererEx} from "./GridCellRendererEx";

function isNullText(cell : DG.Cell) : boolean {
  const bNull : boolean = cell === null || cell === undefined || cell.value === null || cell.value === undefined;
  return bNull;
}

export class ClickableTextRenderer extends GridCellRendererEx {

  render(g : CanvasRenderingContext2D, nX : number, nY : number, nW : number, nH : number, cellGrid : DG.GridCell, style : DG.GridCellStyle) : void {
    const cell : DG.Cell = cellGrid.cell;
    let str = isNullText(cell) ? null : cell.value.toString();
    if (str === null) {
      return;
    }
    str = TextUtils.trimText(str, g, nW);
    const strFont : string  = style.font;
    if (strFont !== null && strFont !== undefined && strFont !== '') {
      g.font = strFont;
    }
    let tm = g.measureText(str);
    const nWLabel = Math.round(tm.width);
    const nYInset = 2;
    tm = g.measureText('W');
    const nAscent = Math.abs(tm.actualBoundingBoxAscent);
    const nDescent = tm.actualBoundingBoxDescent;
    const nHFont : number = nAscent + nDescent + 2 * nYInset;

    const nDeltaY : number = Math.floor((nH - nHFont) / 2);
    const nYY = nY + nDeltaY + nHFont;
    g.fillStyle = 'black';
    const strBaseOld = g.textBaseline;
    g.textBaseline = 'top';
    const nXX = nX + ((nW - nWLabel) / 2);

    g.textAlign = 'start';
    g.fillStyle = 'DodgerBlue';
    g.fillText(str, nXX, nYY - nHFont + nYInset);
    g.textBaseline = strBaseOld;
  }

  isClickable(cellGrid : DG.GridCell, nXOnCell : number, nYOnCell : number) : boolean {
    const cell : DG.Cell = cellGrid.cell;
    const str : string = isNullText(cell) ? null : cellGrid.cell.value.toString();
    if (str === null) {
      return false;
    }
    const eCanvas = cellGrid.grid.overlay;
    const g = eCanvas.getContext('2d');
    if(g === null) {
      return false;
    }

    const strFont : string  = cellGrid.style.font;
    if (strFont !== null && strFont !== undefined && strFont !== '') {
      g.font = strFont;
    }

    let tm = g.measureText(str);
    const nWLabel = Math.round(tm.width);
    const nYInset = 2;
    tm = g.measureText('W');
    const nAscent = Math.abs(tm.actualBoundingBoxAscent);
    const nDescent = tm.actualBoundingBoxDescent;
    const nHFont = nAscent + nDescent + 2 * nYInset;

    const nHRow = GridUtils.getGridRowHeight(cellGrid.grid);
    const nDeltaY = Math.floor((nHRow - nHFont) / 2);
    const nXX = 0 + ((cellGrid.gridColumn.width - nWLabel) / 2);
    const bX = nXX <= 0 || (nXOnCell >= nXX && nXOnCell <= nXX + nWLabel);
    if (!bX) {
      return false;
    }
    const bY = nDeltaY < 0 || (nYOnCell >= nDeltaY && nYOnCell <= nDeltaY + nHFont);
    return bY;
  }


  onMouseUpEx(cellGrid : DG.GridCell, e : MouseEvent, nXOnCell : number, nYOnCell : number) : void {
    if (document.body.style.cursor !== 'auto') {
      document.body.style.cursor = 'auto';
    }
  }

  onMouseMoveEx(cellGrid : DG.GridCell, e : MouseEvent, nXOnCell : number, nYOnCell : number) : void {
    //console.log('Move');

    if(e.buttons !== 0) {
      return;
    }

    const b = this.isClickable(cellGrid, nXOnCell, nYOnCell);
    if (b) {
      if (document.body.style.cursor !== 'pointer') {
        document.body.style.cursor = 'pointer';
      }
    } else {
      if (document.body.style.cursor !== 'auto') {
        document.body.style.cursor = 'auto';
      }
    }
    //this.onMouseMove(cellGrid, e);
  }


  onMouseLeaveEx(cellGrid : DG.GridCell, e : MouseEvent, nXOnCell : number, nYOnCell : number) : void {
    //console.log('Out');
    document.body.style.cursor = 'auto';
    //this.onMouseLeave(cellGrid, e);
  }
}
