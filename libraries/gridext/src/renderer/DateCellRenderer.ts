import {GridCellRendererEx} from "./GridCellRendererEx";
import * as MathUtils from "../utils/MathUtils";
import * as RendUtils from './RendUtils'
import * as DG from 'datagrok-api/dg';
import * as TextUtils from "../utils/TextUtils";

const MONTHS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];

function isNullText(cell : DG.Cell) : boolean {
  const bNull : boolean = cell === null || cell === undefined || cell.value === null || cell.value === undefined;
  return bNull;
}

function formatDate(nTime : number) : string {
  nTime = Math.floor(nTime);

  const date = new Date(nTime);
  const nMonth = date.getMonth() +1;
  const nDay = date.getDate();
  const nYear = date.getFullYear();

  const str = nMonth.toString() + "-" + nDay + "-" + nYear;
  return str;
}

function formatDate2Day(nTime : number) : string {
  nTime = Math.floor(nTime);

  const date = new Date(nTime);
  const nDay = date.getDate();

  const str = nDay.toString();
  return str;
}

function formatDate2Month(nTime : number) : string {
  nTime = Math.floor(nTime);

  const date = new Date(nTime);
  const nMonth = date.getMonth();
  const str = MONTHS[nMonth];
  return str;
}

function formatDate2DayMonth(nTime : number) : string {
  nTime = Math.floor(nTime);

  const date = new Date(nTime);
  const nMonth = date.getMonth() +1;
  const nDay = date.getDate();

  const str = MONTHS[nMonth] + " " + nDay;
  return str;
}

function formatDate2Year(nTime : number) : string {
  nTime = Math.floor(nTime);

  const date = new Date(nTime);
  const nYear = date.getFullYear();

  const str = nYear.toString();
  return str;
}


export class DateCellRenderer extends GridCellRendererEx {
  render(g: CanvasRenderingContext2D, nX: number, nY: number, nW: number, nH: number, cellGrid: DG.GridCell, style: DG.GridCellStyle): void {
    const cell: DG.Cell = cellGrid.cell;
    if(cell.column.type !== DG.COLUMN_TYPE.DATE_TIME) {
      return;
    }

    let nTime = NaN;
    const ob = cell.value;
    if(MathUtils.isNullValue(ob)) {
      nTime = NaN;
    }
    else if (typeof ob !== "number")
      nTime = ob.a;

    if(isNaN(nTime)) {
     return;
    }

    const str = formatDate(nTime);

    const strFont = style.font;
    if (strFont !== null && strFont !== undefined && strFont !== '') {
      g.font = strFont;
    }

    let tm = g.measureText(str);
    const nWLabel = Math.round(tm.width);
    if(nWLabel < nW) {
      RendUtils.renderXYCenteredText(str, g, nX, nY, nW, nH, style.font, 'black');
      return;
    }

    const nYInset = 2;
    tm = g.measureText('W');
    const nAscent = Math.abs(tm.actualBoundingBoxAscent);
    const nDescent = tm.actualBoundingBoxDescent;
    const nHFont : number = nAscent + nDescent;
    const bTwoLevel = nH - 2*nHFont - 3*nYInset >= 0;

    if(bTwoLevel) {
        let strMonthDay = formatDate2DayMonth(nTime);
        let strYear = formatDate2Year(nTime);
        tm = g.measureText(strMonthDay);
        const nWMonthDay = Math.round(tm.width);
        tm = g.measureText(strYear);
        const nWYear = Math.round(tm.width);

        if(nWMonthDay <= nW && nWYear <= nW) {
          const nYOffset = Math.floor((nH - 2*nHFont - 3*nYInset)/2);
          RendUtils.renderXCenteredText(strMonthDay, g, nX, nY + nYOffset, nW, nH, style.font, 'black');
          RendUtils.renderXCenteredText(strYear, g, nX, nY + 2*nYOffset + nHFont, nW, nH, style.font, 'black');
          return;
       }
       else if(nWMonthDay > nW) {
          const bTriLevel = nH - 3*nHFont - 4*nYInset >= 0;
          if (bTriLevel) {
            let strMonth = formatDate2Month(nTime);
            let strDay = formatDate2Day(nTime);
            strDay = TextUtils.trimText(strDay, g, nW);
            strMonth = TextUtils.trimText(strMonth, g, nW);
            strYear = TextUtils.trimText(strYear, g, nW);
            const nYOffset = Math.floor((nH - 3 * nHFont - 4 * nYInset) / 2);
            RendUtils.renderXCenteredText(strMonth, g, nX, nY + nYOffset, nW, nH, style.font, 'black');
            RendUtils.renderXCenteredText(strDay, g, nX, nY + 2 * nYOffset + nHFont, nW, nH, style.font, 'black');
            RendUtils.renderXCenteredText(strYear, g, nX, nY + 3 * nYOffset + 2 * nHFont, nW, nH, style.font, 'black');
            return;
          }
        }
       strMonthDay = TextUtils.trimText(strMonthDay, g, nW);
       strYear = TextUtils.trimText(strYear, g, nW);
       const nYOffset = Math.floor((nH - 2*nHFont - 3*nYInset)/2);
       RendUtils.renderXCenteredText(strMonthDay, g, nX, nY + nYOffset, nW, nH, style.font, 'black');
       RendUtils.renderXCenteredText(strYear, g, nX, nY + 2*nYOffset + nHFont, nW, nH, style.font, 'black');
       return;
    }

    RendUtils.renderXYCenteredText(str, g, nX, nY, nW, nH, style.font, 'black');
  }
}
