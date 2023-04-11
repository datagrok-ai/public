import * as DG from 'datagrok-api/dg';
import {GridCellRendererEx} from "./GridCellRendererEx";
import * as TextUtils from "../utils/TextUtils";

export class AbstractVertLayoutTextRenderer extends GridCellRendererEx {

  private m_bGapToBottom: boolean = true;

  isGapToBottom() : boolean {
    return this.m_bGapToBottom;
  }

  setGapToBottom(bGap: boolean): void {
    this.m_bGapToBottom = bGap;
  }

  fillLabelsAndUI(cellGrid : DG.GridCell, arTextLabels : string[], arTextFonts: string[], arTextColors : string[], arBackColors: string[]) : void {
    throw new Error('Not Implemented');
  }

  render(g: CanvasRenderingContext2D, nX: number, nY: number, nW: number, nH: number, cellGrid: DG.GridCell, style: DG.GridCellStyle): void {
    super.render(g, nX, nY, nW, nH, cellGrid, style);

    const cell : DG.Cell = cellGrid.cell;
    const crBack = DG.Color.toHtml(style.backColor);
    if (crBack !== undefined) {
      g.fillStyle = crBack;
      g.fillRect(nX,nY, nW, nH);
    }

    const arTextLabels : string[] = [];
    const arTextFonts : string[] = [];
    const arTextColors : string[] = [];
    const arBackColors : string[] = [];
    this.fillLabelsAndUI(cellGrid, arTextLabels, arTextFonts, arTextColors, arBackColors);
    this.paintLabels(g, arTextLabels, arTextFonts, arTextColors, arBackColors, nX, nY, nW, nH);
  }

  paintLabels(g: CanvasRenderingContext2D, arTextLabels : string[], arTextFonts: string[], arTextColors: string[], arBackColors: string[], nX: number, nY: number, nW: number, nH: number) : void {
    //g.font = this.getFont();//"12px Arial";
    g.fillStyle = "white";
    g.strokeStyle = "black";
    g.textAlign = "left";
    g.textBaseline = "top";

    g.fillRect(nX, nY, nW, nH);

    const nYInset = 2;

    let nHAvail = nH;
    let nWAvail = nW;

    let nLabelCount = arTextLabels.length;

    let nHFont = -1;
    let nHSum = 0;
    let nHSumExcl0 = 0;
    let nHTmp = 0;
    let nFittedRowCount = 0;
    let nFontOrCrIndexFromEnd = -1;
    let tm = null;

    let cr = null;
    let font = null;

    for (var nLabel = 0; nLabel < nLabelCount; ++nLabel) {
      nFontOrCrIndexFromEnd = nLabel === 0 ? nLabel : arTextFonts.length-1 - (nLabelCount-1 - nLabel);
      font = arTextFonts[nFontOrCrIndexFromEnd];
      g.font = font;
      tm = g.measureText("W");
      nHFont = Math.abs(tm.actualBoundingBoxAscent) + tm.actualBoundingBoxDescent + 2*nYInset;

      if (nLabel === 0)
        nHTmp = Math.floor(nHFont/2);
      else
        nHTmp += nHFont;

      if (nHSum + nHFont > nHAvail)
        break;

      if (nLabel > 0)
        nHSumExcl0 += nHFont;

      nHSum += nHFont;
      ++nFittedRowCount;
    }

    let bHasRoom = false;
    let nFittedHeight = nHSum;
    let nAscent = -1;
    let nDescent = -1;
    let nYY = -1;
    let nWLabel = -1;
    let ob: any = null;
    let str = null;

    //draw first id first, then continue from the end
    let nFittedRow = 0;
    if(nFittedRowCount > 0) {
      g.font = arTextFonts[nFittedRow];
      tm = g.measureText("W");
      nAscent = Math.abs(tm.actualBoundingBoxAscent);
      nDescent = tm.actualBoundingBoxDescent;
      nHFont =  nAscent + nDescent + 2*nYInset;

      ob = arTextLabels[nFittedRow];
      str = ob === null ? "" : ob.toString();
      str = TextUtils.trimText(str, g, nW);
      if (bHasRoom)
        nYY = nY + Math.floor((nH + nHFont)/2);
      else {
        let nDeltaY = Math.floor((nHAvail - nFittedHeight)/2);
        nYY =  nY + nDeltaY + nHFont;
      }
      cr = arBackColors[nFittedRow];
      if (cr !== null) {
        g.fillStyle = cr;
        g.fillRect(nX, nYY - nHFont, nW, nHFont);
      }

      tm = g.measureText(str);
      nWLabel = tm.width;
      cr = arTextColors[nFittedRow];
      g.fillStyle = cr;
      g.fillText(str, nX + ((nW - nWLabel)>>1), nYY - nHFont + nYInset);
    }

    nYY = this.isGapToBottom() ? nY + nH : nYY + nHSumExcl0;

    for (nFittedRow = nLabelCount-1; nFittedRow >= nLabelCount-1 - nFittedRowCount+2; --nFittedRow) {
      nFontOrCrIndexFromEnd = arTextFonts.length-1 - (nLabelCount-1 - nFittedRow);
      font = arTextFonts[nFontOrCrIndexFromEnd];

      g.font = font;
      tm = g.measureText("W");
      nAscent = Math.abs(tm.actualBoundingBoxAscent);
      nDescent = tm.actualBoundingBoxDescent;
      nHFont = nAscent + nDescent + 2*nYInset;

      ob = arTextLabels[nFittedRow];

      if (typeof ob == "number") {
        cr = arBackColors[nFittedRow];
        g.fillStyle = cr;
        str = ob.toString();
      }
      else {
        str = ob == null ? "" : ob.toString();
        cr = arBackColors[nFittedRow];
        if(cr !== null) {
          g.fillStyle = cr;
          g.fillRect(nX, nYY - nHFont, nW, nHFont);
        }
      }

      str = TextUtils.trimText(str, g, nW);
      tm = g.measureText(str);
      nWLabel = tm.width;
      cr = arTextColors[nFontOrCrIndexFromEnd];
      g.fillStyle = cr;
      let nXX = nX + Math.floor((nW - nWLabel)/2);
      g.fillText(str, nXX, nYY - nHFont + nYInset);

      nYY -= nHFont;
    }
  }
}
