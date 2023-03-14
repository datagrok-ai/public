import * as DG from 'datagrok-api/dg';
import * as TextUtils from "../utils/TextUtils";
import {GridCellRendererEx} from "./GridCellRendererEx";

export class AbstractComposedCellRenderer extends GridCellRendererEx {

  private m_renderer: DG.GridCellRenderer | null = null;
  constructor(renderer: DG.GridCellRenderer | null = null) {
    super();

    this.m_renderer = renderer;
  }

  private renderExternalContent(g: CanvasRenderingContext2D, cellGrid: DG.GridCell, nX: number, nY: number, nW: number, nH: number): void {
    this.m_renderer!.render(g, nX, nY, nW, nH, cellGrid, cellGrid.style);
  }

  fillLabelsAndUI(cell: DG.Cell, arTextLabels : string[], arTextFonts: string[], arTextColors : string[], arBackColors: string[]) : void {
    throw new Error('Not Implemented');
  }

  render(g: CanvasRenderingContext2D, nX: number, nY: number, nW: number, nH: number, cellGrid: DG.GridCell, style: DG.GridCellStyle): void {
    if (nW <= 3 || nH <= 3)
      return;

    const cell: DG.Cell = cellGrid.cell;
    const value = cell.value;

    const arTextLabels : string[] = [];
    const arTextFonts : string[] = [];
    const arTextColors : string[] = [];
    const arBackColors : string[] = [];
    this.fillLabelsAndUI(value, arTextLabels, arTextFonts, arTextColors, arBackColors);

    const font = style.font;
    const nHFont = TextUtils.getFontSize(font);
    const nHFontSub = nHFont < 0 ? nHFont : nHFont - 2;

    //const arFonts = ["bold " + font, "bold " + TextUtils.setFontSize(font, nHFontSub)];

    try {
      this.paintLabelsAndExternal(g, cellGrid, arTextLabels, arTextFonts, arTextColors, arBackColors, nX, nY, nW, nH);
    } catch (e) {
      throw e; //provided for a debugger point
    }
  }

  paintLabelsAndExternal(g: CanvasRenderingContext2D, cellGrid: DG.GridCell, arIDs: string[], arFonts: string[], arTextColors: string[], arBackColors: string[],
                         nX: number, nY: number, nW: number, nH: number): void {
    g.fillStyle = "white";
    g.strokeStyle = "black";
    g.textAlign = "center";
    g.textBaseline = "top";

    g.fillRect(nX, nY, nW, nH);

    let nYInset = 2;
    let nHAvail = nH;
    let nWAvail = nW;
    const nIDCount = arIDs === null ? 0 : arIDs.length;

    g.font = arFonts[0];
    const tm = g.measureText("W");
    const nHFont = Math.abs(tm.actualBoundingBoxAscent) + tm.actualBoundingBoxDescent + 2 * nYInset;
    let nFittedRowCount = Math.floor(nHAvail / nHFont);//.getHeight());
    nFittedRowCount = nFittedRowCount < nIDCount ? nFittedRowCount : nIDCount;
    const nFittedHeight = nFittedRowCount * nHFont;//.getHeight();
    const nExtraHeight = nHAvail - nFittedHeight;

    const bFitAllIDs = nFittedRowCount >= nIDCount;
    const nHLabel = Math.floor(bFitAllIDs ? nHFont : (nFittedRowCount === 0 ? 0 : nHAvail / nFittedRowCount));

    const nCPdMinHeight = 21;

    let bCpd = true;
    if (nFittedRowCount < nIDCount || nExtraHeight < nCPdMinHeight)
      bCpd = false;

    if (nHLabel > 0) {

      let cr = null;
      let str = null;
      let nYTemp = null;
      for (let n = 0; n < nFittedRowCount; n++) {

        g.font = arFonts[n];
        str = arIDs[n];
        str = TextUtils.trimText(str, g, nW);

        nYTemp = nH - nHLabel * (nFittedRowCount - n);

        g.translate(0, nYTemp);
        cr = arBackColors[n];
        if (cr !== null) {
          g.fillStyle = cr;
          g.fillRect(nX, nY + nYInset, nW, nHFont);
        }

        g.fillStyle = arTextColors[n];
        g.fillText(str, nX + Math.floor(nW / 2), nY + nYInset /*+ nHFont*/);
        g.translate(0, -nYTemp);
      }
    }

    if (!bCpd)
      return;

    if (nHAvail > 5)
      this.renderExternalContent(g, cellGrid, nX, nY, nWAvail, nHAvail - nHLabel * nFittedRowCount);
  }
/*
  static paintLabels(g: CanvasRenderingContext2D, arIDs: string[], arBackColors: string[], nX: number, nY: number, nW: number, nH: number, crBack: string) {

    //g.font = "bold 13px Dialog";
    g.fillStyle = "black";
    g.strokeStyle = "black";
    g.textAlign = "center";
    g.textBaseline = "top";

    let nYInset = 2;

    let nHAvail = nH;
    let nWAvail = nW;

    const nIDCount = arIDs === null ? 0 : arIDs.length;

    const tm = g.measureText("W");
    const nHFont = Math.abs(tm.actualBoundingBoxAscent) + tm.actualBoundingBoxDescent + 2 * nYInset;
    let nFittedRowCount = Math.floor(nHAvail / nHFont);//.getHeight());
    nFittedRowCount = nFittedRowCount < nIDCount ? nFittedRowCount : nIDCount;
    const nFittedHeight = nFittedRowCount * nHFont;//.getHeight();
    const nExtraHeight = nHAvail - nFittedHeight;

    const bFitAllIDs = nFittedRowCount >= nIDCount;
    const nHLabel = Math.floor(bFitAllIDs ? nHFont : (nFittedRowCount === 0 ? 0 : nHAvail / nFittedRowCount));

    const nCPdMinHeight = 21;

    let bCpd = true;
    if (nFittedRowCount < nIDCount || nExtraHeight < nCPdMinHeight)
      bCpd = false;

    if (nHLabel > 0) {

      let cr = null;
      let str = null;
      let nYTemp = null;
      for (var n = 0; n < nFittedRowCount; n++) {
        str = arIDs[n];
        str = TextUtils.trimText(str, g, nW);

        nYTemp = nH - nHLabel * (nFittedRowCount - n);

        g.translate(0, nYTemp);
        cr = arBackColors[n];
        if (cr !== null) {
          g.fillStyle = cr;
          g.fillRect(nX, nY + nYInset, nW, nHFont);
        }

        g.fillStyle = "black";
        g.fillText(str, nX + Math.floor(nW / 2), nY + nYInset );
        g.translate(0, -nYTemp);
      }
    }

    if (!bCpd)
      return 0;

    return nHAvail - nHLabel * nFittedRowCount;
    //if (nHAvail > 5)
    //   this.paintStructure(g, entity, nX, nY, nWAvail, nHAvail - nHLabel * nFittedRowCount, crBack);
  }*/
}
