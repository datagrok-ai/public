import * as TextUtils from "../utils/TextUtils";

export function renderXYCenteredText(str : string, g : CanvasRenderingContext2D, nX : number, nY : number, nW : number, nH : number, strFont : string, strCrFore : string) : void {
  str = TextUtils.trimText(str, g, nW);

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
  const strBaseOld = g.textBaseline;
  g.textBaseline = 'top';
  const nXX = nX + Math.floor((nW - nWLabel) / 2);

  g.textAlign = 'start';
  g.fillStyle = strCrFore !== undefined ? strCrFore : 'black';
  g.fillText(str, nXX, nYY - nHFont + nYInset);
  g.textBaseline = strBaseOld;
}


export function renderXCenteredText(str : string, g : CanvasRenderingContext2D, nX : number, nY : number, nW : number, nH : number, strFont : string, strCrFore : string) : void {

  if (strFont !== null && strFont !== undefined && strFont !== '') {
    g.font = strFont;
  }
  let tm = g.measureText(str);
  const nWLabel = Math.round(tm.width);
 // const nYInset = 2;
  tm = g.measureText('W');
  const nAscent = Math.abs(tm.actualBoundingBoxAscent);
  const nDescent = tm.actualBoundingBoxDescent;
  const nHFont : number = nAscent + nDescent;// + 2 * nYInset;

  //const nDeltaY : number = Math.floor((nH - nHFont) / 2);
  //const nYY = nY + nDeltaY + nHFont;
  const strBaseOld = g.textBaseline;
  g.textBaseline = 'top';
  const nXX = nX + Math.floor((nW - nWLabel) / 2);

  g.textAlign = 'start';
  g.fillStyle = strCrFore !== undefined ? strCrFore : 'black';
  g.fillText(str, nXX, nY /*- nHFont*/);
  g.textBaseline = strBaseOld;
}
