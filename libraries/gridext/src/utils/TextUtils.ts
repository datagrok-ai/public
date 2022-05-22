
export function trimText(str : string, g : CanvasRenderingContext2D, nWidth : number) : string {
  let tm = g.measureText(str);
  let nW  = tm.width;
  if(nW <= nWidth)
    return str;

  //let nHFont = tm.actualBoundingBoxAscent + tm.actualBoundingBoxDescent;
  let strDots = "...";
  tm = g.measureText(strDots);
  let nWDots  = tm.width;
  if(nWDots > nWidth)
  {
    strDots = "..";
    tm = g.measureText(strDots);
    nWDots = tm.width;
    if(nWDots <= nWidth)
      return  strDots;

    strDots = ".";
    tm = g.measureText(strDots);
    nWDots = tm.width;
    return nWDots <= nWidth ? strDots : "";
  }

  let nWAvail = nWidth - nWDots;
  let strW = "W";
  tm = g.measureText(strW);
  let nWW = tm.width;

  let nCharCount = Math.floor(nWAvail / nWW);
  let strAdj = str.substring(0, nCharCount);
  tm = g.measureText(strAdj);
  if(tm.width > nWAvail)
  {
    for(var n=nCharCount -1; n>=0; --n)
    {
      strAdj = str.substring(0, n);
      tm = g.measureText(strAdj);
      if(tm.width <= nWAvail)
        return strAdj + strDots;
    }
  }
  else
  {
    let strAdjOld = strAdj;
    for(var n=nCharCount +1; n<str.length; ++n)
    {
      strAdj = str.substring(0, n);
      tm = g.measureText(strAdj);
      if(tm.width > nWAvail)
        return strAdjOld + strDots;

      strAdjOld = strAdj;
    }
  }
  return "...";
}
