
export function getFontSize(strFont : string) : number {
  const nIdxPX = strFont.indexOf("px", 0);
  if(nIdxPX < 0) {
    return -1;
  }

  let strSize = "";
  for(let n=nIdxPX-1; n >=0;  --n) {
    if(strFont.charAt(n) === ' ') {
      strSize = strFont.substring(n+1, nIdxPX);
      const nSize = parseInt(strSize);
      return nSize;
    }
  }

  strSize = strFont.substring(0, nIdxPX);
  const nSize = parseInt(strSize);
  return nSize;
}

export function setFontSize(strFont : string, nSize : number) : string {
  const nIdxPX = strFont.indexOf("px", 0);
  if(nIdxPX < 0) {
    return strFont;
  }

  let strSize = "";
  for(let n=nIdxPX-1; n >=0;  --n) {
    if(strFont.charAt(n) === ' ') {
      strSize = strFont.substring(n+1, nIdxPX);
      const str = strFont.replace(strSize, nSize.toString());
      return str;
    }
  }

  strSize = strFont.substring(0, nIdxPX);
  const str = strFont.replace(strSize, nSize.toString());
  return str;
}



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
