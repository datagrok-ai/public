import {ErrorUtils} from "./ErrorUtils";

export class TextUtils
{
    constructor()
    {
        throw new Error("Cannot create instances of this class");
    }
}


TextUtils.buf2Str = function(buffer)
{
   let str = TextUtils.UTF8_DECODER.decode(buffer);
   return str;
}


TextUtils.trimText = function(str, ctx, nWidth) {
    let tm = ctx.measureText(str);
    let nW  = tm.width;
    if(nW <= nWidth)
        return str;

    let nHFont = tm.actualBoundingBoxAscent + tm.actualBoundingBoxDescent;
    let strDots = "...";
    tm = ctx.measureText(strDots);
    let nWDots  = tm.width;
    if(nWDots > nWidth)
    {
        strDots = "..";
        tm = ctx.measureText(strDots);
        nWDots = tm.width;
        if(nWDots <= nWidth)
            return  strDots;

        strDots = ".";
        tm = ctx.measureText(strDots);
        nWDots = tm.width;
        return nWDots <= nWidth ? strDots : "";
    }

    let nWAvail = nWidth - nWDots;
    let strW = "W";
    tm = ctx.measureText(strW);
    let nWW = tm.width;

    let nCharCount = Math.floor(nWAvail / nWW);
    let strAdj = str.substring(0, nCharCount);
    tm = ctx.measureText(strAdj);
    if(tm.width > nWAvail)
    {
        for(var n=nCharCount -1; n>=0; --n)
        {
            strAdj = str.substring(0, n);
            tm = ctx.measureText(strAdj);
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
            tm = ctx.measureText(strAdj);
            if(tm.width > nWAvail)
                return strAdjOld + strDots;

            strAdjOld = strAdj;
        }
     }
}

TextUtils.formatDate = function(nTime) {
    nTime = Math.floor(nTime);

    const date = new Date(nTime);
    const nMonth = date.getMonth() +1;
    const nDay = date.getDate();
    const nYear = date.getFullYear();

    let str = nMonth.toString() + "-" + nDay + "-" + nYear;//moment(new Date(nTime)).format('MM-DD-YYYY');
    return str;
}


TextUtils.splitTextIntoWords = function(strText)
{
    if (strText === null || strText.length === 0)
        return [];

    const maxLength = 4000;

    strText = strText.trim();
    const regx = /[^A-Za-z0-9_\-]+/g; // Non word chars except dash
    // String regx = "[^a-zA-Z0-9_\\-]";

    let s = strText.replaceAll(regx, " ").replaceAll(/\s+/g, " ").trim();
    if (s.length >= maxLength)
        s = s.substring(0, maxLength).trim();

    const words = s.split(" ");
    return words;
}


TextUtils.formatNumber = function(fValue)
{
    ErrorUtils.verifyType(fValue, Number);

    if(fValue === DG.FLOAT_NULL)
        return "";

    if(fValue === Math.floor(fValue))
    {
        return fValue.toString();
    }


    let str = "";
    if(fValue < 1)
    {
        let nPow = 1 + Math.abs(Math.floor(Math.log10(Math.abs(fValue))));
        if(nPow> 3)
            nPow = 3;

        str = fValue.toFixed(nPow);
         //console.log("pow: " + nPow);
    }
    else if(fValue < 1000)
    {
        if(Math.floor(fValue) === fValue)
            str = fValue.toFixed(0);
        else
            str = fValue.toFixed(1);
    }
    else str = fValue.toFixed(0);

    return str;
}


TextUtils.getFontSize = function(strFont)
{
    ErrorUtils.verifyType(strFont, String);

    const nIdxTo = strFont.indexOf("px");
    if(nIdxTo < 0)
        return -1;

    let nFr = nIdxTo-1;
    for(; nFr>=0; --nFr)
    {
        if(strFont.charAt(nFr) === " ")
         break;
    }


    const str = strFont.substring(nFr+1, nIdxTo);
    const nH =  parseInt(str);
    return  nH;
}

TextUtils.setFontSize = function(strFont, nH) {
    ErrorUtils.verifyType(strFont, String);
    ErrorUtils.verifyType(nH, Number);

    const nIdxTo = strFont.indexOf("px");
    if(nIdxTo < 0)
        return false;

    let nFr = nIdxTo-1;
    for(; nFr>=0; --nFr)
    {
        if(strFont.charAt(nFr) === " ")
            break;
    }

    strFont = strFont.substring(0, nFr+1) + nH.toString() +  strFont.substring(nIdxTo);
    return strFont;
}



// Define recursive function to print nested values
TextUtils.printValues = function(obj) {
    for(var k in obj) {
        if(obj[k] instanceof Object) {
            console.log("Node " + k + " " + obj[k]);
            if(k === 'children')
            {
                const arChildren = obj[k];
                for(let n=0; n<arChildren.length; ++n) {
                    const e = arChildren[n].state.element;
                    let strType = "";
                }
            }

            TextUtils.printValues(obj[k]);
        } else {
           console.log(k + " " + obj[k]);
        }
    }
}


TextUtils.UTF8_DECODER = new TextDecoder("utf-8");
