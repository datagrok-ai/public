import {SemEntityRenderer} from "../../../ui/SemEntityRenderer";
import {TextUtils} from "../../../../utils/TextUtils";
import {NotImplementedError} from "../../../../lang/NotImplementedError";
import {MathUtils} from "../../../../utils/MathUtils";
import {DateUtils} from "../../../../utils/DateUtils";

export class AbstractCpdEntityRenderer extends SemEntityRenderer
{

    getPreferredCellWidth()
    {
        const i = this.getInsets();
        return 160 + i.getL() + i.getR();
    }

    getPreferredCellHeight()
    {
        const i = this.getInsets();
        return 95 + i.getT() + i.getB();
    }


    paintStructure(g, entity, arIDs, nX, nY, nW, nH, crBack) {throw new NotImplementedError();}

    paint(g, entity, nX, nY, nW, nH, crBack) {

        super.paint(g, entity, nX, nY, nW, nH, crBack);

        const strSmiles = entity.getSmiles();
        const strGnfId = entity.getGnfRegId();
        const strNVPId = entity.getSampleId();


        let strConceptId = entity.getConceptId();
        if(strConceptId === undefined)
           strConceptId = "";

        const strSIDId = "SID" + entity.getCpdSidId().toString();
        const nTime = entity.getLastTested();

        const strLastTested = MathUtils.isNullValue(nTime) || nTime < 0 ? "" : TextUtils.formatDate(nTime);

        const arBackColors = [null, DateUtils.getRecentDateColor(nTime)];

        const font = this.getFont();
        const nHFont = TextUtils.getFontSize(font);
        const nHFontSub = nHFont < 0 ? nHFont : nHFont -2;

        const arFonts = ["bold " + font, "bold " + TextUtils.setFontSize(font, nHFontSub)];

        try{this.paintLabelsAndStruct(g, entity, [strConceptId, strLastTested], arFonts, arBackColors, nX, nY, nW, nH, crBack);}
        catch(e)
        {
            throw e; //provided for a debugger point
        }
    }



    paintLabelsAndStruct(g, entity, arIDs, arFonts, arBackColors, nX, nY, nW, nH, crBack) {


        //const font = this.getFont();
        //g.font = "bold " + font;// "bold 13px Dialog";
        g.fillStyle = "black";
        g.strokeStyle = "black";
        g.textAlign = "center";
        g.textBaseline = "top";

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

            let cr  = null;
            let str = null;
            let nYTemp = null;
            for (let n = 0; n < nFittedRowCount; n++) {

                g.font = arFonts[n];
                str = arIDs[n];
                str = TextUtils.trimText(str, g, nW);

                nYTemp = nH - nHLabel * (nFittedRowCount - n);

                g.translate(0, nYTemp);
                cr = arBackColors[n];
                if(cr !== null)
                {
                    g.fillStyle = cr;
                    g.fillRect(nX, nY + nYInset, nW, nHFont);
                }

                g.fillStyle = "black";
                g.fillText(str, nX + Math.floor(nW / 2), nY + nYInset /*+ nHFont*/);
                g.translate(0, -nYTemp);
            }
        }

        if (!bCpd)
            return;

        if (nHAvail > 5)
            this.paintStructure(g, entity, nX, nY, nWAvail, nHAvail - nHLabel * nFittedRowCount, crBack);
    }
}

AbstractCpdEntityRenderer.paintLabels = function(g, entity, arIDs, arBackColors, nX, nY, nW, nH, crBack) {

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

        let cr  = null;
        let str = null;
        let nYTemp = null;
        for (var n = 0; n < nFittedRowCount; n++) {
            str = arIDs[n];
            str = TextUtils.trimText(str, g, nW);

            nYTemp = nH - nHLabel * (nFittedRowCount - n);

            g.translate(0, nYTemp);
            cr = arBackColors[n];
            if(cr !== null)
            {
                g.fillStyle = cr;
                g.fillRect(nX, nY + nYInset, nW, nHFont);
            }

            g.fillStyle = "black";
            g.fillText(str, nX + Math.floor(nW / 2), nY + nYInset /*+ nHFont*/);
            g.translate(0, -nYTemp);
        }
    }

    if (!bCpd)
        return 0;


    return nHAvail - nHLabel * nFittedRowCount;
    //if (nHAvail > 5)
     //   this.paintStructure(g, entity, nX, nY, nWAvail, nHAvail - nHLabel * nFittedRowCount, crBack);
}

