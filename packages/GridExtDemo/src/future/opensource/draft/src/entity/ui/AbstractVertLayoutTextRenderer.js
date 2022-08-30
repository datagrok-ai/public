import {SemEntityRenderer} from "./SemEntityRenderer";
import {TextUtils} from "../../utils/TextUtils";
import {NotImplementedError} from "../../lang/NotImplementedError";

export class AbstractVertLayoutTextRenderer extends SemEntityRenderer
{
    fillLabelsAndUI(entity, arTextLabels, arTextFonts, arTextColors, arBackColors)
    {
        throw new NotImplementedError();
    }

    paint(g, entity, nX, nY, nW, nH, crBack) {

        super.paint(g, entity, nX, nY, nW, nH, crBack);

        if(entity === null)
        {
            return;
        }

        const arTextLabels = [];
        const arTextFonts = [];
        const arTextColors = [];
        const arBackColors = [];

        if(crBack !== undefined)
        {
            g.fillStyle = crBack;
            g.fillRect(nX,nY, nW, nH);
        }

        this.fillLabelsAndUI(entity, arTextLabels, arTextFonts, arTextColors, arBackColors);
        this.paintLabels(g, arTextLabels, arTextFonts, arTextColors, arBackColors, nX, nY, nW, nH);
    }


    paintLabels(g, arTextLabels, arTextFonts, arTextColors, arBackColors, nX, nY, nW, nH) {
        g.font = this.getFont();//"12px Arial";
        g.fillStyle = "black";
        g.strokeStyle = "black";
        g.textAlign = "left";
        g.textBaseline = "top";

        const nYInset = 2;

        let nHAvail = nH;
        let nWAvail = nW;

        let nLabelCount = arTextLabels.length;

        let nHFont = -1;
        let nHSum = 0;
        let nHTmp = 0;
        let nFittedRowCount = 0;
        let nFontOrCrIndexFromEnd = -1;
        let tm = null;

        let cr = null;
        let font = null;

        for(var nLabel = 0; nLabel < nLabelCount; ++nLabel)
        {
            nFontOrCrIndexFromEnd = nLabel === 0 ? nLabel : arTextFonts.length-1 - (nLabelCount-1 - nLabel);
            font = arTextFonts[nFontOrCrIndexFromEnd];
            g.font = font;
            tm = g.measureText("W");
            nHFont = Math.abs(tm.actualBoundingBoxAscent) + tm.actualBoundingBoxDescent + 2*nYInset;

            if(nLabel === 0)
                nHTmp = Math.floor(nHFont/2);
            else
                nHTmp += nHFont;

            if(nHSum + nHFont> nHAvail)
                break;

            nHSum += nHFont;
            ++nFittedRowCount;
        }

        let bHasRoom = false;
        let nFittedHeight = nHSum;
        let nAscent = -1;
        let nDescent = -1;
        let nYY = -1;
        let nWLabel = -1;
        let ob = null;
        let str = null;

        //draw first id first, then continue from the end
        let nFittedRow = 0;
        if(nFittedRowCount > 0)
        {
            g.font = arTextFonts[nFittedRow];
            tm = g.measureText("W");
            nAscent = Math.abs(tm.actualBoundingBoxAscent);
            nDescent = tm.actualBoundingBoxDescent;
            nHFont =  nAscent + nDescent + 2*nYInset;

            ob = arTextLabels[nFittedRow];
            str = ob === null ? "" : ob.toString();
            str = TextUtils.trimText(str, g, nW);
            if(bHasRoom)
                nYY = nY + Math.floor((nH /2)+(nHFont/2));
            else
            {
                let nDeltaY = Math.floor((nHAvail - nFittedHeight)/2);
                nYY =  nY + nDeltaY + nHFont;
            }
            cr = arBackColors[nFittedRow];
            if(cr !== null)
            {
                g.fillStyle = cr;
                g.fillRect(nX, nYY - nHFont, nW, nHFont);
            }

            tm = g.measureText(str);
            nWLabel = tm.width;
            cr = arTextColors[nFittedRow];
            g.fillStyle = cr;
            g.fillText(str, nX + ((nW - nWLabel)>>1), nYY - nHFont + nYInset);
        }

        nYY = nY + nH;

        for(nFittedRow = nLabelCount-1; nFittedRow >= nLabelCount-1 - nFittedRowCount+2; --nFittedRow)
        {
            nFontOrCrIndexFromEnd = arTextFonts.length-1 - (nLabelCount-1 - nFittedRow);
            font = arTextFonts[nFontOrCrIndexFromEnd];

            g.font = font;
            tm = g.measureText("W");
            nAscent = Math.abs(tm.actualBoundingBoxAscent);
            nDescent = tm.actualBoundingBoxDescent;
            nHFont = nAscent + nDescent + 2*nYInset;

            ob = arTextLabels[nFittedRow];

            if(typeof ob == "number")
            {
                cr = arBackColors[nFittedRow];
                g.fillStyle = cr;
                let fPctEff = ob;

                if(fPctEff <= 100.0)
                {
                    let fK = fPctEff/100.0;
                    g.fillRect(nX, nYY - nHFont, Math.floor(nW*fK), nHFont);
                }
                else
                {
                    g.fillRect(nX, nYY - nHFont, nW, nHFont);
                    g.fillStyle = "white";
                    let fK = 100.0/fPctEff;
                    g.strokeLine(nX + Math.floor(nW*fK), nYY - nHFont, nX + Math.floor(nW*fK), nYY);
                }

                str = ob  + "% Eff";
            }
            else
            {
                str = ob == null ? "" : ob.toString();
                cr = arBackColors[nFittedRow];
                if(cr !== null)
                {
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
