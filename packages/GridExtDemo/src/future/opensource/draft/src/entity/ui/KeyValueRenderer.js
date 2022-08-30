import {SemEntityRenderer} from "./SemEntityRenderer";
import {TextUtils} from "../../utils/TextUtils";
import {NotImplementedError} from "../../lang/NotImplementedError";
import {ErrorUtils} from "../../utils/ErrorUtils";

export class KeyValueRenderer extends SemEntityRenderer {

    constructor() {
        super();

        this.m_fKeys2ValsWidthRatio = 0.5;
    }

    preparePairsAndUI(entity, arKeys, arVals, arKeysFonts, arKeysForeColors, arKeysBackColots,
                      arValsFonts, arValsForeColors, arValsBackColors)
    {
        throw new NotImplementedError();
    }

    getKeys2ValsWidthRatio()
    {
        return this.m_fKeys2ValsWidthRatio;
    }

    setKeys2ValsWidthRatio(fRatio)
    {
        ErrorUtils.verifyType(fRatio, Number);
        if(fRatio < 0 || fRatio > 1.0)
            throw new Error("Ratio " + fRatio + " is out of range [0,1].");

        this.m_fKeys2ValsWidthRatio = fRatio;
    }

    paint(g, entity, nX, nY, nW, nH, crBack) {

        super.paint(g, entity, nX, nY, nW, nH, crBack);

        const arKeys = [];
        const arVals = [];
        const arKeysFonts = [];
        const arKeysForeColors = [];
        const arKeysBackColots = [];
        const arValsFonts = [];
        const arValsForeColors = [];
        const arValsBackColors = [];

        this.preparePairsAndUI(entity, arKeys, arVals, arKeysFonts, arKeysForeColors, arKeysBackColots,
            arValsFonts, arValsForeColors, arValsBackColors);
        this.paintPairs(g, arKeys, arVals, arKeysFonts, arKeysForeColors, arKeysBackColots,
            arValsFonts, arValsForeColors, arValsBackColors, nX, nY, nW, nH);
    }

    paintPairs(g, arKeys, arVals, arKeysFonts, arKeysForeColors, arKeysBackColots,
                          arValsFonts, arValsForeColors, arValsBackColors, nX, nY, nW, nH)
    {
        if(arKeys.length === 0)
            return;

        g.font = this.getFont();
        g.textAlign = "left";


        const nYInset = 2;
        //const nW05 = Math.floor(nW * this.m_fKeys2ValsWidthRatio);
        const nWKets = Math.floor(nW * this.m_fKeys2ValsWidthRatio);
        const nWVals = nW - nWKets;

        let nWAvailHalf = nW;
        let nHAvail = nH;
        ;

        let tm = g.measureText("W");
        let nHFont = tm.actualBoundingBoxAscent + tm.actualBoundingBoxDescent + 2 * nYInset;

        let nFittedRowCount = nHFont === 0 ? 10000 : Math.floor(nHAvail / nHFont);
        nFittedRowCount = nFittedRowCount < arKeys.length ? nFittedRowCount : arKeys.length;
        const nFittedHeight = nFittedRowCount * nHFont;
        const nExtraHeight = nHAvail - nFittedHeight;

        let str = null;
        let cr = null;
        let nYY = -1;
        for (let nFittedRow = 0; nFittedRow < nFittedRowCount; ++nFittedRow) {

            //Key
            g.fillStyle = arKeysBackColots[nFittedRow];
            g.fillRect(nX, nY + (nExtraHeight / 2) + nFittedRow * nHFont, nWKets, nHFont);
            g.strokeStyle = "white";
            g.strokeRect(nX, nY + Math.floor(nExtraHeight / 2) + nFittedRow * nHFont, nWKets, nHFont);

            nYY = Math.floor(nExtraHeight / 2) + (nFittedRow + 1) * nHFont - tm.actualBoundingBoxDescent;
            str = TextUtils.trimText(arKeys[nFittedRow], g, nWKets);

            g.font = arKeysFonts[nFittedRow];
            g.fillStyle = arKeysForeColors[nFittedRow];
            g.fillText(str, nX, nY + nYY - nYInset);


            //Value
            str = arVals[nFittedRow];
            str = TextUtils.trimText(str, g, nWVals);
            g.fillStyle = arValsBackColors[nFittedRow];
            g.fillRect(nX + nWKets, nY + Math.floor(nExtraHeight / 2) + nFittedRow * nHFont, nWVals, nHFont);
            g.strokeStyle = "white";
            g.strokeRect(nX + nWKets, nY + Math.floor(nExtraHeight / 2) + nFittedRow * nHFont, nWVals, nHFont);

            g.font = arValsFonts[nFittedRow];
            g.fillStyle = arValsForeColors[nFittedRow];
            g.fillText(str, nX + nWKets, nY + nYY - nYInset);
        }
    }
}
