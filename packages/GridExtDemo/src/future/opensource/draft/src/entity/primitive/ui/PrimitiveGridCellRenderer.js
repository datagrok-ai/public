//import * as DG from "datagrok-api/dg";
import {GridCellRenderer} from "../../ui/GridCellRenderer";
import {MathUtils} from "../../../utils/MathUtils";
import {TextUtils} from "../../../utils/TextUtils";
import {ErrorUtils} from "../../../utils/ErrorUtils";

export class PrimitiveGridCellRenderer extends GridCellRenderer
{
    constructor(rendererMol)
    {
        super();

        this.m_rendererMol = rendererMol === undefined ? null : rendererMol;
    }

    getPreferredCellHeight()
    {
        const i = this.getInsets();
        return 25 + i.getT() + i.getB();
    }

    paint(g, cell, nX, nY, nW, nH, crBack)
    {
        super.paint(g, cell, nX, nY, nW, nH, crBack);

        const style = cell.style;
        const nRow = cell.tableRow.idx;
        const col = cell.tableColumn;
        const value = col.get(nRow);
        const type = col.type;

        if(MathUtils.isNullValue(value))
            return false;

        const font = style.font;

        g.font = font;

        const nYInset = 2;
        let str = "";

        if(type === DG.COLUMN_TYPE.FLOAT)
         str = TextUtils.formatNumber(value);
        else if(type === DG.COLUMN_TYPE.QNUM)
        {
         str = DG.Qnum.toString(value)
        }

        else if(type === DG.COLUMN_TYPE.DATE_TIME) {
            const nTime = value.a;

            if (!isNaN(nTime))
                str = TextUtils.formatDate(nTime);
        }

        else if(type === DG.COLUMN_TYPE.STRING)
        {
            ErrorUtils.verifyType(value, String);
            const typeSemDG = col.semType;
            if(typeSemDG !== null && typeSemDG === DG.SEMTYPE.MOLECULE && this.m_rendererMol !== null)
            {
                window.devicePixelRatio = 1;
                g.save();
                //cellstyle.backColor = DG.Color.white;
                this.m_rendererMol.render(g, nX, nY, nW, nH, cell, cell.style);
                g.restore();
                return;
            }



            const ar  = value.split("\n");
            if(ar.length > 0)
            {
                let fValue = NaN;
                for(let n=0; n<ar.length; ++n)
                {
                    fValue = Number(ar[n]);
                    if(!isNaN(fValue))
                     ar[n] = TextUtils.formatNumber(fValue);
                }

                PrimitiveGridCellRenderer.paintVertLayout(g, ar, nX, nY, nW, nH, font);
            }

        }


        else str = value.toString();

        //str = "test " + nRow;
        str = TextUtils.trimText(str, g, nW);
        let tm = g.measureText(str);
        const nWLabel = tm.width;


        tm = g.measureText("W");
        const nAscent = Math.abs(tm.actualBoundingBoxAscent);
        const nDescent = tm.actualBoundingBoxDescent;
        const nHFont =  nAscent + nDescent + 2*nYInset;

        const nDeltaY = Math.floor((nH - nHFont)/2);
        const nYY =  nY + nDeltaY + nHFont;
        g.fillStyle = "black";
        g.textBaseline = "top";
        g.fillText(str, nX + ((nW - nWLabel)>>1), nYY - nHFont + nYInset);
        g.textBaseline;
        return false;
    }
}


PrimitiveGridCellRenderer.paintVertLayout = function(g, arLabels, nX, nY, nW, nH, font)
{
    if(arLabels.length === 0)
        return;

    g.font = font;
    g.textAlign = "left";
    g.fillStyle = "black";

    const nYInset = 2;
    let nW05 = nW;//Math.floor(nW / 2);

    let nWAvailHalf = nW;
    let nHAvail = nH;
    ;

    let tm = g.measureText("W");
    const nHFont = tm.actualBoundingBoxAscent + tm.actualBoundingBoxDescent + 2 * nYInset;

    let nFittedRowCount = nHFont === 0 ? 10000 : Math.floor(nHAvail / nHFont);
    nFittedRowCount = nFittedRowCount < arLabels.length ? nFittedRowCount : arLabels.length;
    const nFittedHeight = nFittedRowCount * nHFont;
    const nExtraHeight = nHAvail - nFittedHeight;

    let nWTextMax = -1;
    let str = null;
    let nFittedRow = 0;
    for (nFittedRow = 0; nFittedRow < nFittedRowCount; ++nFittedRow) {

        str = arLabels[nFittedRow];
        str = TextUtils.trimText(str, g, nW05);
        arLabels[nFittedRow] = str;
        tm = g.measureText(str);

        if(nWTextMax === -1 || tm.width > nWTextMax)
            nWTextMax = tm.width;
    }

    const nXoffset = Math.floor((nW - nWTextMax) / 2);
    let nYY = -1;
    for (nFittedRow = 0; nFittedRow < nFittedRowCount; ++nFittedRow) {

        nYY = Math.floor(nExtraHeight / 2) + (nFittedRow + 1) * nHFont - tm.actualBoundingBoxDescent;
        str =  arLabels[nFittedRow];
        g.fillText(str, nX + nXoffset, nY + nYY - nYInset);
    }
}
