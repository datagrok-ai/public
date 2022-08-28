import {GridUtils} from "../../utils/GridUtils";
import {AbstractRenderer} from "./AbstractRenderer";
import {TextUtils} from "../../utils/TextUtils";

/**
 * The FlowTextRenderer class provides implementation for a renderer that layouts
 * a long text into readable lines to be fit into a size-limited drawing area.
 */
export class FlowTextRenderer extends AbstractRenderer
{
    constructor()
    {
        super();

        this.m_bVertText = false;
        this.m_nGapLines = 2;
    }


    getLinesGap() {return this.m_nGapLines;}

    isVerticalTextMode()
    {
        return this.m_bVertText;
    }

    setVerticalTextMode(bVertical)
    {
        this.m_bVertText = bVertical;
    }


    paint(g, strText, nX, nY, nW, nH, crBack)
    {
        const arLines = [];

        g.fillStyle = "black";
        g.strokeStyle = "black";
        g.textAlign = "left";
        g.textBaseline = "top";
        g.font = this.getFont();

        const nYInsets = this.getLinesGap();

        let nXTmp = nX;
        let nYTmp = nY;
        let nWTmp = nW;
        let nHTmp = nH;

        let bSkipLayout = false;
        if(this.m_bVertText)
        {
            const tm = g.measureText("W");
            const nHFont = Math.abs(tm.actualBoundingBoxAscent) + tm.actualBoundingBoxDescent + nYInsets;
            if(2*nHFont > nW)   //only one line can be fit
            {
             const strLine = TextUtils.trimText(strText, g, nH);
             arLines.push(strLine);
              bSkipLayout = true;
            }

            g.save();
            nXTmp = 0;
            nYTmp = 0;
            nWTmp = nH;
            nHTmp = nW;
            g.translate(nX, nY + nH);
            g.rotate(-Math.PI/2.0);
        }

        if(!bSkipLayout)
         GridUtils.calcWordsCellLayout(g, arLines, strText, nWTmp, nHTmp);

        this.paintLines(g, arLines, nXTmp, nYTmp, nWTmp, nHTmp, nYInsets);

        if(this.m_bVertText)
          g.restore();

        return true;
    }

    /**
     * Paint the textual content of a cell whose words have a specified layout.
     * @param g the canvas' graphics context. Should be initialized with font and alignment flags before calling this function.
     * @param arLines the specified array that contains the words' layout.
     * @param nX the starting x-coordinate to paint text
     * @param nY the starting x-coordinate to paint text
     * @param nW he width of the cell's textual area
     * @param nH the height of the cell's textual area
     * @param nYInsets the specified gap between the text lines. Usually 1-2 pixels.
     */
    paintLines(g, arLines, nX, nY, nW, nH, nYInsets) {
        g.font = this.getFont();

        const tm = g.measureText("W");
        const nHFont = Math.abs(tm.actualBoundingBoxAscent) + tm.actualBoundingBoxDescent + nYInsets;

        let str = "";
        let nWLine = 0;
        let nHText = 0;
        for(let nLine=0; nLine<arLines.length; ++nLine)
        {
            if(nHText + nHFont> nH)
                break;

            nWLine = nX;
            for (let nWord = 0; nWord < arLines[nLine].length; ++nWord)
            {
                str = arLines[nLine][nWord];
                g.fillText(str, nWLine, nY + nYInsets + nHFont*(nLine/* +1*/));
                nWLine += g.measureText(str).width;
            }
            nHText += nHFont;
        }
    }
}
