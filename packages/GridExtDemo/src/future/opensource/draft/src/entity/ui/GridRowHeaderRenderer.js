import {AbstractRenderer} from "./AbstractRenderer";
import {TextUtils} from "../../utils/TextUtils";

export class GridRowHeaderRenderer extends AbstractRenderer
{

    getPreferredCellWidth()
    {
        const i = this.getInsets();
        return 55 + i.getL() + i.getR();
    }

    getPreferredCellHeight()
    {
        const i = this.getInsets();
        return 20 + i.getT() + i.getB();
    }


    /**
     * Exports the rendered content into html DOM tree. The default implementation
     * create a single canvas element and output the graphics content onto it. Subclasses
     * can override this method to create a DOM tree consisting of various html elements.
     * @param eParent a parent DOM node to which the output DOM tree will be appended.
     * @param nRowGrid the index of the grid's row identifying the cell to be painted.
     * @param grid the grid viwer.
     * @param nW the width o the drawing area.
     * @param nH the height o the drawing area.
     * @param crBack the textual CSS representation for the background color ("white", "red")
     * @returns {HTMLElement} a reference to the root element of the output DOM tree.
     */
    toHtml(eParent, nRowGrid, grid, nW, nH, crBack)
    {
        const eCanvas = AbstractRenderer.toHtml(eParent, nW, nH, crBack)
        const g = eCanvas.getContext("2d");

        this.paint(g, nRowGrid, grid, 0, 0, nW, nH);

        return eCanvas;
    }



    paint(g, nRowGrid, nRowTable, grid, nX, nY, nW, nH, crBack) {
        this.paintBackground(g, nX, nY, nW, nH);
        this.paintBorder(g, nX, nY, nW, nH);
        this.paintContent(g, nRowGrid, nRowTable, grid, nX, nY, nW, nH, crBack);
    }

    paintContent(g, nRowGrid, nRowTable, grid, nX, nY, nW, nH, crBack)
    {
        let str = nRowGrid.toString();
        str = TextUtils.trimText(str, g, nW);
        g.font = "13px Dialog";
        const tm = g.measureText(str);
        const nWLabel = tm.width;

        const nAscent = Math.abs(tm.actualBoundingBoxAscent);
        const nDescent = tm.actualBoundingBoxDescent;
        const nHFont =  nAscent + nDescent;// + 2*nYInset;

        g.textBaseline = "alphabetic";
        g.fillStyle = "Gray";
        g.fillText(str, nX + ((nW - nWLabel)>>1), nY + (nH - nHFont)/2 + nAscent);
    }
}
