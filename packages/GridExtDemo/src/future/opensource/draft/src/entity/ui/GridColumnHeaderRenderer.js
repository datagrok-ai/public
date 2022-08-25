import {AbstractRenderer} from "./AbstractRenderer";
import {DGApp} from "../app/DGApp";

export class GridColumnHeaderRenderer extends AbstractRenderer
{
    constructor()
    {
        super();

        this.m_bPaintBorder = true;
        this.m_bPaintBackground = true;
    }

    getPaintBorder() {return this.m_bPaintBorder;}
    setPaintBorder(bPaint)
    {
        this.m_bPaintBorder = bPaint;
    }

    getPaintBackground() {return this.m_bPaintBackground;}
    setPaintBackground(bPaint)
    {
        this.m_bPaintBackground = bPaint;
    }

    adjustColumn(cell)
    {
        return cell.tableColumn;
    }

    createTootipContent(cell)
    {
        const nW = Math.floor(1.3*this.getPreferredCellWidth());
        const nH = Math.floor(1.3*this.getPreferredCellHeight());
        const eCanvas = ui.canvas(nW, nH);
        const g = eCanvas.getContext('2d');
        this.paint(g, cell.gridColumn, 0,0, nW, nH);
        return eCanvas;
    }


    /**
     * Exports the rendered content into html DOM tree. The default implementation
     * create a single canvas element and output the graphics content onto it. Subclasses
     * can override this method to create a DOM tree consisting of various html elements.
     * @param eParent a parent DOM node to which the output DOM tree will be appended.
     * @param column the grid's column.
     * @param nW the width o the drawing area.
     * @param nH the height o the drawing area.
     * @param crBack the textual CSS representation for the background color ("white", "red")
     * @returns {HTMLElement} a reference to the root element of the output DOM tree.
     */
    toHtml(eParent, column, nW, nH, crBack)
    {
        const eCanvas = AbstractRenderer.toHtml(eParent, nW, nH, crBack)
        const g = eCanvas.getContext("2d");

        this.paint(g, column, 0, 0, nW, nH);

        return eCanvas;
    }

    paint(g, colGrid, nX, nY, nW, nH, crBack)
    {
        const col = colGrid === null ? null : colGrid.column;
        if(col !== null && !DGApp.isVirtual(col)) {

            const crBackPrimi = colGrid.column.getTag(AbstractRenderer.TAG_PRIMITIVE_COL_BACKGROUND);
            if(crBackPrimi !== null && crBackPrimi !== undefined)
                crBack = crBackPrimi;
        }

        if(this.m_bPaintBackground)
         this.paintBackground(g, nX, nY, nW, nH, crBack);

        if(this.m_bPaintBorder)
         this.paintBorder(g, nX, nY, nW, nH);
    }

}

