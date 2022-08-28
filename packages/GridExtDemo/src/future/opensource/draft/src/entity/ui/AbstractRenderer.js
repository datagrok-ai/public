import {Insets} from "../../geom/Insets";

export class AbstractRenderer
{
    constructor()
    {
       this.m_font = "13px Dialog";
    }


    getInsets() {return AbstractRenderer.ZERO_INSETS;}

    getFont()
    {
        return this.m_font;
    }

    setFont(font)
    {
        this.m_font = font;
    }


    /**
     * Returns the preferred width of the drawing area.This method is called by the application
     * framework to determine the size of the visual componnt for which the renderer is to be used.
     * @returns {number} the preferred drawing area's width in pixels.
     */
    getPreferredCellWidth()
    {
        const i = this.getInsets();
        return 75 + i.getL() + i.getR();
    }

    /**
     * Returns the preferred height of the drawing area.This method is called by the application
     * framework to determine the size of the visual componnt for which the renderer is to be used.
     * @returns {number} the preferred drawing area's width in pixels.
     */
    getPreferredCellHeight()
    {
        const i = this.getInsets();
        return 75 + i.getT() + i.getB();
    }


    getPopupGroupCount(cell) {return 0;}
    getPopupGroupName(cell, nGroup) {return "";}
    toPopupGrouptemsArray(cell, nGroup) {return [];}
    onPopupGroupAction(cell, nGroup, nItem) {}

    toPopupItemsArray(cell) {return [];}
    onPopupAction(cell, nItem) {}

    /**
     * Called by the application framework when the mouse is pressed on the area
     * rendered by this renderer. By default does nothing. Subclasses should provide
     * implementation to enable features like shouwing tooltips or property panels based on mouse hittest results.
     * @param cell the grid cell where the mouse was pressed.
     * @param nX the x-coordinate of the mouse relaative to the cell's origin.
     * @param nY the y-coordinate of the mouse relaative to the cell's origin.
     */
    onMousePressed(cell, nX, nY, nButton) {}

    /**
     * Called by the application framework when the mouse is released on the area
     * rendered by this renderer. By default does nothing. Subclasses should provide
     * implementation to enable features like shouwing tooltips or property panels based on mouse hittest results.
     * @param cell the grid cell where the mouse was pressed.
     * @param nX the x-coordinate of the mouse relaative to the cell's origin.
     * @param nY the y-coordinate of the mouse relaative to the cell's origin.
     */
    onMouseReleased(cell, nX, nY, nButton) {}

    /**
     * Called by the application framework when the mouse is clicked on the area
     * rendered by this renderer. By default does nothing. Subclasses should provide
     * implementation to enable features like shouwing tooltips or property panels based on mouse hittest results.
     * @param cell the grid cell where the mouse was pressed.
     * @param nX the x-coordinate of the mouse relaative to the cell's origin.
     * @param nY the y-coordinate of the mouse relaative to the cell's origin.
     */
    onMouseClicked(cell, nX, nY, nButton) {}

    /**
     * Called by the application framework when the mouse is moved on the area
     * rendered by this renderer. By default does nothing. Subclasses should provide
     * implementation to enable features like shouwing tooltips or property panels based on mouse hittest results.
     * @param cell the grid cell where the mouse was pressed.
     * @param nX the x-coordinate of the mouse relaative to the cell's origin.
     * @param nY the y-coordinate of the mouse relaative to the cell's origin.
     */
    onMouseMoved(cell, nX, nY, nButton) {}

    onResizeWidth(colGrid, grid, nW, bAdjusting)
    {
        //const arColRowIdxs = new Array(4);
        //GridUtils.fillVisibleGridCells(arColRowIdxs, grid);
       //console.log(colGrid.name + " " + nW + " " + bAdjusting);
    }
    onResizeHeight(colGrid, grid, nH, bAdjusting)
    {
        //console.log(nH + " " + bAdjusting);
    }

    paintBackground(g, nX, nY, nW, nH, crBack) {}
    paintBorder(g, nX, nY, nW, nH) {}

}

AbstractRenderer.ZERO_INSETS = new Insets(0,0,0,0);
AbstractRenderer.LIGHT_COLORS = ["LightBlue", "LightPink", "LightGreen", "LightSalmon", "Khaki", "PaleGreen"];
AbstractRenderer.USE_LIGHT_COLORS_COUNT = 0;
AbstractRenderer.generateLightColor = function()
{
    const n = AbstractRenderer.USE_LIGHT_COLORS_COUNT % AbstractRenderer.LIGHT_COLORS.length;
    ++AbstractRenderer.USE_LIGHT_COLORS_COUNT;
    return AbstractRenderer.LIGHT_COLORS[n];
}

AbstractRenderer.TAG_PRIMITIVE_COL_BACKGROUND = "PRIMITIVE_COL_BACKGROUND";



AbstractRenderer.toHtml = function(eParent, nW, nH, crBack)
{

   const eCanvas = ui.canvas();
   eCanvas.width  = nW;
   eCanvas.height = nH;
   eParent.appendChild(eCanvas);

   const g = eCanvas.getContext("2d");
   if(crBack === undefined)
     crBack = 'white';

   g.fillStyle = crBack;
   g.fillRect(0, 0, nW, nH);

   return eCanvas;
}
