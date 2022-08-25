/**
 * The SemEntityRenderer is na abstract class that defines a top level interface for semantic entities renderers
 * used across various visual components of the DataDrok platform such as viewers, tooltips, legends, property panels.
 * Subclasses must implement the actual rendering part
 */

import {AbstractRenderer} from "./AbstractRenderer";

export class SemEntityRenderer extends AbstractRenderer
{


    /**
     * Returns the preferred width of the drawing area.This method is called by the application
     * framework to determine the size of the visual componnt for which the renderer is to be used.
     * @returns {number} the preferred drawing area's width in pixels.
     */
    getPreferredCellWidth() {return 75;}

    /**
     * Returns the preferred height of the drawing area.This method is called by the application
     * framework to determine the size of the visual componnt for which the renderer is to be used.
     * @returns {number} the preferred drawing area's width in pixels.
     */
    getPreferredCellHeight() {return 75;}

    /**
     * Called by the application framework when the mouse is pressed on the area
     * rendered by this renderer. By default does nothing. Subclasses should provide
     * implementation to enable features like shouwing tooltips or property panels based on mouse hittest results.
     * @param nX the x-coordinate of top left corner of the drawing area.
     * @param nY the x-coordinate of top left corner of the drawing area.
     * @param nW the width of the drawing area.
     * @param nH the height of the drawing area.
     * @param e the DOM mouse pressed event.
     */
    onMousePressedOnHeader(nX, nY, nW, nH, e) {}

    /**
     * Called by the application framework when the mouse is released on the area
     * rendered by this renderer. By default does nothing. Subclasses should provide
     * implementation to enable features like shouwing tooltips or property panels based on mouse hittest results.
     * @param nX the x-coordinate of top left corner of the drawing area.
     * @param nY the x-coordinate of top left corner of the drawing area.
     * @param nW the width of the drawing area.
     * @param nH the height of the drawing area.
     * @param e the DOM mouse pressed event.
     */
    onMouseReleasedOnHeader(nX, nY, nW, nH, e) {}

    /**
     * Exports the rendered content into html DOM tree. The default implementation
     * create a single canvas element and output the graphics content onto it. Subclasses
     * can override this method to create a DOM tree consisting of various html elements.
     * @param eParent a parent DOM node to which the output DOM tree will be appended.
     * @param entity the entity object to be rendererd.
     * @param nW the width o the drawing area.
     * @param nH the height o the drawing area.
     * @param crBack the textual CSS representation for the background color ("white", "red")
     * @returns {HTMLElement} a reference to the root element of the output DOM tree.
     */
    toHtml(eParent, entity, nW, nH, crBack)
    {
        const eCanvas = AbstractRenderer.toHtml(eParent, nW, nH, crBack)
        const g = eCanvas.getContext("2d");

        this.paint(g, entity, 0, 0, nW, nH);

        return eCanvas;
    }

    /**
     *
     * @param g the graphics context to where the rendered content will be output.
     * @param entity the entity object to be rendererd.
     * @param nX the x-coordinate of the top left corner of the drawing area.
     * @param nY the y-coordinate of the top left corner of the drawing area.
     * @param nW the width o the drawing area.
     * @param nH the height o the drawing area.
     * @param crBack the textual CSS representation for the background color ("white", "red")
     */
    paint(g, entity, nX, nY, nW, nH, crBack)
    {
         if(crBack === null)
            return;

        if(crBack === undefined)
            crBack = "white";

        g.fillStyle = crBack;
        g.fillRect(nX, nY, nW, nH);

    }
}
