import {AbstractRenderer} from "./AbstractRenderer";
import {GridUtils} from "../../utils/GridUtils";
import {SemType} from "../SemType";
import {DGApp} from "../app/DGApp";
import {FiltersUtils} from "../../ui/filters/FiltersUtils";
import {GridRowHeader, GridRowHeaderEmpty} from "../../ui/grid/GridRowHeader";
import * as DG from "datagrok-api/dg";
import {TextUtils} from "../../utils/TextUtils";

export class GridCellRenderer extends AbstractRenderer
{
    adjustColumnName(strColName)
    {
        return strColName;
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
        this.paint(g, cell, 0,0, nW, nH);
        return eCanvas;
    }


    getPopupGroupCount(cell) {return 2;}
    getPopupGroupName(cell, nGroup)
    {
        return nGroup === 0 ? "Sort" : "Columns Names Display";
    }

    toPopupGrouptemsArray(cell, nGroup)
    {
        if(nGroup === 0) {
            const col = this.adjustColumn(cell);
            const typeSem = GridUtils.getSemType(col);
            const arItems = SemType.toSortOptionNamesArray(typeSem);
            return arItems;
        }
        else
        {
            return [DGApp.USE_ALT_COL_NAMES ? "Use Default Column Names" : "Use Alt Column Names"];
        }
    }

    onPopupGroupAction(cell, nGroup, nItem)
    {
        if(nGroup === 0) {
            const nOption = nItem;
            const nRecord = cell.tableRowIndex;

            GridUtils.sort(cell, nRecord, nOption);
        }
        else
        {
            DGApp.USE_ALT_COL_NAMES = !DGApp.USE_ALT_COL_NAMES;
            const grid = cell.grid;
            GridUtils.repaint(grid);
        }
    }

    toPopupItemsArray(cell)
    {
     const grid = cell.grid;
     const col = this.adjustColumn(cell);
     const headerRows = GridUtils.getRowsHeader(grid, cell.gridColumn.name);

     const bTooltipOn = GridUtils.isTooltipEnabled(grid);

     let strTooltipOnOff = "Turn Tooltip On";
     if(bTooltipOn)
         strTooltipOnOff = "Turn Tooltip Off";


     let strFreezeUnfreeze = "Unfreeze";
     if(headerRows === null)
         strFreezeUnfreeze = "Freeze";
     /*
     else
     {
         const strColName = headerRows.getColumnName();
         if(strColName !== col.name)
          strFreezeUnfreeze = "Freeze";
     } */

      const arItems = [strTooltipOnOff, strFreezeUnfreeze, "Filter"]; //there is always a freeze/unfreeze filter


        const viewerFilter = FiltersUtils.getFiltersViewer();
        if(viewerFilter !== null) {
            const panelFilters = FiltersUtils.getFiltersPanel(col.dataFrame);
            const filter = panelFilters.getFilter(col.name);
            if (filter !== null)
                arItems.push("Close Filter");
        }

        const typeSem = GridUtils.getSemType(col);
        if(typeSem instanceof SemType && DGApp.isVirtual(col)) {
            arItems.push("Expand to Primitives");
            return arItems;
        }

        const strParentColName = col.getTag(DGApp.VIRTUAL_PARENT_TAG_NAME);
        if(strParentColName !== null && strParentColName !== undefined) {
            arItems.push("Collapse to Object");
        }

        return arItems;
    }

    onPopupAction(cell, nItem)
    {
        const grid = cell.grid;
        if(nItem === 0)
        {
            const bTooltipOn = GridUtils.isTooltipEnabled(grid);
            GridUtils.enableTooltip(grid, !bTooltipOn);
        }

        else if(nItem === 1)  //Freeze/Unfreeze
        {
            const strColName = cell.gridColumn.name;
            let headerRows = GridUtils.getRowsHeader(cell.grid, strColName);
            if(headerRows === null) {
                headerRows = new GridRowHeader(cell.grid.view, strColName);
/*
                let ar = DG.Viewer.getViewerTypes();
                let v = null;

                try{

                    v = DG.Viewer.fromType("emptyviewer", cell.grid.dataFrame);
                    cell.grid.view.addViewer(v);}
                catch(e) {//my changes
                    let hkhj = 0;
                }*/
            }
            else
            {
                headerRows.close();
            }

            /*
            const strColName = headerRows.getColumnName();
            if(strColName === cell.gridColumn.name)
                headerRows.setColumnName(null);
            else
                headerRows.setColumnName(cell.gridColumn.name);*/
            return;
        }

        if(nItem === 2)
        {
            //alert("Show FIlter for column " + colGrid.name);
            FiltersUtils.floatFilterToColumn(this, cell);
            return;
        }

        if(nItem === 3)
        {
            const col = this.adjustColumn(cell);
            const panelFilters = FiltersUtils.getFiltersPanel(col.dataFrame);
            const filter = panelFilters.getFilter(col.name);
            if(filter !== null) {

                panelFilters.removeFilter(col.name);
                return;
            }
        }

        const colGrid = cell.gridColumn;
        const col = this.adjustColumn(cell);//cell.tableColumn;
        const type = col.type;
        const typeSem = GridUtils.getSemType(col);
        if(typeSem instanceof SemType && DGApp.isVirtual(col)) {
            let crBack = col.getTag(AbstractRenderer.TAG_PRIMITIVE_COL_BACKGROUND);
            if(crBack === null || crBack === undefined) {
                crBack = AbstractRenderer.generateLightColor();
                col.setTag(AbstractRenderer.TAG_PRIMITIVE_COL_BACKGROUND, crBack);
            }

            const scrollH = grid.horzScroll;
            const scrollV = grid.vertScroll;
            const nXMin = scrollH.min;
            const nYMin = scrollV.min;

            let colGPrim = null;
            const arPrimCols = col.getTag(DGApp.PRIMITIVE_COLS_TAG_NAME);
            const nChildCount = arPrimCols.length;
            for(let nChild=0; nChild < nChildCount; ++nChild)
            {
                colGPrim = grid.columns.byName(arPrimCols[nChild]);
                colGPrim.visible = true;
                colGPrim.column.setTag(AbstractRenderer.TAG_PRIMITIVE_COL_BACKGROUND, crBack);
            }

            colGrid.visible = false;
            grid.scrollToPixels(nXMin, nYMin);

            return;
        }

        const strParentColName = col.getTag(DGApp.VIRTUAL_PARENT_TAG_NAME);
        if(strParentColName !== null && strParentColName !== undefined)
        {
            const scrollH = grid.horzScroll;
            const scrollV = grid.vertScroll;
            const nXMin = scrollH.min;
            const nYMin = scrollV.min;

                const colGridParent = grid.columns.byName(strParentColName);
                const arPrimCols = colGridParent.column.getTag(DGApp.PRIMITIVE_COLS_TAG_NAME);
                const nChildCount = arPrimCols.length;
                let colGPrim = null;
                for(var nChild=0; nChild < nChildCount; ++nChild)
                {
                    colGPrim = grid.columns.byName(arPrimCols[nChild]);
                    colGPrim.visible = false;
                }

                //if(colGridParent.column !== null && GridUtils.getSemType(colGridParent.column) instanceof CpdSemType)
                //{
                    //let nColCpd = MCBUtils.findGridColumnBySemType()
                   // grid.columns.byIndex(0).visible = true;
                //}
                //else
                colGridParent.visible = true;
            grid.scrollToPixels(nXMin, nYMin);
        }
    }




    /**
     * Exports the rendered content into html DOM tree. The default implementation
     * create a single canvas element and output the graphics content onto it. Subclasses
     * can override this method to create a DOM tree consisting of various html elements.
     * @param eParent a parent DOM node to which the output DOM tree will be appended.
     * @param cell the grid's cell to be rendererd.
     * @param nW the width o the drawing area.
     * @param nH the height o the drawing area.
     * @param crBack the textual CSS representation for the background color ("white", "red")
     * @returns {HTMLElement} a reference to the root element of the output DOM tree.
     */
    toHtml(eParent, cell, nW, nH, crBack)
    {
        const eCanvas = AbstractRenderer.toHtml(eParent, nW, nH, crBack)
        const g = eCanvas.getContext("2d");

        this.paint(g, cell, 0, 0, nW, nH);

        return eCanvas;
    }

    paint(g, cell, nX, nY, nW, nH, crBack)
    {
        if(crBack === null)
            return;

        let cr = null;
        cr = GridUtils.getCellBackgroundColor(cell);
        if(cr !== undefined && cr !== null && cr !== "#ffffff" && cr !== "#000000")
         crBack = cr;

        if(crBack === undefined)
            crBack = "white";

        g.fillStyle = crBack;
        g.fillRect(nX, nY, nW, nH);
    }
}
