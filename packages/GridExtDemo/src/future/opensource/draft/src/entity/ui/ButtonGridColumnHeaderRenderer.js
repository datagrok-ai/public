import {GridColumnHeaderRenderer} from "./GridColumnHeaderRenderer";
import {FlowTextRenderer} from "./FlowTextRenderer";
import {Insets} from "../../geom/Insets";
import {SemType} from "../SemType";
import {GridUtils, SortOptionCtx} from "../../utils/GridUtils";
import * as ui from 'datagrok-api/ui';
import {FiltersUtils} from "../../ui/filters/FiltersUtils";
import {DGApp} from "../app/DGApp";
import {TextUtils} from "../../utils/TextUtils";

export class ButtonGridColumnHeaderRenderer extends GridColumnHeaderRenderer
{
    constructor()
    {
        super();
        this.m_renderer = new FlowTextRenderer();

        this.m_bFilterEnabled = true;
        this.m_nColLastPressed = -1;
        this.m_nXLastPressed = -1;
        this.m_nYLastPressed = -1;
        this.m_nButtonLastPressed = -1;

        this.m_bPaintLabel = true;
    }

    getPaintLabel() {return this.m_bPaintLabel;}
    setPaintLabel(bPaint)
    {
      this.m_bPaintLabel = bPaint;
    }

    toPopupItemsArray(cell)
    {
        const arItems = ["Filter"];

        const viewerFilter = FiltersUtils.getFiltersViewer();
        if(viewerFilter !== null) {
            const col = this.adjustColumn(cell);
            const panelFilters = FiltersUtils.getFiltersPanel(col.dataFrame);
            const filter = panelFilters.getFilter(col.name);
            if (filter !== null)
                arItems.push("Close Filter");
        }
        return arItems;
    }

    onPopupAction(cell, nItem)
    {
        if(nItem === 0)
         FiltersUtils.floatFilterToColumn(this, cell);
        else
        {
            const col = this.adjustColumn(cell);
            const panelFilters = FiltersUtils.getFiltersPanel(col.dataFrame);
            panelFilters.removeFilter(col.name);
        }
    }


    getInsets() {return ButtonGridColumnHeaderRenderer.BORDER_INSETS;}
    isFilterEnabled() {return this.m_bFilterEnabled;}
    setFilterEnabled(bEnable)
    {
        this.m_bFilterEnabled = bEnable;
    }

    toHtml(eParent, column, nW, nH, crBack)
    {
        let eCanvas = ui.canvas();
        eCanvas.width  = nW;
        eCanvas.height = nH;
        eParent.appendChild(eCanvas);

        var g = eCanvas.getContext("2d");
        if(crBack === undefined)
            crBack = 'white';

        g.fillStyle = crBack;
        g.fillRect(0, 0, nW, nH);

        this.paint(g, column, dframe, 0, 0, nW, nH);

        return eCanvas;
    }


    paintBackground(g, nX, nY, nW, nH, crBack)
    {
        ButtonGridColumnHeaderRenderer.paintBackground(g, nX, nY, nW, nH, crBack);
    }

    paintBorder(g, nX, nY, nW, nH)
    {
        ButtonGridColumnHeaderRenderer.paintBorder(g, nX, nY, nW, nH, ButtonGridColumnHeaderRenderer.GAP);
    }

    paintFilter(g, nX, nY, nW, nH)
    {
        ButtonGridColumnHeaderRenderer.paintFilter(g, nX, nY, nW, nH, this.getInsets());
    }

    paintColLabel(g, strLabel, nX, nY, nW, nH)
    {
        const i = this.getInsets();

        const nWAvail = nW - ButtonGridColumnHeaderRenderer.FILTER_WIDTH - i.getL() - i.getR();
        strLabel = TextUtils.trimText(strLabel, g, nW);

        g.fillStyle = "blue";
        g.textBaseline = "bottom";
        g.fillText(strLabel, nX + i.getL(), nY + nH - i.getB());
        //ButtonGridColumnHeaderRenderer.paintFilter(g, nX, nY, nW, nH, this.getInsets());
    }


    paint(g, column, nX, nY, nW, nH, crBack)
    {
     super.paint(g, column, nX, nY, nW, nH, crBack);

     if(this.m_bFilterEnabled)
        this.paintFilter(g, nX, nY, nW, nH);


        const font = this.getFont();
        this.m_renderer.setFont(font);

        g.font = this.m_renderer.getFont();
        const tm = g.measureText("W");
        const nGapLines = this.m_renderer.getLinesGap();
        const nHFont = Math.abs(tm.actualBoundingBoxAscent) + tm.actualBoundingBoxDescent + nGapLines;
        const bVertText = nW < nHFont * 4;


     if(!bVertText && this.m_bPaintLabel) {
         const arColPri = column === null || column.column === null ? null : column.column.getTag(DGApp.PRIMITIVE_COLS_TAG_NAME);
         if (arColPri !== undefined && arColPri !== null && arColPri.length > 0) {

             //const font = this.m_renderer.getFont();
             const nHFont = TextUtils.getFontSize(font);
             const nHFontSub = nHFont < 0 ? nHFont : nHFont - 2;
             const fontSub = TextUtils.setFontSize(font, nHFontSub);

             g.font = fontSub;
             this.paintColLabel(g, arColPri.length.toString(), nX, nY, nW, nH);
         }

         const strLabel = column === null || column.column === null ? null : column.column.getTag(DGApp.PRIMITIVE_COL_LABEL_TAG_NAME);
         if (strLabel !== undefined && strLabel !== null) {

             //const font = this.m_renderer.getFont();
             const nHFont = TextUtils.getFontSize(font);
             const nHFontSub = nHFont < 0 ? nHFont : nHFont - 2;
             const fontSub = TextUtils.setFontSize(font, nHFontSub);

             g.font = fontSub;
             this.paintColLabel(g, strLabel, nX, nY, nW, nH);
         }
     }

     const i = this.getInsets();

     let strColName = this.adjustColumnName(column === null ? "" : column.name);

     if(DGApp.USE_ALT_COL_NAMES) {
         const strAltColName = column === null || column.column === null ? null : column.column.getTag(DGApp.VIRTUAL_COL_ALT_NAME_TAG_NAME);
         if(strAltColName !== null && strAltColName !== undefined)
             strColName = strAltColName;
     }

     nW -= (i.getL()+i.getR());
     nH -= (i.getT()+i.getB())

     //const font = this.getFont();
     this.m_renderer.setFont(font);

     /*
     g.font = this.m_renderer.getFont();
     const tm = g.measureText("W");
     const nGapLines = this.m_renderer.getLinesGap();
     const nHFont = Math.abs(tm.actualBoundingBoxAscent) + tm.actualBoundingBoxDescent + nGapLines;
     const bVertText = nW < nHFont * 4;*/
     this.m_renderer.setVerticalTextMode(bVertText);

     this.m_renderer.paint(g, strColName, nX+i.getL(), nY+i.getT(), nW, nH, crBack);

     return true;
    }

    adjustColumnName(strColName)
    {
        return strColName;
    }


    adjustColumn(cell)
    {
        return cell.tableColumn;
    }

    hitTestFilter(cell, nX, nY)
    {
      const rc = cell.bounds;
      if(nX < 0 + rc.width -1 - ButtonGridColumnHeaderRenderer.GAP - ButtonGridColumnHeaderRenderer.FILTER_WIDTH - ButtonGridColumnHeaderRenderer.FILTER_GAP ||
         nX > 0 + rc.width -1 - ButtonGridColumnHeaderRenderer.GAP - ButtonGridColumnHeaderRenderer.FILTER_GAP)
        return false;

      const nH = cell.grid.getOptions().look.colHeaderHeight;
      rc.height = nH;

        if(nY < rc.height -10- 1- ButtonGridColumnHeaderRenderer.GAP - ButtonGridColumnHeaderRenderer.FILTER_HEIGHT - ButtonGridColumnHeaderRenderer.FILTER_GAP  ||
           nY > rc.height - 1- ButtonGridColumnHeaderRenderer.GAP - ButtonGridColumnHeaderRenderer.FILTER_GAP)
         return false;

        return true;
    }


    onMouseMoved(cell, nX, nY, nButton)
    {
        if(!this.isFilterEnabled())
            return;

        const eCanvas = cell.grid.canvas;
        const g = eCanvas.getContext("2d");
        const rc = cell.bounds;
        const nH = cell.grid.getOptions().look.colHeaderHeight;
        rc.height = nH;
        rc.y = 0;
        const b = this.hitTestFilter(cell, nX, nY);
        ButtonGridColumnHeaderRenderer.paintFilter(g, rc.x, rc.y, rc.width, rc.height, b ? "orange" : "white");
    }


    onMousePressed(cell, nX, nY, nButton)
    {
      const col = cell.gridColumn;
      const strName = col.name;
      const nColGrid = col.idx;

      this.m_nColLastPressed = nColGrid;
      this.m_nXLastPressed = nX;
      this.m_nYLastPressed = nY;
      this.m_nButtonLastPressed = nButton;

      const b = this.hitTestFilter(cell, nX, nY);
      //console.log("Pressed: " + nColGrid);
   }


    async onMouseClicked(cell, nX, nY, nButton)
    {
        const colGrid = cell.gridColumn;
        const strName = colGrid.name;
        const nColGrid = colGrid.idx;

        if(this.m_nColLastPressed !== nColGrid || this.m_nXLastPressed !== nX ||
           this.m_nYLastPressed !== nY || this.m_nButtonLastPressed !== nButton || nButton !== 0)
        {
            this.m_nColLastPressed = -1;
            this.m_nXLastPressed = -1;
            this.m_nYLastPressed = -1;
            this.m_nButtonLastPressed = -1;
            console.log("Clicked Exit: " + nColGrid);
            return;
        }

        const bOld = this.hitTestFilter(cell, this.m_nXLastPressed, this.m_nYLastPressed);
        const b = this.hitTestFilter(cell, nX, nY);
        if(bOld && b)
        {
            this.onFilterIconClicked(cell);
        }
        else
        {
            const colT = this.adjustColumn(cell);
            const typeSem = GridUtils.getSemType(colT);
            if(typeSem instanceof SemType)
            {
                const grid = cell.grid;
                const nRecordCount = colT.length;
                const arData = new Array(nRecordCount);
                const nIdxSort = typeSem.getDefaultSortOptionIndex();


                let bAscend = true;
                if(grid.dart.m_ctxSort !== null && grid.dart.m_ctxSort.m_nIdxColGrid === nColGrid &&  grid.dart.m_ctxSort.m_nSortOption === nIdxSort)
                    bAscend = !grid.dart.m_ctxSort.m_bSortAscend;

                await SemType.sort(arData, colT, -1, nIdxSort, bAscend);
                GridUtils.setRowOrder(grid, arData);  //temporary fixed the DG bug
                //temporary work around until DG provides support
                grid.dart.m_ctxSort = new SortOptionCtx(nColGrid, nIdxSort, bAscend);
            }

            //console.log("Clicked: " + nColGrid);
        }
    }

    onFilterIconClicked(cell)
    {
        FiltersUtils.floatFilterToColumn(this, cell);
    }
}

ButtonGridColumnHeaderRenderer.BORDER_RADIUS = 5;
ButtonGridColumnHeaderRenderer.GAP = 0;     //2
ButtonGridColumnHeaderRenderer.PADDING = 2;//3
ButtonGridColumnHeaderRenderer.BORDER_INSETS = new Insets(ButtonGridColumnHeaderRenderer.GAP+ButtonGridColumnHeaderRenderer.PADDING,
    ButtonGridColumnHeaderRenderer.GAP+ButtonGridColumnHeaderRenderer.PADDING,
    ButtonGridColumnHeaderRenderer.GAP+ButtonGridColumnHeaderRenderer.PADDING,
    ButtonGridColumnHeaderRenderer.GAP+ButtonGridColumnHeaderRenderer.PADDING+10);


ButtonGridColumnHeaderRenderer.FILTER_WIDTH = 10;
ButtonGridColumnHeaderRenderer.FILTER_HEIGHT= 12;
ButtonGridColumnHeaderRenderer.FILTER_GAP= 1;
ButtonGridColumnHeaderRenderer.FILTER_X_COORDS = [0,9,6,6,3,3,0];
ButtonGridColumnHeaderRenderer.FILTER_Y_COORDS = [0,0,5,11,11,5,0];


ButtonGridColumnHeaderRenderer.paintBackground = function(g, nX, nY, nW, nH, crBack)
{
    //const crBtn = "rgb(214, 217, 223)";

    const crBtn = crBack !== undefined ? crBack : "rgb(222,221,220)";

    const gradient = g.createLinearGradient(nX,nY, nX,nY+nH-1);
    gradient.addColorStop(0, 'white');
    gradient.addColorStop(.5, crBtn);
    gradient.addColorStop(.7, crBtn);
    gradient.addColorStop(1, 'white');
    g.fillStyle = gradient;
    g.fillRect(nX+1, nY+1, nW-2, nH-2);
}


ButtonGridColumnHeaderRenderer.paintBorder = function(g, nX, nY, nW, nH, nGap)
{
    if(nGap === undefined)
        nGap = ButtonGridColumnHeaderRenderer.GAP;
       const nRadius = ButtonGridColumnHeaderRenderer.BORDER_RADIUS;

    const crBorder = "rgba(105,109, 114)";
    g.beginPath();
    g.lineWidth = 1;
    g.strokeStyle = crBorder;
    g.arc(nX + nGap + nRadius, nY + nGap + nRadius, nRadius, Math.PI, 1.5*Math.PI);
    g.stroke();

    g.beginPath();
    g.arc(nX + nW -1 -nGap - nRadius, nY + nGap + nRadius, nRadius, 1.5*Math.PI, 2.0*Math.PI);
    g.stroke();

    g.beginPath();
    g.arc(nX + nW -1- nGap - nRadius, nY + nH - 1 - nRadius -nGap, nRadius, 0, Math.PI/2);
    g.stroke();

    g.beginPath();
    g.arc(nX + nGap + nRadius, nY + nH - 1 - nRadius -nGap, nRadius, Math.PI/2.0, Math.PI);
    g.stroke();


    //Draw Lines
    g.beginPath();        //Horz Top
    g.moveTo(nX + nGap + nRadius, nY + nGap);
    g.lineTo(nX + nW - 1 - nGap - nRadius, nY + nGap);
    g.stroke();

    g.beginPath(); //Vert Left
    g.moveTo(nX + nGap, nY+nRadius);
    g.lineTo(nX + nGap, nY+ nH -1 - nGap -nRadius);
    g.stroke();

    g.beginPath();    //Vert Right
    g.moveTo(nX +nW -1 - nGap, nY+nRadius);
    g.lineTo(nX +nW -1 - nGap, nY+ nH - 1-nRadius);
    g.stroke();

    g.beginPath(); //Horz Bottom
    g.moveTo(nX +nGap + nRadius, nY + nH-1-nGap);
    g.lineTo(nX + nW-1 - nGap -nRadius, nY + nH-1-nGap);
    g.stroke();

        /*
               g.strokeStyle = "red";
    g.beginPath();
    g.moveTo(nX, nY);
    g.lineTo(nX, nY+ nH -1);
    g.stroke();

    g.strokeStyle = "blue";
    g.beginPath();    //Vert Right
    g.moveTo(nX +nW -1, nY);
    g.lineTo(nX +nW -1, nY+ nH - 1);
    g.stroke();
          */

}




ButtonGridColumnHeaderRenderer.paintFilter = function(g, nX, nY, nW, nH, crFill)
{
    const nWIcon = ButtonGridColumnHeaderRenderer.FILTER_WIDTH;
    const nHIcon = ButtonGridColumnHeaderRenderer.FILTER_HEIGHT;
    const nGapFilter = ButtonGridColumnHeaderRenderer.FILTER_GAP;

    const nWAvail = nW - 2*(ButtonGridColumnHeaderRenderer.GAP + nGapFilter);
    if(nWAvail < nWIcon)
        return;

    const x = nX + nW  - 1 - nWIcon -(ButtonGridColumnHeaderRenderer.GAP+nGapFilter);
    const y = nY + nH  - 1 - nHIcon -(ButtonGridColumnHeaderRenderer.GAP+nGapFilter);

    g.strokeStyle = "rgba(105,109, 114)";
    g.fillStyle = crFill === undefined ? 'white' : crFill;

    const arXs = ButtonGridColumnHeaderRenderer.FILTER_X_COORDS;
    const arYs = ButtonGridColumnHeaderRenderer.FILTER_Y_COORDS;

    g.translate(x, y);
    g.beginPath();
    g.moveTo(arXs[0],arXs[0]);

    const nLength = arXs.length;
    for(var n=1; n<nLength; ++n)
    {
        g.lineTo(arXs[n],arYs[n]);
    }

    g.fill();
    g.stroke();
    g.translate(-x, -y);
}

