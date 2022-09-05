import {ButtonGridRowHeaderRenderer} from "../../../ui/ButtonGridRowHeaderRenderer";
import {GridUtils, SortOptionCtx} from "../../../../utils/GridUtils";
import {CpdSemType} from "../CpdSemType";
import {SemType} from "../../../SemType";
import {DGApp} from "../../../app/DGApp";
import {AbstractRenderer} from "../../../ui/AbstractRenderer";
import {FiltersUtils} from "../../../../ui/filters/FiltersUtils";

export class CpdRowHeaderRenderer extends ButtonGridRowHeaderRenderer
{
    getPreferredCellWidth()
    {
        const i = this.getInsets();
        return 160 + i.getL() + i.getR();
    }

    getPreferredCellHeight()
    {
        const i = this.getInsets();
        return 95 + i.getT() + i.getB();
    }
          /*
    getPopupItemCount(cell)
    {
        const dframe = cell.grid.dataFrame;
        const nColCpd = GridUtils.findColumnBySemType(dframe, CpdSemType);

        const col = dframe.columns.byIndex(nColCpd);
        const typeSem = GridUtils.getSemType(col);
        const nSortCount =  typeSem.getSortOptionCount();

        if(DGApp.isVirtual(col))
            return nSortCount + 1;

        const strParentColName = col.getTag(DGApp.VIRTUAL_PARENT_TAG_NAME);
        if(strParentColName !== null && strParentColName !== undefined)
            return nSortCount + 1;

        return nSortCount +1;//filter
    }       */


    toPopupItemsArray(cell)
    {
        const dframe = cell.grid.dataFrame;
        const nColCpd = GridUtils.findColumnBySemType(dframe, CpdSemType);
        const col = dframe.columns.byIndex(nColCpd);
        const typeSem = GridUtils.getSemType(col);
        const arSortOptions = SemType.toSortOptionNamesArray(typeSem);

        let str = null;
        if(DGApp.hasPrimitiveColumns(col))
            str = "Expand to Primitives";
        else
        {
            const strParentColName = col.getTag(DGApp.VIRTUAL_PARENT_TAG_NAME);
            if (strParentColName !== null && strParentColName !== undefined)
                str = "Collapse to Object";
        }

        if(str !== null)
            arSortOptions.push(str);

        arSortOptions.push("Filter");

        return arSortOptions;
    }


    async onPopupAction(cell, nItem)
    {
        const grid = cell.grid;
        const nOption = nItem;//SemType.sortOptionIndexFromName(typeSem, strItem);
        const dframe = grid.dataFrame;
        const colGrid = cell.gridColumn;
        const nColCpd = GridUtils.findColumnBySemType(dframe, CpdSemType);
        const col = dframe.columns.byIndex(nColCpd);
        const typeSem = GridUtils.getSemType(col);

        const nSortCount =  typeSem.getSortOptionCount();
        //the sort part
        if(nItem < nSortCount) {
            const nRecordCount = col.length;
            const arData = new Array(nRecordCount);

            let bAscend = true;
            const ctxSort = GridUtils.getSortOptionContext(grid);
            if (typeSem.supportsSortDirectionSwitch(nOption) && ctxSort !== null && ctxSort.m_nIdxColGrid === colGrid.idx && ctxSort.m_nSortOption === nOption)
                bAscend = !ctxSort.m_bSortAscend;

            const nRecord = cell.tableRowIndex;
            await SemType.sort(arData, col, nRecord, nOption, bAscend);
            GridUtils.setRowOrder(grid, arData);
            //grid.setRowOrder(arData);
            grid.scrollToCell(0, 0);

            GridUtils.setSortOptionContext(grid, new SortOptionCtx(colGrid.idx, nOption, bAscend));
            return;
        }

        //Expand/Collapse
        if(nItem === nSortCount && DGApp.hasPrimitiveColumns(col))
        {
            let crBack = col.getTag(AbstractRenderer.TAG_PRIMITIVE_COL_BACKGROUND);
            if(crBack === null || crBack === undefined) {
                crBack = AbstractRenderer.generateLightColor();
                col.setTag(AbstractRenderer.TAG_PRIMITIVE_COL_BACKGROUND, crBack);
            }

            let colGPrim = null;
            const arPrimCols = col.getTag(DGApp.PRIMITIVE_COLS_TAG_NAME);
            const nChildCount = arPrimCols.length;
            for(var nChild=0; nChild < nChildCount; ++nChild)
            {
                colGPrim = grid.columns.byName(arPrimCols[nChild]);
                colGPrim.visible = true;
                colGPrim.column.setTag(AbstractRenderer.TAG_PRIMITIVE_COL_BACKGROUND, crBack);
            }

            //from table header colGrid.visible = false;
            return;
        }

        const strParentColName = col.getTag(DGApp.VIRTUAL_PARENT_TAG_NAME);
        if(nItem === nSortCount && strParentColName !== null && strParentColName !== undefined)
        {
            const colGridParent = grid.columns.byName(strParentColName);
            const arPrimCols = colGridParent.column.getTag(DGApp.PRIMITIVE_COLS_TAG_NAME);
            const nChildCount = arPrimCols.length;
            let colGPrim = null;
            for(let nChild=0; nChild < nChildCount; ++nChild)
            {
                colGPrim = grid.columns.byName(arPrimCols[nChild]);
                colGPrim.visible = false;
            }

            if(colGridParent.column !== null && GridUtils.getSemType(colGridParent.column) instanceof CpdSemType)
            {
                //let nColCpd = MCBUtils.findGridColumnBySemType()
                //gridd.columns.byIndex(0).visible = true;
            }
            else
                colGridParent.visible = true;
        }


        //FiltersUtils.getFilterToColumn(this, cell);
    }
}
