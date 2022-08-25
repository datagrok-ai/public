import {PrimitiveGridCellRenderer} from "../../../primitive/ui/PrimitiveGridCellRenderer";
import {GridUtils} from "../../../../utils/GridUtils";
import {SemType} from "../../../SemType";
import {DGApp} from "../../../app/DGApp";

export class CpdPrimitiveGridCellRenderer extends PrimitiveGridCellRenderer
{
    onPopupAction(cell, nItem)
    {
        const grid = cell.grid;
        const colGrid = cell.gridColumn;

        const col = cell.tableColumn;
        const typeSem = GridUtils.getSemType(col);
        if(typeSem instanceof SemType && DGApp.isVirtual(col)) {
            let colGPrim = null;
            const arPrimCols = col.getTag(DGApp.PRIMITIVE_COLS_TAG_NAME);
            const nChildCount = arPrimCols.length;
            for(let nChild=0; nChild < nChildCount; ++nChild)
            {
                colGPrim = grid.columns.byName(arPrimCols[nChild]);
                colGPrim.visible = true;
            }

            colGrid.visible = false;
            return;
        }

        const strParentColName = col.getTag(DGApp.VIRTUAL_PARENT_TAG_NAME);
        if(strParentColName !== null && strParentColName !== undefined)
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

            //grid.columns.byIndex(0).visible = true;//0 - rows header

        }
    }
}
