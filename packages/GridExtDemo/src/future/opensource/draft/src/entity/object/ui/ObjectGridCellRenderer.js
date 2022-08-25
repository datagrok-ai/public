import {ContainerGridCellRenderer} from "../../ui/ContainerGridCellRenderer";
import {ObjectEntityRenderer} from "./ObjectEntityRenderer";
import {GridUtils} from "../../../utils/GridUtils";
import {ErrorUtils} from "../../../utils/ErrorUtils";
import {ObjectSemType} from "../ObjectSemType";
import {DGApp} from "../../app/DGApp";
import * as DG from "datagrok-api/dg";

export class ObjectGridCellRenderer extends ContainerGridCellRenderer
{
    constructor()
    {
        super(new ObjectEntityRenderer());
    }


    paint(g, cell, nX, nY, nW, nH, crBack)
    {
        const col = cell.tableColumn;
        const dframe = col.dataFrame;
        const typeSem = GridUtils.getSemType(col);
        ErrorUtils.verifyClass(typeSem,ObjectSemType);
        const arColNames = col.getTag(DGApp.PRIMITIVE_COLS_TAG_NAME);
        const arColsGroup = GridUtils.names2Cols(dframe, arColNames);

        const colObj = arColsGroup[0]; //0 must exist
        const cellGridObjPrim = cell.grid.cell(colObj.name, cell.gridRow);
        const cr = GridUtils.getCellBackgroundColor(cellGridObjPrim);
        if(cr !== undefined && cr !== null && cr !== "#ffffff")
        {
            crBack = cr;
        }

        let s = cellGridObjPrim.style;

        super.paint(g, cell, nX, nY, nW, nH, crBack);
    }
}
