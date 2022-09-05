import {ContainerGridCellRenderer} from "../../ui/ContainerGridCellRenderer";
import {QNumDateEntityRenderer} from "./QNumDateEntityRenderer";
import {GridUtils} from "../../../utils/GridUtils";
import {ErrorUtils} from "../../../utils/ErrorUtils";
import {DGApp} from "../../app/DGApp";
import {SemType} from "../../SemType";
import * as DG from "datagrok-api/dg";
import {QNumDateSemType} from "../QNumDateSemType";

export class QNumDateGridCellRenderer extends ContainerGridCellRenderer
{
    constructor()
    {
        super(new QNumDateEntityRenderer());
    }


    paint(g, cell, nX, nY, nW, nH, crBack)
    {
        const col = cell.tableColumn;
        const dframe = col.dataFrame;
        const typeSem = GridUtils.getSemType(col);
        ErrorUtils.verifyClass(typeSem,QNumDateSemType);
        const arColNames = col.getTag(DGApp.PRIMITIVE_COLS_TAG_NAME);
        const arColsGroup = GridUtils.names2Cols(dframe, arColNames);

        const arTypeKeys = QNumDateSemType.CDF_TYPES;
        const arColKeys = SemType.col2CDFPrimTypes(arColsGroup);
        const arHitColIdxs = new Array(arTypeKeys.length);
        const nHitCount = SemType.fillCompatibilityLevelHits(arHitColIdxs, arTypeKeys, arColKeys);
        const colQNumy = arColsGroup[arHitColIdxs[0]]; //0 must exist
        const cellGridQNumPrim = cell.grid.cell(colQNumy.name, cell.gridRow);
        const cr = GridUtils.getCellBackgroundColor(cellGridQNumPrim);
        if(cr !== undefined && cr !== null && cr !== "#ffffff")
        {
            crBack = cr;
        }

        super.paint(g, cell, nX, nY, nW, nH, crBack);
    }
}
