import {ContainerGridCellRenderer} from "../../ui/ContainerGridCellRenderer";
import {TextDateEntityRenderer} from "./TextDateEntityRenderer";
import {GridUtils} from "../../../utils/GridUtils";
import {ErrorUtils} from "../../../utils/ErrorUtils";
import {TextDateSemType} from "../TextDateSemType";
import {DGApp} from "../../app/DGApp";
import {SemType} from "../../SemType";
import * as DG from "datagrok-api/dg";

export class TextDateGridCellRenderer extends ContainerGridCellRenderer
{
    constructor()
    {
        super(new TextDateEntityRenderer());
    }


    paint(g, cell, nX, nY, nW, nH, crBack)
    {
        const col = cell.tableColumn;
        const dframe = col.dataFrame;
        const typeSem = GridUtils.getSemType(col);
        ErrorUtils.verifyClass(typeSem,TextDateSemType);
        const arColNames = col.getTag(DGApp.PRIMITIVE_COLS_TAG_NAME);
        const arColsGroup = GridUtils.names2Cols(dframe, arColNames);

        const arTypeKeys = TextDateSemType.CDF_TYPES;
        const arColKeys = SemType.col2CDFPrimTypes(arColsGroup);
        const arHitColIdxs = new Array(arTypeKeys.length);
        const nHitCount = SemType.fillCompatibilityLevelHits(arHitColIdxs, arTypeKeys, arColKeys);
        const colText = arColsGroup[arHitColIdxs[0]]; //0 must exist
        const cellGridTextPrim = cell.grid.cell(colText.name, cell.gridRow);
        const cr = GridUtils.getCellBackgroundColor(cellGridTextPrim);
        if(cr !== undefined && cr !== null && cr !== "#ffffff")
        {
            crBack = cr;
        }

        super.paint(g, cell, nX, nY, nW, nH, crBack);
    }
}
