import {GridRowHeaderRenderer} from "../../../../ui/GridRowHeaderRenderer";
import {GridUtils} from "../../../../../utils/GridUtils";
import {CpdSemType} from "../../CpdSemType";
import {OCLCpdEntityRenderer} from "./OCLCpdEntityRenderer";

export class OCLCpdGridRowHeaderRenderer extends GridRowHeaderRenderer
{
    constructor()
    {
        super();
        this.m_rendererEntity = new OCLCpdEntityRenderer();
    }

    paint(g, nRowGrid, nRowTable, grid, nX, nY, nW, nH, crBack) {
        const dframe = grid.dataFrame;
        const nColCpd = GridUtils.findColumnBySemType(dframe, CpdSemType);
        const colCpd = dframe.columns.byIndex(nColCpd);

        const entity = colCpd.get(nRowTable);
        this.m_rendererEntity.paint(g, entity, nX, nY, nW, nH, crBack);
    }
}