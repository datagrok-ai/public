import {CpdSemType} from "../../CpdSemType";
import {GridUtils} from "../../../../../utils/GridUtils";
import {RDKitCpdEntityRenderer} from "./RDKitCpdEntityRenderer";
import {CpdRowHeaderRenderer} from "../CpdRowHeaderRenderer";

export class RDKitCpdGridRowHeaderRenderer extends CpdRowHeaderRenderer
{
    constructor(rendererEntity)
    {
        super();
        this.m_rendererEntity = rendererEntity;
    }

    paintContent(g, nRowGrid, nRowTable, grid, nX, nY, nW, nH, crBack)
    {

    }


    paint(g, nRowGrid, nRowTable, grid, nX, nY, nW, nH, crBack)
    {
        super.paint(g, nRowGrid, nRowTable, grid, nX, nY, nW, nH, crBack);

        const dframe = grid.dataFrame;
        const nColCpd = GridUtils.findColumnBySemType(dframe, CpdSemType);
        const colCpd = dframe.columns.byIndex(nColCpd);

        const i = this.getInsets();
        const entity = colCpd.get(nRowTable);
        this.m_rendererEntity.paint(g, entity, nX + i.getR(), nY+i.getT(), nW-(i.getL()+i.getR()), nH-(i.getT() + i.getB()), crBack);
    }
}

RDKitCpdGridRowHeaderRenderer.create = async function()
{
    const rendererEntity = await RDKitCpdEntityRenderer.create();
    const renderer = new RDKitCpdGridRowHeaderRenderer(rendererEntity);

    return renderer;
}
