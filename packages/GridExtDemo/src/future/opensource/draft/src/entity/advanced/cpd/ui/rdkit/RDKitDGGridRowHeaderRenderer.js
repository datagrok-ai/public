import {CpdRowHeaderRenderer} from "../CpdRowHeaderRenderer";
import {AbstractCpdEntityRenderer} from "../AbstractCpdEntityRenderer";
import {GridUtils} from "../../../../../utils/GridUtils";
import {CpdSemType} from "../../CpdSemType";
import {MathUtils} from "../../../../../utils/MathUtils";
import {TextUtils} from "../../../../../utils/TextUtils";
import {DateUtils} from "../../../../../utils/DateUtils";

export class RDKitDGGridRowHeaderRenderer extends CpdRowHeaderRenderer
{
    constructor(rendererRDGDKit)
    {
        super();
        this.m_rendererDGRDKit = rendererRDGDKit;
    }

    paintContent(g, nRowGrid, nRowTable, grid, nX, nY, nW, nH, crBack) {}

    paint(g, nRowGrid, nRowTable, grid, nX, nY, nW, nH, crBack)
    {
        super.paint(g, nRowGrid, nRowTable, grid, nX, nY, nW, nH, crBack);

        const dframe = grid.dataFrame;
        const nColCpd = GridUtils.findColumnBySemType(dframe, CpdSemType);
        const colCpd = dframe.columns.byIndex(nColCpd);
        //const colStruct = dframe.columns.byName("Structure");
        const cell = grid.cell("Structure", nRowGrid);

        const i = this.getInsets();
        const entity = colCpd.get(nRowTable);
        //this.m_rendererEntity.paint(g, entity, nX + i.getR(), nY+i.getT(), nW-(i.getL()+i.getR()), nH-(i.getT() + i.getB()), crBack);

        if (crBack === undefined)
            crBack = "white";

        g.fillStyle = crBack;

        nX += i.getR();
        nY += i.getT();
        nW -= (i.getL()+i.getR());
        nH -= (i.getT() + i.getB());

        g.fillRect(nX, nY, nW, nH);

        let strConceptId = entity.getConceptId();
        if(strConceptId === undefined)
            strConceptId = "";

        const nTime = entity.getLastTested();
        const strLastTested = MathUtils.isNullValue(nTime) || nTime < 0 ? "" : TextUtils.formatDate(nTime);
        const arBackColors = [null, DateUtils.getRecentDateColor(nTime)];
        const nHStruct = AbstractCpdEntityRenderer.paintLabels(g, entity, [strConceptId, strLastTested], arBackColors, nX, nY, nW, nH, crBack)
        if(nHStruct > 10)
        {
            window.devicePixelRatio = 1;
            g.save();
            this.m_rendererDGRDKit.render(g, nX, nY, nW, nHStruct, cell, cell.style);//_drawMolecule(nX, nY, nW, nH, g.canvas, entity.getSmiles(), "", false, false, false);//render(g, nX, nY, nW-5, nH-5, cell, cell.style);
            g.restore();
        }
    }


}

RDKitDGGridRowHeaderRenderer.create = async function()
{
    let rendererDGRDKit = null;
    try{rendererDGRDKit = await grok.functions.call("chem:rdkitCellRenderer");}
    catch(e)
    {
        return null;
    }

    const renderer = new RDKitDGGridRowHeaderRenderer(rendererDGRDKit);
    return renderer;
}
