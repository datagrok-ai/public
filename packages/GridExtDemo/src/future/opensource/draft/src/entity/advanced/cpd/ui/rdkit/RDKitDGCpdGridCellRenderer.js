import {ContainerGridCellRenderer} from "../../../../ui/ContainerGridCellRenderer";
import {MathUtils} from "../../../../../utils/MathUtils";
import {TextUtils} from "../../../../../utils/TextUtils";
import {DateUtils} from "../../../../../utils/DateUtils";
import {AbstractCpdEntityRenderer} from "../AbstractCpdEntityRenderer";

export class RDKitDGCpdGridCellRenderer extends ContainerGridCellRenderer {
    constructor(rendererEntity) {
        super(rendererEntity);
    }

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

    paint(g, cell, nX, nY, nW, nH, crBack) {
        const nRowGrid = cell.gridRow;
        const colCpd = cell.tableColumn;
        const nCol = cell.gridColumn.idx;
        const nRowTable = cell.tableRow.idx;

        const i = this.getInsets();
        const entity = colCpd.get(nRowTable);
        const cellStruct = cell.grid.cell("Structure", nRowGrid);


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

        g.font = cell.style.font;
        const nHStruct = AbstractCpdEntityRenderer.paintLabels(g, entity, [strConceptId, strLastTested], arBackColors, nX, nY, nW, nH, crBack)
        if(nHStruct > 10)
        {
            window.devicePixelRatio = 1;
            g.save();
            cellStruct.style.backColor = DG.Color.white;
            try{this.m_renderer.render(g, nX, nY, nW, nHStruct, cellStruct, cellStruct.style);}
            catch(e)
            {
             let rrr = 0;
             console.error("RDFIT DG Renderer " + nX + " " +  nY + " " +  nW + " " + nHStruc);
            }//_drawMolecule(nX, nY, nW, nH, g.canvas, entity.getSmiles(), "", false, false, false);//render(g, nX, nY, nW-5, nH-5, cell, cell.style);
            g.restore();
        }




       /*

        window.devicePixelRatio = 1;
        g.save();
        this.m_renderer.render(g, nX, nY, nW, nH, cell, cell.style);
        g.restore();*/
    }
}

RDKitDGCpdGridCellRenderer.create = async function()
{
    let rendererDGRDKit = null;
    try{rendererDGRDKit = await grok.functions.call("chem:rdkitCellRenderer");}
    catch(e)
    {
        return null;
    }
    const renderer = new RDKitDGCpdGridCellRenderer(rendererDGRDKit);

    return renderer;
}
