import {IC50SemType} from "../../IC50SemType";
import {IC50EntityRenderer} from "../IC50EntityRenderer";

export class DGIC50GridCellRenderer extends DG.GridCellRenderer
{
    constructor()
    {
        super();
        this.m_renderer = new IC50EntityRenderer();
    }

    get cellType() { return IC50SemType.Instance; }
    render(g, nX, nY, nW, nH, cell, style)
    {

        if(cell.cell.value === null)
            return;

        //let ctx = new SemEntityRendererContext(cell.cell, cell.cell.column, cell.cell.column.dataFrame);

        try{this.m_renderer.paint(g, ctx, nX, nY, nW, nH, "white");}
        catch(e)
        {
            let ccc = cell.cell.value;
            throw e;
        }

        //g.fillStyle = cell.cell.value.color;
        //g.fillText(cell.cell.value.name, x + 10, y + 10);
    }
}
