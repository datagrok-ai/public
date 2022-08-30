import {IC50EntityRenderer} from "../IC50EntityRenderer";

export class DGIC50CanvasRenderer extends DG.GridCellRenderer
{
    constructor()
    {
        super();
        this.m_renderer = new IC50EntityRenderer();
    }

    render(g, nX, nY, nW, nH, entity, context)
    {
        try{this.m_renderer.paint(g, ctx, nX, nY, nW, nH, "white");}
        catch(e)
        {
            let ccc = cell.cell.value;
            throw e;
        }

        //g.fillText(fruit.name, x + 10, y + 10);
    }
}