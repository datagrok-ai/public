import {GridCellRenderer} from "./GridCellRenderer";

export class ContainerGridCellRenderer extends GridCellRenderer
{
    constructor(renderer)
    {
        super();
        this.m_renderer = renderer;
    }

    getPreferredCellWidth()
    {
        return this.m_renderer.getPreferredCellWidth();
    }

    getPreferredCellHeight()
    {
        return this.m_renderer.getPreferredCellHeight();
    }

    paint(g, cell, nX, nY, nW, nH, crBack)
    {
        super.paint(g, cell, nX, nY, nW, nH, crBack);

        const nRow = cell.gridRow;
        const nCol = cell.gridColumn.idx;
        const font = cell.style.font;
        this.m_renderer.setFont(font);

        const entity = cell.cell.value;
        this.m_renderer.paint(g, entity, nX, nY, nW, nH, crBack);
    }
}
