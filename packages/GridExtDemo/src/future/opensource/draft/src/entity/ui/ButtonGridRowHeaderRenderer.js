import {GridRowHeaderRenderer} from "./GridRowHeaderRenderer";
import {ButtonGridColumnHeaderRenderer} from "./ButtonGridColumnHeaderRenderer";

export class ButtonGridRowHeaderRenderer extends GridRowHeaderRenderer
{

    getInsets()
    {
        const i = ButtonGridColumnHeaderRenderer.BORDER_INSETS.clone();
        i.m_nR -= 10;

        return i;
    }

    paintBackground(g, nX, nY, nW, nH)
    {
        ButtonGridColumnHeaderRenderer.paintBackground(g, nX, nY, nW, nH);
    }

    paintBorder(g, nX, nY, nW, nH)
    {
        ButtonGridColumnHeaderRenderer.paintBorder(g, nX, nY, nW, nH, ButtonGridColumnHeaderRenderer.GAP, "white");
    }

    paintFilter(g, nX, nY, nW, nH)
    {
        ButtonGridColumnHeaderRenderer.paintFilter(g, nX, nY, nW, nH, this.getInsets());
    }

    paint(g, nRowGrid, nRowTable, grid, nX, nY, nW, nH, crBack)
    {
        super.paint(g, nRowGrid, nRowTable, grid, nX, nY, nW, nH, crBack);
    }
}
