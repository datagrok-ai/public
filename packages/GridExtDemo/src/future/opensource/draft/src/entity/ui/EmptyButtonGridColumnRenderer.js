import {ButtonGridColumnHeaderRenderer} from "./ButtonGridColumnHeaderRenderer";

export class EmptyButtonGridColumnRenderer extends ButtonGridColumnHeaderRenderer
{
    adjustColumnName(strColName)
    {
        return "";
    }

    paintFilter(g, nX, nY, nW, nH)
    {
        const c = 0;
    }
}