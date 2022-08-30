import {EmptyButtonGridColumnRenderer} from "../../../ui/EmptyButtonGridColumnRenderer";
import {GridUtils} from "../../../../utils/GridUtils";
import {CpdSemType} from "../CpdSemType";

export class CpdColumnHeaderRenderer extends EmptyButtonGridColumnRenderer
{
    adjustColumn(cell)
    {
        const grid = cell.grid;

        const dframe = grid.dataFrame;
        const nColCpd = GridUtils.findColumnBySemType(dframe, CpdSemType);
        const colT = dframe.columns.byIndex(nColCpd);
        return colT;
    }
}