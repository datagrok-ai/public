import {FilterGroup} from "../../../ui/filters/FilterGroup";
import {TextDateSemType} from "../TextDateSemType";
import {EntityColumnCtx} from "../../ui/EntityColumnCtx";
import * as DG from "datagrok-api/dg";
import {RecentTimeFilter} from "../../../ui/filters/RecentTimeFilter";
import {IncompatibleFilterError} from "../../../ui/filters/IncompatibleFilterError";
import {CheckBoxesFilter} from "../../../ui/filters/CheckBoxesFilter";
import {MultilineTextFilter} from "../../../ui/filters/MultilineTextFilter";

export class TextDateFilter extends FilterGroup
{
    constructor(ctxCol) {
        super(ctxCol);
    }

    isCompatible(ctxCol) {
        const typeSem = ctxCol.getSemType();
        const b = typeSem instanceof TextDateSemType;

        return b;
    }

    createFilters(ctxColEmtity)
    {
        let ctxCol = new TextColumnCtx(ctxColEmtity);
        let filterText = null;
        try{filterText = new CheckBoxesFilter(ctxCol);}
        catch(e)
        {
            if(e instanceof IncompatibleFilterError)
                filterText = new MultilineTextFilter(ctxCol);
        }

        ctxCol = new TimeColumnCtx(ctxColEmtity);
        let filterTime = null;
        try{filterTime = new RecentTimeFilter(ctxCol);}
        catch(e)
        {
            if(e instanceof IncompatibleFilterError)
                filterTime = new CheckBoxesFilter(ctxCol);
        }

        return [filterText, filterTime];
    }


}


class TextColumnCtx extends EntityColumnCtx
{
    constructor(ctxCol)
    {
        super(ctxCol);
    }

    getName() {return "Value";}
    getType() {return DG.COLUMN_TYPE.STRING;}
    getChildValue(entity)
    {
        return entity.getText();
    }
}

class TimeColumnCtx extends EntityColumnCtx
{
    constructor(ctxCol) {
        super(ctxCol);
    }

    getName() {return "Last Modified";}
    getType() {return DG.COLUMN_TYPE.DATE_TIME;}
    getChildValue(entity)
    {
        const nTime = entity.getTime();
        return nTime;
    }
}