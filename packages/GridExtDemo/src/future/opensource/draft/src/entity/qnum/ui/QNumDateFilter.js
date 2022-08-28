import {FilterGroup} from "../../../ui/filters/FilterGroup";
import {QNumDateSemType} from "../QNumDateSemType";
import {EntityColumnCtx} from "../../ui/EntityColumnCtx";
import * as DG from "datagrok-api/dg";
import {RangeSliderFilter} from "../../../ui/filters/RangeSliderFilter";
import {RecentTimeFilter} from "../../../ui/filters/RecentTimeFilter";
import {IncompatibleFilterError} from "../../../ui/filters/IncompatibleFilterError";
import {CheckBoxesFilter} from "../../../ui/filters/CheckBoxesFilter";

export class QNumDateFilter extends FilterGroup
{
    constructor(ctxCol) {
        super(ctxCol);
    }

    isCompatible(ctxCol) {
        const typeSem = ctxCol.getSemType();
        const b = typeSem instanceof QNumDateSemType;

        return b;
    }

    createFilters(ctxColEmtity)
    {
        let ctxCol = new ValueColumnCtx(ctxColEmtity);
        const sliderValue = new RangeSliderFilter(ctxCol);

        ctxCol = new TimeColumnCtx(ctxColEmtity);
        let filter = null;
        try{filter = new RecentTimeFilter(ctxCol);}
        catch(e)
        {
            if(e instanceof IncompatibleFilterError)
             filter = new CheckBoxesFilter(ctxCol);
        }

        return [sliderValue, filter];
    }


}


class ValueColumnCtx extends EntityColumnCtx
{
    constructor(ctxCol)
    {
        super(ctxCol);
    }

    getName() {return "Value";}
    getType() {return DG.COLUMN_TYPE.FLOAT;}
    getChildValue(entity)
    {
        return entity.getValue();
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