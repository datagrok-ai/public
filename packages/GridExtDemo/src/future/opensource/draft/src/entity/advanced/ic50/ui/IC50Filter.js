import * as DG from "datagrok-api/dg";
import {FilterGroup} from "../../../../ui/filters/FilterGroup";
import {IC50SemType} from "../IC50SemType";
import {RangeSliderFilter} from "../../../../ui/filters/RangeSliderFilter";
import {EntityColumnCtx} from "../../../ui/EntityColumnCtx";
import {RecentTimeFilter} from "../../../../ui/filters/RecentTimeFilter";

export class IC50Filter extends FilterGroup
{
    constructor(ctxCol)
    {
        super(ctxCol);
    }

    isCompatible(ctxCol)
    {
        const typeSem = ctxCol.getSemType();
        const b = typeSem instanceof IC50SemType;

        return b;
    }


    createFilters(ctxColEmtity)
    {
        let ctxCol = new MedianIC50ColumnCtx(ctxColEmtity);
        const sliderMedIC50 = new RangeSliderFilter(ctxCol);

        ctxCol = new PctEffColumnCtx(ctxColEmtity);
        const sliderPctEff = new RangeSliderFilter(ctxCol);

        ctxCol = new LastTestedColumnCtx(ctxColEmtity);
        const sliderLastTested = new RecentTimeFilter(ctxCol);

        return [sliderMedIC50, sliderPctEff, sliderLastTested];
    }
}

class MedianIC50ColumnCtx extends EntityColumnCtx
{
    constructor(ctxCol)
    {
        super(ctxCol);
    }

    getName() {return "Median IC50";}
    getType() {return DG.COLUMN_TYPE.FLOAT;}
    getChildValue(entity)
    {
       return entity.getMedianIC50();
    }
}


class PctEffColumnCtx extends EntityColumnCtx
{
    constructor(ctxCol)
    {
        super(ctxCol);
    }

    getName() {return "Pct Efficacy";}
    getType() {return DG.COLUMN_TYPE.FLOAT;}
    getChildValue(entity)
    {
        return entity.getPctEff();
    }
}


class LastTestedColumnCtx extends EntityColumnCtx
{
    constructor(ctxCol) {
        super(ctxCol);
    }

    getName() {return "Last Tested";}
    getType() {return DG.COLUMN_TYPE.DATE_TIME;}
    getChildValue(entity)
    {
        const nTime = entity.getLastTested();
        return nTime;
       //return DG.DateTime.fromMillisecondsSinceEpoch();
    }
}