import * as DG from "datagrok-api/dg";

import {CpdSemType} from "../CpdSemType";
import {RecentTimeFilter} from "../../../../ui/filters/RecentTimeFilter";
import {EntityColumnCtx} from "../../../ui/EntityColumnCtx";
import {MultilineTextFilter} from "../../../../ui/filters/MultilineTextFilter";
import {FilterGroup} from "../../../../ui/filters/FilterGroup";

export class CpdFilter extends FilterGroup
{
    constructor(ctxCol)
    {
        super(ctxCol);
    }

    isCompatible(ctxCol)
    {
        const typeSem = ctxCol.getSemType();
        const b = typeSem instanceof CpdSemType;

        return b;
    }

    createFilters(ctxColEntity)
    {
        let ctxCol = new AnyIDColumnCtx(ctxColEntity);
        const filterSampleID = new MultilineTextFilter(ctxCol);

        ctxCol = new LastTestedColumnCtx(ctxColEntity);
        const sliderLastTested = new RecentTimeFilter(ctxCol);

        return [filterSampleID, sliderLastTested];
    }

}


class AnyIDColumnCtx extends EntityColumnCtx
{
    constructor(ctxCol)
    {
        super(ctxCol);
    }

    getName() {return "Any Compond ID";}
    getType() {return DG.COLUMN_TYPE.STRING;}
    getChildValue(entity)
    {
        let s = "";
        if(entity.getSampleId() !== null && entity.getSampleId() !== "")
            s += entity.getSampleId();

        if(entity.getConceptId() !== null && entity.getConceptId() !== "")
            s += " " + entity.getConceptId();


        if(entity.getCpdSidId() !== null && entity.getCpdSidId() !== "")
            s += " " + entity.getCpdSidId();

        if(entity.getGnfRegId() !== null && entity.getGnfRegId() !== "")
            s += " " + entity.getGnfRegId();

      return s;
    }
}


class LastTestedColumnCtx extends EntityColumnCtx
{
    constructor(ctxCol) {
        super(ctxCol);
    }

    getName() {return "Last Published Date";}
    getType() {return DG.COLUMN_TYPE.DATE_TIME;}
    getChildValue(entity)
    {
        const nTime = entity.getLastTested();
        return nTime;
        //return DG.DateTime.fromMillisecondsSinceEpoch();
    }
}
