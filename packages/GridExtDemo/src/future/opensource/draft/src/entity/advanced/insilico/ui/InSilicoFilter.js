import {FilterGroup} from "../../../../ui/filters/FilterGroup";
import {RangeSliderFilter} from "../../../../ui/filters/RangeSliderFilter";
import {InSilicoSemType} from "../InSilicoSemType";
import {CheckBoxesFilter} from "../../../../ui/filters/CheckBoxesFilter";
import {EntityColumnCtx} from "../../../ui/EntityColumnCtx";

export class InSilicoFilter extends FilterGroup
{
    constructor(ctxCol)
    {
        super(ctxCol);
    }

    isCompatible(ctxCol)
    {
        const typeSem = ctxCol.getSemType();
        const b = typeSem instanceof InSilicoSemType;

        return b;
    }


    createFilters(ctxColEntity)
    {
        let ctxCol = new LipinskiColumnCtx(ctxColEntity);
        const cbViolations = new CheckBoxesFilter(ctxCol);

        ctxCol = new MWColumnCtx(ctxColEntity);
        const sliderMW = new RangeSliderFilter(ctxCol);

        return [cbViolations, sliderMW];
    }
}

class LipinskiColumnCtx extends EntityColumnCtx
{
    constructor(ctxCol)
    {
        super(ctxCol);
    }

    getName() {return "Lipinski Violations";}
    getType() {return DG.COLUMN_TYPE.INT;}
    getChildValue(entity)
    {
       return entity.getLipinskiViolations();
    }
}

class MWColumnCtx extends EntityColumnCtx
{
    constructor(ctxCol)
    {
        super(ctxCol);
    }

    getName() {return "Mol Weight";}
    getType() {return DG.COLUMN_TYPE.FLOAT;}

    getChildValue(entity)
    {
        return entity.getMW();
    }
}

