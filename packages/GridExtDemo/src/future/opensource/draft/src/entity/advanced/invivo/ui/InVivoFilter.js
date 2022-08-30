import {FilterGroup} from "../../../../ui/filters/FilterGroup";
import {InVivoSemType} from "../InVivoSemType";
import {CheckBoxesFilter} from "../../../../ui/filters/CheckBoxesFilter";
import {EntityColumnCtx} from "../../../ui/EntityColumnCtx";

export class InVivoFilter extends FilterGroup
{
    constructor(ctxCol)
    {
        super(ctxCol);
    }

    isCompatible(ctxCol)
    {
        const typeSem = ctxCol.getSemType();
        const b = typeSem instanceof InVivoSemType;

        return b;
    }

    createFilters(ctxColEntity)
    {
        let ctxCol = new LipinskiColumnCtx(ctxColEntity);
        const cbValue = new CheckBoxesFilter(ctxCol);
        return [cbValue];
    }
}

class LipinskiColumnCtx extends EntityColumnCtx
{
    constructor(ctxCol)
    {
        super(ctxCol);
    }

    getName() {return "Value";}
    getType() {return DG.COLUMN_TYPE.BOOL;}
    getChildValue(entity)
    {
        return entity.hasValue();
    }
}