import {FilterGroup} from "../../../ui/filters/FilterGroup";
import {EntityColumnCtx} from "../../ui/EntityColumnCtx";
import * as DG from "datagrok-api/dg";
import {IncompatibleFilterError} from "../../../ui/filters/IncompatibleFilterError";
import {CheckBoxesFilter} from "../../../ui/filters/CheckBoxesFilter";
import {MultilineTextFilter} from "../../../ui/filters/MultilineTextFilter";
import {ObjectSemType} from "../ObjectSemType";

export class ObjectFilter extends FilterGroup
{
    constructor(ctxCol) {
        super(ctxCol);
    }

    isCompatible(ctxCol) {
        const typeSem = ctxCol.getSemType();
        const b = typeSem instanceof ObjectSemType;

        return b;
    }

    createFilters(ctxColEmtity)
    {
        let ctxCol = new ValueColumnCtx(ctxColEmtity);
        let filterText = null;
        try{filterText = new CheckBoxesFilter(ctxCol);}
        catch(e)
        {
            if(e instanceof IncompatibleFilterError)
                filterText = new MultilineTextFilter(ctxCol);
        }

        return [filterText];
    }


}


class ValueColumnCtx extends EntityColumnCtx
{
    constructor(ctxCol)
    {
        super(ctxCol);
    }

    getName() {return "Value";}
    getType() {return DG.COLUMN_TYPE.STRING;}
    getChildValue(entity)
    {
        return entity.getValue();
    }
}
