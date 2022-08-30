import {ColumnCtx} from "../../ui/filters/ColumnCtx";
import {PrimitiveSemType} from "../primitive/PrimitiveSemType";
import {MathUtils} from "../../utils/MathUtils";
import {NotImplementedError} from "../../lang/NotImplementedError";

export class EntityColumnCtx extends ColumnCtx
{
    constructor(ctxCol)
    {
        super();

        this.m_ctxCol = ctxCol;

        const setCategories = new Set();
        const nRecordCount = ctxCol.getLength();

        const ar = new Array(nRecordCount);


        let obValue = -1;
        for(var n=0; n<nRecordCount; ++n)
        {
            obValue = this.getValue(n);
            setCategories.add(obValue);

            ar[n] = {v: obValue, idx: n};
        }

        ar.sort((obOne, obTwo) => {

            const bNullOne = MathUtils.isNullValue(obOne.v);
            const bNullTwo = MathUtils.isNullValue(obTwo.v);

            if(bNullOne && !bNullTwo)
                return -1;

            if(!bNullOne && bNullTwo)
                return 1;

            if(bNullOne && bNullTwo)
                return 0;

            return obOne.v < obTwo.v ? -1 : obOne.v > obTwo ? 1 : 0;
        });

        for(var n=0; n<nRecordCount; ++n) {
            ar[n] = ar[n].idx;
        }

        this.m_arCategories = Array.from(setCategories);
        this.m_arCategories = this.m_arCategories.sort();

        this.m_arSortedOrder = ar;
    }

    getSemType() {return PrimitiveSemType.Instance;}
    getCategories()
    {
       return this.m_arCategories;
    }

    getSortedOrder()
    {
        return this.m_arSortedOrder;
    }


    getLength()
    {
        return this.m_ctxCol.getLength();
    }

    getValue(nRecord) {

        const entity = this.m_ctxCol.getValue(nRecord);
        if(entity === null)
            return null;

        return this.getChildValue(entity);
    }

    getChildValue(entity)
    {
        throw new NotImplementedError();
    }

    getDataFrame()
    {
        return this.m_ctxCol.getDataFrame();
    }
}