import {GridUtils} from "../../utils/GridUtils";
import {MathUtils} from "../../utils/MathUtils";
import {NotImplementedError} from "../../lang/NotImplementedError";

export class ColumnCtx
{
    getName()
    {
        throw new NotImplementedError();
    }

    getMin()
    {
        throw new NotImplementedError();
    }

    getMax()
    {
        throw new NotImplementedError();
    }

    getCategories()
    {
        throw new NotImplementedError();
    }

    getSortedOrder()
    {
        throw new NotImplementedError();
    }

    getLength()
    {
        throw new NotImplementedError();
    }

    getValue(nRecord)
    {
        throw new NotImplementedError();
    }

    getType()
    {
        throw new NotImplementedError();
    }

    getSemType()
    {
        throw new NotImplementedError();
    }

    getDataFrame()
    {
        throw new NotImplementedError();
    }
}


ColumnCtx.getFirstNonNullSortedValueIndex = function(ctxCol)
{
    const arSortedIdxs = ctxCol.getSortedOrder();
    const nLength = arSortedIdxs.length;
    let fValue = null;
    for(var n=0; n<nLength; ++n)
    {
        fValue = ctxCol.getValue(arSortedIdxs[n]);
        if(!MathUtils.isNullValue(fValue))
            return n;
    }

    return -1;
}


export class DefaultColumnCtx extends ColumnCtx
{
    constructor(col)
    {
        super();

        this.m_col = col;

        const setCategories = new Set();
        const nRecordCount = this.m_col.length;
        /*
        for(var n=0; n<nRecordCount; ++n)
        {
            setCategories.add(this.m_col.get(n));
        } */



        const ar = new Array(nRecordCount);


        let obValue = -1;
        for(var n=0; n<nRecordCount; ++n)
        {
            obValue = this.getValue(n);
            setCategories.add(obValue);

            ar[n] = {v: obValue, idx: n};
        }

        const arCategories = Array.from(setCategories);
        this.m_arCategories = arCategories.sort();

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

        this.m_arSortedOrder = ar;

    }

    getName() {return this.m_col.name;}

    getMin()
    {
        return this.m_col.min;
    }

    getMax()
    {
        return this.m_col.max;
    }

    getLength()
    {
        return this.m_col.length;
    }

    getCategories()
    {
        return this.m_arCategories;
    }

    getSortedOrder()
    {
        return this.m_arSortedOrder;//this.m_col.getSortedOrder();
    }


    getValue(nRecord)
    {
        return this.m_col.get(nRecord);
    }

    getType()
    {
        return this.m_col.type;
    }

    getSemType()
    {
        const typeSem = GridUtils.getSemType(this.m_col);
        return typeSem;
    }

    getDataFrame()
    {
        const dframe = this.m_col.dataFrame;
        return dframe;
    }
}