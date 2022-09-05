import {SemType} from "../SemType";
import {SemSorter} from "../SemSorter";
import {QNumDateEntity} from "../qnum/QNumDateEntity";

export class PrimitiveSemType extends SemType
{
    constructor(strIconFilePath)
    {
        super(strIconFilePath, [new PrimitiveSemType.ValueSorter()]);
    }

    getEntityClass()
    {
        return object;
    }
}


PrimitiveSemType.ValueSorter = class extends SemSorter {
    getName() {return "Value";}
   /*
    isNullSortValue(obValue)
    {
        if(super.isNullSortValue(obValue))
            return true;

        const fValue = obValue.getValue();
        return super.isNullSortValue(fValue);
    }*/

    comparer(ctxOne, ctxTwo)
    {
        const obOne = ctxOne.m_obValue;
        const obTwo = ctxTwo.m_obValue;
        const valOne = obOne;
        const valTwo = obTwo;

        const nResult =  valOne < valTwo ? -1 : valOne > valTwo ? 1 : 0;
        return nResult;
    }
}


PrimitiveSemType.Instance = new PrimitiveSemType();
