import {SemType} from "../SemType";
import {SemSorter} from "../SemSorter";
import {DateEntity} from "./DateEntity";
import {CDFGenericDateType} from "../../service/nx/cdf/composite/CDFGenericDateType";
import {MathUtils} from "../../utils/MathUtils";

export class DateSemType extends SemType {

    constructor(strIconFilePath, arSorters) {
        super(strIconFilePath, arSorters === undefined ? [new DateSemType.DateSorter()] : arSorters);
    }

    getEntityClass()
    {
        return DateEntity;
    }

    getCompatibilityLevel(arCols)
    {
        const arKeys = [CDFGenericDateType.dValue];
        const arHitColIdxs = new Array(arKeys.length);
        let nHitCount = SemType.fillCompatibilityLevelHits(arKeys, arHitColIdxs, arCols);
        return arHitColIdxs[0]  >= 0 ? nHitCount : 0;
   }

    createEntity(arCols, nRow)
    {
        const arKeys = [CDFGenericDateType.dValue];
        const arHitColIdxs = new Array(arKeys.length);
        let nHitCount = SemType.fillCompatibilityLevelHits(arKeys, arHitColIdxs, arCols);
        if(arHitColIdxs[0] < 0)
            throw new Error("Qualified number was not identified.");

        let col = arCols[arHitColIdxs[0]]; //0 must exist

        let nTime = NaN;
        if(arHitColIdxs[0] >= 0)
        {
            col = arCols[arHitColIdxs[0]];
            let date = col.get(nRow);
            if(date !== null && date !== undefined)
                nTime = date.a;
        }

        return new DateEntity(fValue, sign, nTime);
    }
}


DateSemType.DateSorter = class extends SemSorter
{
    getName() {return "Date";}

    isNullSortValue(ob)
    {
        if(super.isNullSortValue(ob))
            return true;

        const nTime = ob.getTime();
        return MathUtils.isNullValue(nTime);
    }

    comparer(ctxOne, ctxTwo)
    {
        let obOne = ctxOne.m_obValue;
        let obTwo = ctxTwo.m_obValue;
        let nTimeOne = obOne.getTime();
        let nTimeTwo = obTwo.geTime();

        let nResult =  nTimeOne < nTimeTwo ? 1 : nTimeOne > nTimeTwo ? -1 : 0;
        return nResult;
    }
}

