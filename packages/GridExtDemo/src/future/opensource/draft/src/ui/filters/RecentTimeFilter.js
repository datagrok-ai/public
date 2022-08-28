import {CombinedRangeSlider} from "./CombinedRangeSlider";
import {NumberOfYearsCondition, NumberOfMonthsCondition, NumberOfWeeksCondition} from "./RangeCondition";
import {ColumnCtx} from "./ColumnCtx";

export class RecentTimeFilter extends CombinedRangeSlider
{
     constructor(ctxCol)
     {
         super(ctxCol, RecentTimeFilter.createConditions(ctxCol));
     }

    isCompatible(ctxCol)
    {
        const b = ctxCol.getType() === DG.COLUMN_TYPE.DATE_TIME;
        if(!b)
            return false;

        const nFirstNonNullIndex = ColumnCtx.getFirstNonNullSortedValueIndex(ctxCol);
        return  nFirstNonNullIndex >= 0;
    }
}


RecentTimeFilter.createConditions = function(ctxCol)
{

    const nFirstNonNullIndex = ColumnCtx.getFirstNonNullSortedValueIndex(ctxCol);
    if(nFirstNonNullIndex < 0)
        return [];

    const nRowMin = nFirstNonNullIndex;
    const nRowMax = ctxCol.getLength() -1;

    const arOrder = ctxCol.getSortedOrder();
    //const f = ctxCol.getValue(arOrder[nRowMin]);
    let fMinimum = ctxCol.getValue(arOrder[nRowMin]);

    if(typeof fMinimum !== "number")
        fMinimum = fMinimum.a;

    let fMaximum = ctxCol.getValue(arOrder[nRowMax]);
    if(typeof fMaximum !== "number")
        fMaximum = fMaximum.a;


    const nTimeNow = new Date().getTime();

    const arCnds = [];

    const fYearCount = ((nTimeNow - fMinimum)/(360*86400000));
    if(fYearCount < 1.0)
    {
        const fMonthCount = ((nTimeNow - fMinimum)/(30*86400000));
        if(fMonthCount < 1.0)
        {
            const fWeekCount = ((nTimeNow - fMinimum)/(604800000));
            if(fWeekCount < 1.0)
                arCnds.push(new NumberOfWeeksCondition(1));
            else//fWeekCount >= 1.0
            {
                const nWeekCount = fWeekCount <= 3.0 ? (fWeekCount === 1.0 || fWeekCount === 2.0 || fWeekCount === 3.0 ? 0 : 1) + Math.floor(fWeekCount) : 3;
                for(var nWeek=1; nWeek<=nWeekCount; ++nWeek)
                {
                    arCnds.push(new NumberOfWeeksCondition(nWeek));
                }
                if(fWeekCount > 3.0)
                    arCnds.push(new NumberOfMonthsCondition(1));
            }
        }
        else //fMonthCount >= 1.0
        {
            arCnds.push(new NumberOfWeeksCondition(1));
            arCnds.push(new NumberOfWeeksCondition(2));
            arCnds.push(new NumberOfWeeksCondition(3));
            const nMonthCount = fMonthCount <= 3.0 ? (fMonthCount === 1.0 || fMonthCount === 2.0 || fMonthCount === 3.0 ? 0 : 1) + Math.floor(fMonthCount) : 3;
            for(var nMonth=1; nMonth<=nMonthCount; ++nMonth)
            {
                arCnds.push(new NumberOfMonthsCondition(nMonth));
            }
            if(fMonthCount > 3.0)
                arCnds.push(new NumberOfMonthsCondition(1 + Math.floor(fMonthCount)));
        }
    }
    else//fYearCount >= 1.0
    {
        arCnds.push(new NumberOfWeeksCondition(1));
        arCnds.push(new NumberOfWeeksCondition(2));
        arCnds.push(new NumberOfWeeksCondition(3));
        arCnds.push(new NumberOfMonthsCondition(1));
        arCnds.push(new NumberOfMonthsCondition(2));
        arCnds.push(new NumberOfMonthsCondition(3));

        const nYearCount = fYearCount <= 3.0 ? (fYearCount === 1.0 || fYearCount === 2.0 || fYearCount === 3.0 ? 0 : 1) + Math.floor(fYearCount) : 3;
        for(var nYear=1; nYear<=nYearCount; ++nYear)
        {
            arCnds.push(new NumberOfYearsCondition(nYear));
        }
        if(fYearCount > 3.0)
            arCnds.push(new NumberOfYearsCondition(1 + Math.floor(fYearCount)));
    }


    return arCnds;
}