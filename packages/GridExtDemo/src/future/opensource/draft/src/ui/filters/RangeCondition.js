import {NotImplementedError} from "../../lang/NotImplementedError";

export class RangeCondition
{
    constructor() {
        this.m_bIncludeNulls = false;
    }

    fillRangeBounds(fMinimum, fMaximum, arMinMaxs)
    {
        throw new NotImplementedError();
    }
}


export class AllRangeCondition extends RangeCondition
{
   fillRangeBounds(fMinimum, fMaximum, arMinMaxs)
    {
        arMinMaxs[0] = Number.NEGATIVE_INFINITY;
        arMinMaxs[1] = Number.POSITIVE_INFINITY;
        return true;
    }

    toString() {return "All";}
}


export class NumberOfTimeUnitsCondition extends RangeCondition
{
    constructor(nUnitCount, nRangeFromNow, strUnitName)
    {
    super();

    if(nUnitCount < 0)
     throw new Error("The number of units cannot be negative.");

    this.m_nUnitCount = nUnitCount;
    this.m_strUnitName = strUnitName;
    this.m_nRangeFromNow = nRangeFromNow;
    }

    fillRangeBounds(fMinimum, fMaximum, arMinMaxs)
    {
    const nTimeMin = new Date().getTime() - this.m_nRangeFromNow*this.m_nUnitCount;
    arMinMaxs[0] = nTimeMin;
    if(fMaximum < arMinMaxs[0])
    {
        arMinMaxs[0] = NaN;
        arMinMaxs[1] = NaN;
        return false;
    }

    arMinMaxs[0] = Math.max(fMinimum, arMinMaxs[0]);
    arMinMaxs[1] = fMaximum;

    return true;
    }

    toString()
    {
        return this.m_nUnitCount + " " + this.m_strUnitName + (this.m_nUnitCount === 1 ? "" : "s");
    }

}

export class NumberOfWeeksCondition extends NumberOfTimeUnitsCondition
{
    constructor(nWeekCount)
{
    super(nWeekCount, 604800000, "Week");
}
}

export class NumberOfMonthsCondition extends NumberOfTimeUnitsCondition
{
    constructor(nMonthCount)
{
    super(nMonthCount, 30*86400000, "Month");
}
}

export class NumberOfYearsCondition extends NumberOfTimeUnitsCondition
{
    constructor(nYearCount)
{
    super(nYearCount, 360*86400000, "Year");
}
}
