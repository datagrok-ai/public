import {TypesUtils} from "./TypesUtils";

export class DateUtils
{
    constructor()
    {
        throw new Error("Cannot create instances of this class");
    }
}

DateUtils.isValidDate = function (dt)
{
  return dt instanceof Date && !isNaN(dt);
}


DateUtils.getHourCount = function(nTimeSpan)
{
    const nSecCount = Math.floor(nTimeSpan/1000);
    const nMinCount = Math.floor(nSecCount/60);
    const nHourCount = Math.floor(nMinCount/60);

    return nHourCount;
}

DateUtils.getMinCount = function(nTimeSpan)
{
    const nSecCount = Math.floor(nTimeSpan/1000);
    const nMinCount = Math.floor(nSecCount/60);

    return nMinCount;
}

DateUtils.getSecCount = function(nTimeSpan)
{
    const nSecCount = Math.floor(nTimeSpan/1000);
    return nSecCount;
}


DateUtils.isRecent = function(nTime, nDayCount)
{
    //const nTime = entity === null ? Number.MIN_SAFE_INTEGER : entity.getLastTested();
    let nDiff = new Date().getTime() - nTime;
    nDiff /= DateUtils.MILLIS_IN_DAY;
    return nDiff < nDayCount;
}


DateUtils.getRecentDateColor = function(nTime)
{
    if(isNaN(nTime))
        return null;

    const cr = DateUtils.isRecent(nTime, DateUtils.LATEST_DAY_COUNT) ? DateUtils.LATEST_DAY_COLOR :
               DateUtils.isRecent(nTime, DateUtils.LATER_DAY_COUNT) ? DateUtils.LATER_DAY_COLOR :
               DateUtils.isRecent(nTime, DateUtils.ONE_DAY_COUNT) ? DateUtils.ONE_DAY_COLOR : null;
    return cr;
}


DateUtils.getTimeFromDGDate = function(dtDG)
{
    let nTime = NaN;

    if(TypesUtils.isString(dtDG))
    {
        if(dtDG === "") {
            nTime = NaN;
        }
        else
        {
            let ghjghj = 0;
        }

        return nTime;
    }


    else if(dtDG !== null && dtDG !== undefined)
        nTime = dtDG.a;

    if(!TypesUtils.isNumber(nTime))
    {
        nTime = NaN;
    }

    return nTime;
}



DateUtils.MILLIS_IN_DAY = 1000 * 60 * 60 * 24;
DateUtils.LATER_DAY_COUNT = 14;
DateUtils.LATEST_DAY_COUNT = 7;
DateUtils.ONE_DAY_COUNT = 1;

DateUtils.LATER_DAY_COLOR = "rgba(255, 243, 21, 0.15)";
DateUtils.LATEST_DAY_COLOR = "rgba(255, 250, 161, 1.0)";
DateUtils.ONE_DAY_COLOR = "rgba(255, 243, 21, 0.85)";
