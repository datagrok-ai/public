---
title: "Binning functions"
---

Binning is a way to group a number of more or less continuous values into a smaller number of "bins"
. For example, if you have data about a group of people, you might want to arrange their ages into a smaller number of
age intervals.

*Function List:*

- [BinByDateTime](#binbydatetime)
- [BinBySpecificLimits](#binbyspecificlimits)

## <a name="binbydatetime"></a>BinByDateTime(`dt`, `levels`, `levelIndex`)

Groups the values into bins based on a datetime hierarchy.

`dt` is the DateTime value to bin. `levels` (string) is the definition of the levels in the DateTime hierarchy. The
hierarchy levels should be written in the form of a string containing the desired date parts, separated by dots, for
example "yy.qq.mm". `levelIndex` (int) is the pruning level which specifies the level of the hierarchy to
display. `levelIndex` numbering starts at 0.

Valid arguments for `levels` are combinations of:

- "yy" - The year.
- "qq" - The quarter.
- "mm" - The month.
- "dy" - The day of year.
- "dd" - The day.
- "wk" - The week.
- "dw" - The weekday.
- "hh" - The hour.
- "mi" - The minute.
- "ss" - The second.
- "ms" - The millisecond.

```javascript
BinByDateTime(${Age}, "yy.qq.mm", 2)
BinByDateTime(Date(1970, 11, 17), "qq.mm.dd", 0)    // 4
```

## <a name="binbyspecificlimits"></a>BinBySpecificLimits(`x`, `limits`)

Groups the values in the column by defined limits for the bins. `x` is the checked value to bin and the following
arguments are the `limits` for the bins. `limits` here is a list of numbers.

```javascript
BinBySpecificLimits(${Age}, [18, 30, 45, 60, 75])
BinBySpecificLimits(20, [18, 30, 45])    // "18 < x <= 30"
```
