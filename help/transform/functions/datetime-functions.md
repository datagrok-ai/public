---
title: "Date and Time functions"
---

*Function List:*

- [Date](#date)
- [DateAdd](#dateadd)
- [DateDiff](#datediff)
- [DateNow](#datenow)
- [DateParse](#dateparse)
- [DateTime](#datetime)
- [DayOfMonth](#dayofmonth)
- [DayOfWeek](#dayofweek)
- [DayOfYear](#dayofyear)
- [Hour](#hour)
- [Millisecond](#millisecond)
- [Minute](#minute)
- [Month](#month)
- [Quarter](#quarter)
- [Second](#second)
- [Time](#time)
- [TimeParse](#timeparse)
- [Today](#today)
- [Weeknum](#weeknum)
- [Year](#year)

## Date(`year`, `month`, `day`) {#date}

Returns a date composed of the specified `year`, `month` and `day`.

```javascript
Date(2002, 3, 15)
```

## DateAdd(`dt`, `ts`) {#dateadd}

Returns the date `dt` with the time span `ts` appended to it.

```javascript
DateAdd(${StartDate}, ${Duration})
```

## DateDiff(`dt1`, `dt2`) {#datediff}

Returns the difference (time span) between two dates `dt1` and `dt2`.

```javascript
DateDiff(${FirstDate}, ${SecondDate})
```

## DateNow() {#datenow}

Returns the current date.

```javascript
DateAdd(DateNow(), ${Duration})
```

## DateParse(`s`) {#dateparse}

Constructs and returns a date based on string pattern `s`.

```javascript
DateParse("20120227T132700")    // 2012-02-27 13:27:00.000
```

## DateTime(`year`, `month`, `day`, `hours`, `minutes`, `seconds`, `milliseconds`) {#datetime}

Returns a date composed of the specified parameters.

```javascript
DateTime(2002, 3, 15, 23, 59, 45, 999)
```

## DayOfMonth(`dt`) {#dayofmonth}

Returns the day of the month for the date `dt`.

```javascript
DayOfMonth(Date(2021, 6, 14))    // 14
```

## DayOfWeek(`dt`) {#dayofweek}

Returns the day of the week for the date `dt`. The days of the week are numbered from 1 to 7.

```javascript
DayOfWeek(Date(2020, 12, 31))    // 4
```

## DayOfYear(`dt`) {#dayofyear}

Returns the day of the year for the date `dt`.

```javascript
DayOfYear(Date(2021, 2, 25))    // 56 (31 + 25)
```

## Hour(`dt`) {#hour}

Returns the hour of the date `dt`.

```javascript
Hour(DateTime(2002, 3, 15, 23, 59, 45, 999))    // 23
```

## Millisecond(`dt`) {#millisecond}

Returns the millisecond of the date `dt`.

```javascript
Millisecond(DateTime(2002, 3, 15, 23, 59, 45, 999))    // 999
```

## Minute(`dt`) {#minute}

Returns the minute of the date `dt`.

```javascript
Minute(DateTime(2002, 3, 15, 23, 59, 45, 999))    // 59
```

## Month(`dt`) {#month}

Returns the month of the date `dt`.

```javascript
Month(DateTime(2002, 3, 15, 23, 59, 45, 999))    // 3
```

## Quarter(`dt`) {#quarter}

Returns the quarter of the date `dt`.

```javascript
Quarter(Date(2002, 3, 15))    // 1
```

## Second(`dt`) {#second}

Returns the second of the date `dt`.

```javascript
Second(DateTime(2002, 3, 15, 23, 59, 45, 999))    // 45
```

## Time(`hours`, `minutes`, `seconds`, `milliseconds`) {#time}

Returns a time composed of the specified parameters.

The time is actually a DateTime with a special insignificant date `0001-01-01`.

```javascript
Time(23, 59, 45, 999)    // 0001-01-01 23:59:45.999
```

## TimeParse(`s`) {#timeparse}

Constructs and returns a time based on string pattern `s`.

The time is actually a DateTime with a special insignificant date `0001-01-01`.

```javascript
TimeParse("13:27")    // 0001-01-01 13:27:00.000
```

## Today() {#today}

Returns today's date.

```javascript
Today()    // On January 1, 2021, would return the value 01-01-2021 00:00:00.000
```

## Weeknum(`dt`) {#weeknum}

Returns the week number of the date `dt`.

```javascript
Weeknum(Date(2021, 2, 3))    // 5
```

## Year(`dt`) {#year}

Returns the year of the date `dt`.

```javascript
Year(Date(2021, 2, 3))    // 2021
```
