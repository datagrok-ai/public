---
title: "TimeSpan functions"
---

TimeSpan (or duration) functions work with time intervals. The internal representation of time intervals is the amount
of time in milliseconds. For example, an interval of 1 hour would be stored as the number 3600000 (the number of
milliseconds in one hour). Time spans include days, hours, minutes, seconds and milliseconds.

*Function List:*

- [InDays](#indays)
- [InHours](#inhours)
- [InMinutes](#inminutes)
- [InSeconds](#inseconds)
- [InMilliseconds](#inmilliseconds)
- [TimeSpan](#timespan)
- [TimeSpanParse](#timespanparse)
- [TotalDays](#totaldays)
- [TotalHours](#totalhours)
- [TotalMinutes](#totalminutes)
- [TotalSeconds](#totalseconds)
- [TotalMilliseconds](#totalmilliseconds)

## InDays(`ts`) {#indays}

Returns the number of full days in time span `ts`.

```javascript
InDays(${Duration})
```

## InHours(`ts`) {#inhours}

Returns the number of full hours days in time span `ts`.

```javascript
InHours(${Duration})
```

## InMinutes(`ts`) {#inminutes}

Returns the number of full minutes days in time span `ts`.

```javascript
InMinutes(${Duration})
```

## InSeconds(`ts`) {#inseconds}

Returns the number of full seconds days in time span `ts`.

```javascript
InSeconds(${Duration})
```

## InMilliseconds(`ts`) {#inmilliseconds}

Returns the number of milliseconds days in time span `ts`.

```javascript
InMilliseconds(${Duration})
```

## TimeSpan(`days`, `hours`, `minutes`, `seconds`, `milliseconds`) {#timespan}

Creates a TimeSpan from specified parameters.

```javascript
TimeSpan(0, 1, 30, 0, 0)    // Time interval of 1.5 hours
```

## TimeSpanParse(`s`) {#timespanparse}

Constructs and returns a TimeSpan based on string pattern `s`.

```javascript
TimeSpanParse("01:00")    //  Creates a time interval of one hour
```

## TotalDays(`ts`) {#totaldays}

Returns the number of full and fractional parts of days in the TimeSpan `ts`.

```javascript
TotalDays(TimeSpan(1, 12, 0, 0, 0))    // Time interval of 1.5 days
```

## TotalHours(`ts`) {#totalhours}

Returns the number of full and fractional parts of hours in the TimeSpan `ts`.

```javascript
TotalHours(TimeSpan(0, 1, 30, 0, 0))    // Time interval of 1.5 hours
```

## TotalMinutes(`ts`) {#totalminutes}

Returns the number of full and fractional parts of minutes in the TimeSpan `ts`.

```javascript
TotalMinutes(TimeSpan(0, 0, 2, 30, 0))    // Time interval of 2.5 minutes
```

## TotalSeconds(`ts`) {#totalseconds}

Returns the number of full and fractional parts of seconds in the TimeSpan `ts`.

```javascript
TotalSeconds(TimeSpan(0, 0, 0, 7, 500))    // Time interval of 7.5 seconds
```

## TotalMilliseconds(`ts`) {#totalmilliseconds}

Returns the number of full and fractional parts of milliseconds in the TimeSpan `ts`.

```javascript
TotalMilliseconds(TimeSpan(0, 0, 0, 7, 500))    // Time interval of 7500 milliseconds
```
