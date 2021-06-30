<!-- TITLE: Conversion functions -->
<!-- SUBTITLE: -->

# Conversion functions

*Function List:*

- [Boolean](#boolean)
- [DateParse](#dateparse)
- [TimeParse](#TimeParse)
- [TimeSpanParse](#timespanparse)

## <a name="boolean"></a>Boolean(`x`)

Represents `x` as a boolean value.

```javascript
Boolean(10)       // true
Boolean("abc")    // true
Boolean(true)     // true
Boolean(0)        // false
Boolean(null)     // false
Boolean("")       // false
Boolean(false)    // false
```

## <a name="dateparse"></a>DateParse(`s`)

Constructs and returns a date based on string pattern `s`.

The function parses a subset of ISO 8601 which includes the subset accepted by [RFC 3339](https://datatracker.ietf.org/doc/html/rfc3339).

The accepted inputs are currently:

- A date: A signed four-to-six digit year, two digit month and two digit day, optionally separated by - characters. Examples: "19700101", "-0004-12-24", "81030-04-01".

- An optional time part, separated from the date by either T or a space. The time part is a two digit hour, then optionally a two digit minutes value, then optionally a two digit seconds value, and then optionally a '.' or ',' followed by at least a one digit second fraction. The minutes and seconds may be separated from the previous parts by a ':'. Examples: "12", "12:30:24.124", "12:30:24,124", "123010.50".

- An optional time-zone offset part, possibly separated from the previous by a space. The time zone is either 'z' or 'Z', or it is a signed two digit hour part and an optional two digit minute part. The sign must be either "+" or "-", and can not be omitted. The minutes may be separated from the hours by a ':'. Examples: "Z", "-10", "+01:30", "+1130".

```javascript
DateParse("20120227T132700")    // 2012-02-27 13:27:00.000
```

## <a name="timeparse"></a>TimeParse(`s`)

Constructs and returns a time based on string pattern `s`.

The time is actually a DateTime with a special insignificant date `0001-01-01`.

```javascript
TimeParse("13:27")    // 0001-01-01 13:27:00.000
```

## <a name="timespanparse"></a>TimeSpanParse(`s`)

Constructs and returns a TimeSpan based on string pattern `s`.

```javascript
TimeSpanParse("01:00")    //  Creates a time interval of one hour
```
