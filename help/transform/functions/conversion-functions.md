---
title: "Conversion functions"
---

*Function List:*

- [Boolean](#booleanx)
- [DateParse](#dateparses)
- [ParseQnum](#parseqnum)
- [Qnum](#qnum)
- [QnumToString](#qnumtostring)
- [TimeSpanParse](#timespanparses)
- [ToString](#tostringx)

## <a name="boolean"></a>Boolean(`x`)

Represents `x` as a boolean value.

```javascript
Boolean(10)         // true
Boolean(true)       // true
Boolean("true")     // true
Boolean("y")        // true

Boolean(0)          // false
Boolean(false)      // false
Boolean("false")    // false
Boolean("abc")      // false
Boolean("n")        // false
Boolean(null)       // false
Boolean("")         // false
```

## <a name="dateparse"></a>DateParse(`s`)

Constructs and returns a date based on string pattern `s`.

```javascript
DateParse("20120227T132700")    // 2012-02-27 13:27:00.000
```

## <a name="parseqnum"></a>ParseQnum(`s`)

Parses a qualified number from string `s`.

```javascript
ParseQnum("10")
ParseQnum("<10")
ParseQnum(" > 10")
```

## <a name="qnum"></a>Qnum(`x`, `q`)

Returns a qualified number with value `x` and qualifier `q`.

```javascript
Qnum(1, "=")
Qnum(2, "<")
Qnum(3, ">")
```

## <a name="qnumtostring"></a>QnumToString(`x`)

Converts a qualified number to a string representation.

```javascript
QnumToString(Qnum(1.5, "<"))  // <1.5
```

## <a name="timespanparse"></a>TimeSpanParse(`s`)

Constructs and returns a TimeSpan based on string pattern `s`.

```javascript
TimeSpanParse("01:00")    //  Creates a time interval of one hour
```

## <a name="tostring"></a>ToString(`x`)

Returns a string representation of number `x`.

```javascript
ToString(3.14)    //  "3.14"
```
