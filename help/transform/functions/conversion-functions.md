---
title: "Conversion functions"
---

*Function List:*

- [Boolean](#boolean)
- [DateParse](#dateparse)
- [ParseQnum](#parseqnum)
- [Qnum](#qnum)
- [QnumToDouble](#qnumtodouble)
- [QnumToString](#qnumtostring)
- [Qualifier](#qualifier)
- [TimeSpanParse](#timespanparse)
- [ToString](#tostring)

## Boolean

`Boolean(`x`)`: Represents `x` as a boolean value.

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

## DateParse

`DateParse(`s`)`: Constructs and returns a date based on string pattern `s`.

```javascript
DateParse("20120227T132700")    // 2012-02-27 13:27:00.000
```

## ParseQnum

`ParseQnum(`s`)`: Parses a qualified number from string `s`.

```javascript
ParseQnum("10")
ParseQnum("<10")
ParseQnum(" > 10")
```

## Qnum

`Qnum(`x`, `q`)`: Returns a qualified number with value `x` and qualifier `q`.

```javascript
Qnum(1, "=")
Qnum(2, "<")
Qnum(3, ">")
```

## QnumToDouble

`QnumToDouble(`x`)`: Converts a qualified number to a double precision floating point value.

```javascript
QnumToDouble(Qnum(1.5, "<"))  // 1.5
```

## QnumToString

`QnumToString(`x`)`: Converts a qualified number to a string representation.

```javascript
QnumToString(Qnum(1.5, "<"))  // <1.50
```

## Qualifier

Qualifier(`x`)`: Returns the qualifier character from a qualified number.

```javascript
Qualifier(Qnum(1.5, "<"))  // "<"
Qualifier(Qnum(1.5, ">"))  // ">"
Qualifier(Qnum(1.5, "="))  // "="
```

## TimeSpanParse

`TimeSpanParse(`s`)`: Constructs and returns a TimeSpan based on string pattern `s`.

```javascript
TimeSpanParse("01:00")    //  Creates a time interval of one hour
```

## ToString

`ToString(`x`)`: Returns a string representation of number `x`.

```javascript
ToString(3.14)    //  "3.14"
```
