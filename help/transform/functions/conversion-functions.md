<!-- TITLE: Conversion functions -->
<!-- SUBTITLE: -->

# Conversion functions

*Function List:*

- [Boolean](#booleanx)
- [DateParse](#dateparses)
- [TimeSpanParse](#timespanparses)

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
Boolean("n")        // false

Boolean(null)       // null
Boolean("")         // null
Boolean("abc")      // null
```

## <a name="dateparse"></a>DateParse(`s`)

Constructs and returns a date based on string pattern `s`.

```javascript
DateParse("20120227T132700")    // 2012-02-27 13:27:00.000
```

## <a name="timespanparse"></a>TimeSpanParse(`s`)

Constructs and returns a TimeSpan based on string pattern `s`.

```javascript
TimeSpanParse("01:00")    //  Creates a time interval of one hour
```
