<!-- TITLE: Conversion functions -->
<!-- SUBTITLE: -->

# Conversion functions

*Function List:*

- [Conversion functions](#conversion-functions)
  - [<a name="boolean"></a>Boolean(`x`)](#booleanx)
  - [<a name="dateparse"></a>DateParse(`s`)](#dateparses)
  - [<a name="timespanparse"></a>TimeSpanParse(`s`)](#timespanparses)


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

```javascript
DateParse("20120227T132700")    // 2012-02-27 13:27:00.000
```

## <a name="timespanparse"></a>TimeSpanParse(`s`)

Constructs and returns a TimeSpan based on string pattern `s`.

```javascript
TimeSpanParse("01:00")    //  Creates a time interval of one hour
```
