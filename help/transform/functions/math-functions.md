<!-- TITLE: Math functions -->
<!-- SUBTITLE: -->

# Math functions

This type of function includes math, trigonometric, and logic functions.

As parameters of the function, you can pass numeric scalars, numeric functions, [math constants](constants.md), boolean
scalars, boolean functions, or a column name. To pass a column cell, you can use the syntax `${columnName}`. Other ways
to use parameters: `true`
, `false`, `PI`, `E` etc.

> By default, real numbers in a column show only 2 digits after the integer part.
> You can change this behavior by setting the appropriate format for displaying the data in the column.
> To do this, click on the column heading and in the opened context menu select "Format" option.

*Function List:*

- [Abs](#abs)
- [Acos](#acos)
- [Add](#add)
- [And](#and)
- [Asin](#asin)
- [Atan](#atan)
- [Atan2](#atan2)
- [Avg](#avg)
- [Ceil](#ceil)
- [Cos](#cos)
- [Div](#div)
- [Eq](#eq)
- [Exp](#exp)
- [Floor](#floor)
- [Greater](#greater)
- [Ln](#ln)
- [Log](#log)
- [Log10](#log10)
- [Max](#max)
- [Median](#median)
- [Min](#min)
- [Mod](#mod)
- [Mul](#mul)
- [Neg](#neg)
- [Not](#not)
- [NotEq](#noteq)
- [NotGreater](#notgreater)
- [NotSmaller](#notsmaller)
- [Or](#or)
- [Percentile](#percentile)
- [Pow](#pow)
- [RandBetween](#randbetween)
- [Rnd](#rnd)
- [Round](#round)
- [Round10](#round10)
- [Sin](#sin)
- [Smaller](#smaller)
- [Sqrt](#sqrt)
- [Sub](#sub)
- [Tan](#tan)
- [Xor](#xor)

## <a name="abs"></a>Abs(`x`)

Returns the absolute value of a number [x].

```javascript
Abs(-10)    // 10
```

## <a name="acos"></a>Acos(`x`)

Returns the arccosine of a number [x] as an angle expressed in radians in the interval 0..PI. `x`
must be in the interval -1..1.

```javascript
Acos(0.5)    // 1.047197580337524
```

## <a name="add"></a>Add(`x`, `y`)

Returns the sum of two numbers `x` and `y`.

```javascript
Add(24.06, 100)    // 124.06
```

## <a name="and"></a>And(`a`, `b`)

Returns logical conjunction of boolean `a` and `b`.

```javascript
And(true, false)        // false
And(5 == 5, 10 < 20)    // true
```

## <a name="asin"></a>Asin(`x`)

Returns the arcsine of a number `x` as an angle expressed in radians in the interval -PI/2..PI/2. `x` must be in the
interval -1..1.

```javascript
Asin(0.5)    // 0.479425549507141
```

## <a name="atan"></a>Atan(`x`)

Returns the arc tangent of a number `x` as an angle expressed in radians in the interval -PI/2..PI/2.

```javascript
Atan(1)    // 0.785398185253143
```

## <a name="atan2"></a>Atan2(`x`, `y`)

Returns the angle in radians between the positive X-axis and the vector (`y`,`x`). The result is in the range -PI..PI.
If `y` is positive, this is the same as `Atan(x / y)`.

The result is negative when `x` is negative (including when `x` is -0.0). If `x` is equal to zero, the vector (`y`,`x`)
is considered parallel to the X-axis, even if `y` is also equal to zero. The sign of `y` determines the direction of the
vector along the X-axis.

```javascript
Atan2(0, -0)    // 3.141592741012573
```

## <a name="avg"></a>Avg([`x1`, `x2`, `x3`...])

Returns the average (arithmetic mean) of the numbers. Arguments are a set of numbers enclosed in square brackets.

```javascript
Avg([-1, 3.5, 6.5])    // 3
Avg([1, 2, 3, 4])      // 2.5
```

## <a name="ceil"></a>Ceil(`x`)

Returns the least integer no smaller than `x`.

```javascript
Ceil(3.5)     // 4
Ceil(-3.5)    // -3
```

## <a name="cos"></a>Cos(`x`)

Returns the cosine of a number [x].

```javascript
Cos(0)         // 1
Cos(PI / 3)    // 0.5
```

## <a name="div"></a>Div(`x`, `y`)

Returns the result of dividing `x` by `y`.

```javascript
Div(7.5, 2)    // 3.75
```

## <a name="eq"></a>Eq(`x`, `y`)

Returns true if `x` equal to `y` and false otherwise.

```javascript
Eq(5, 5)           // true
Eq(true, false)    // false
```

## <a name="exp"></a>Exp(`x`)

Returns the natural exponent (E) to the power `x`.

```javascript
Exp(2)    // 7.389056205749512
```

## <a name="floor"></a>Floor(`x`)

Returns the greatest integer no greater than `x`.

```javascript
Floor(3.5)     // 3
Floor(-3.5)    // 4
```

## <a name="greater"></a>Greater(`x`, `y`)

Returns true if `x` is greater than `y` and false otherwise.

```javascript
Greater(5, 5)    // false
Greater(5, 4)    // true
```

## <a name="ln"></a>Ln(`x`)

Returns the natural logarithm of `x`.

```javascript
Ln(0)    // 0
Ln(E)    // 1
```

## <a name="log"></a>Log(`x`, `base`)

Returns the logarithm of `x` expressed in the base specified by `base`.

```javascript
Log(25, 5)    // 2
```

## <a name="log10"></a>Log10(`x`)

Returns the 10-based logarithm of `x`.

```javascript
Log10(100)    // 2
```

## <a name="max"></a>Max([`x`, `y` ...])

Returns the maximum of `x` and `y`.

```javascript
Max([15, 21])    // 21
```

## <a name="median"></a>Median([`x1`, `x2`, `x3`...])

Calculates the median of the numbers. Arguments are a set of numbers enclosed in square brackets.

```javascript
Median([0, 2, 5])       // 2
Median([0, 2, 5, 9])    // 2.5
```

## <a name="min"></a>Min([`x`, `y` ...])

Returns the minimum of `x` and `y`.

```javascript
Min([15, 21])    // 15
```

## <a name="mod"></a>Mod(`x`, `y`)

Returns the remainder of dividing `x` by `y`.

```javascript
Mod(8, 3)    // 2
```

## <a name="mul"></a>Mul(`x`, `y`)

Returns the product of `x` and `y`.

```javascript
Mul(10, 1.5)    // 15
```

## <a name="neg"></a>Neg(`x`)

Returns `x` with opposite sign.

```javascript
Neg(-5)    // 5
Neg(12)    // -12
```

## <a name="not"></a>Not(`a`)

Returns logical negation of the `a`.

```javascript
Not(true)     // false
Not(5 < 1)    // true
```

## <a name="noteq"></a>NotEq(`x`, `y`)

Returns false if `x` equal to `y` and true otherwise.

```javascript
NotEq(5, 5)           // false
NotEq(true, false)    // true
```

## <a name="notgreater"></a>NotGreater(`x`, `y`)

Returns true if `x` is less than or equal to `y` and false otherwise.

```javascript
NotGreater(5, 5)    // true
NotGreater(6, 5)    // false
```

## <a name="notsmaller"></a>NotSmaller(`x`, `y`)

Returns true if `x` is greater than or equal to `y`  and false otherwise.

```javascript
NotSmaller(5, 5)    // true
NotSmaller(5, 6)    // false
```

## <a name="or"></a>Or(`a`, `b`)

Returns logical disjunction of boolean `a` and `b`.

```javascript
Or(true, false)        // true
Or(5 == 6, 20 < 10)    // false
```

## <a name="percentile"></a>Percentile(`nums`, `percentage`)

Returns a value from `nums` below which a given `percentage` of values fall. `percentage` is in the range 0..1

```javascript
Percentile(([1, 2, 3, 4], 0.25))     // 2
Percentile(([1, 2, 3, 4], 0.75))     // 4
```

## <a name="pow"></a>Pow(`x`, `exponent`)

Returns `x` to the power of `exponent`.

```javascript
Pow(2, 3)     // 8
Pow(2, -2)    // 0.25
```

## <a name="randbetween"></a>RandBetween(`m`, `n`)

Returns a random integer number within the range from `m`, inclusive, to `n`, exclusive.

```javascript
RandBetween(5, 7)    // Randomly returns 5 or 6
```

## <a name="rnd"></a>Rnd(`limit`)

Returns a random integer number within the range from 0, inclusive, to `n`, exclusive.

```javascript
Rnd(2)    // Randomly returns 0 or 1
```

## <a name="round"></a>Round(`x`)

Returns the integer closest to `x`. Function rounds away from zero when there is no closest integer.

```javascript
Round(3.4)     // 3
Round(3.5)     // 4
Round(-3.5)    // -4
```

## <a name="round10"></a>Round10(`x`, `decimalPlaces`)

Returns the number rounded up `x` to the number of decimal places specified by `decimalPlaces`.

`decimalPlaces` can be negative to round to even 10s, 100s, etc.

```javascript
Round10(PI, 2)      // 3.14
Round10(25, -1)     // 30
```

## <a name="sin"></a>Sin(`x`)

Returns the sine of the `x`.

```javascript
Sin(0)         // 0
Sin(PI / 6)    // 0.5
```

## <a name="smaller"></a>Smaller(`x`, `y`)

Returns true if `x` is less than `y` and false otherwise.

```javascript
Smaller(5, 5)    // false
Smaller(5, 6)    // true
```

## <a name="sqrt"></a>Sqrt(`x`)

Returns the square root of the `x`.

```javascript
Sqrt(6.25)    // 2.5
```

## <a name="sub"></a>Sub(`x`, `y`)

Returns the difference between `x` and `y`.

```javascript
Sub(10, 3)    // 7
```

## <a name="tan"></a>Tan(`x`)

Returns the tangent of the `x`.

```javascript
Tan(0)         // 0
Tan(PI / 6)    // 0.5
```

## <a name="xor"></a>Xor(`a`, `b`)

Returns logical exclusive disjunction of boolean `a` and `b`.

```javascript
Xor(true, false)        // true
Xor(5 == 5, 10 < 20)    // false
```
