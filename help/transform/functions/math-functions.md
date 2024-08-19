---
title: "Math functions"
---

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
- [Fixed](#fixed)
- [Floor](#floor)
- [FormatFloat](#formatfloat)
- [Greater](#greater)
- [If](#if)
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
- [Qualifier](#qualifier)
- [RandBetween](#randbetween)
- [Rnd](#rnd)
- [Round](#round)
- [RoundFloat](#roundfloat)
- [Sin](#sin)
- [Smaller](#smaller)
- [Sqrt](#sqrt)
- [Sub](#sub)
- [Tan](#tan)
- [Xor](#xor)

## Abs(`x`) {#abs}

Returns the absolute value of a number [x].

```javascript
Abs(-10)    // 10
```

## Acos(`x`) {#acos}

Returns the arccosine of a number [x] as an angle expressed in radians in the interval 0..PI. `x`
must be in the interval -1..1.

```javascript
Acos(0.5)    // 1.047197580337524
```

## Add(`x`, `y`) {#add}

Returns the sum of two numbers `x` and `y`.

```javascript
Add(24.06, 100)    // 124.06
```

## And(`a`, `b`) {#and}

Returns logical conjunction of boolean `a` and `b`.

```javascript
And(true, false)        // false
And(5 == 5, 10 < 20)    // true
```

## Asin(`x`) {#asin}

Returns the arcsine of a number `x` as an angle expressed in radians in the interval -PI/2..PI/2. `x` must be in the
interval -1..1.

```javascript
Asin(0.5)    // 0.479425549507141
```

## Atan(`x`) {#atan}

Returns the arc tangent of a number `x` as an angle expressed in radians in the interval -PI/2..PI/2.

```javascript
Atan(1)    // 0.785398185253143
```

## Atan2(`x`, `y`) {#atan2}

Returns the angle in radians between the positive X-axis and the vector (`y`,`x`). The result is in the range -PI..PI.
If `y` is positive, this is the same as `Atan(x / y)`.

The result is negative when `x` is negative (including when `x` is -0.0). If `x` is equal to zero, the vector (`y`,`x`)
is considered parallel to the X-axis, even if `y` is also equal to zero. The sign of `y` determines the direction of the
vector along the X-axis.

```javascript
Atan2(0, -0)    // 3.141592741012573
```

## Avg([`x1`, `x2`, `x3`...]) {#avg}

Returns the average (arithmetic mean) of the numbers. Arguments are a set of numbers enclosed in square brackets.

```javascript
Avg([-1, 3.5, 6.5])    // 3
Avg([1, 2, 3, 4])      // 2.5
```

## Ceil(`x`) {#ceil}

Returns the least integer no smaller than `x`.

```javascript
Ceil(3.5)     // 4
Ceil(-3.5)    // -3
```

## Cos(`x`) {#cos}

Returns the cosine of a number [x].

```javascript
Cos(0)         // 1
Cos(PI / 3)    // 0.5
```

## Div(`x`, `y`) {#div}

Returns the result of dividing `x` by `y`.

```javascript
Div(7.5, 2)    // 3.75
```

## Eq(`x`, `y`) {#eq}

Returns true if `x` equal to `y` and false otherwise.

```javascript
Eq(5, 5)           // true
Eq(true, false)    // false
```

## Exp(`x`) {#exp}

Returns the natural exponent (E) to the power `x`.

```javascript
Exp(2)    // 7.389056205749512
```

## Fixed(`x`, `decimalPlaces`) {#fixed}

Returns the rounded number `x` to the specified number of `decimalPlaces`.

```javascript
Fixed(3.5, 2)   // 3.50
Fixed(-3.5, 0)  // -4
```

## Floor(`x`) {#floor}

Returns the greatest integer no greater than `x`.

```javascript
Floor(3.5)     // 3
Floor(-3.5)    // -4
```

## FormatFloat(`x`, `format`) {#formatfloat}

Returns the number `x` formatted according to the specified format.

```javascript
FormatFloat(12345.12345, 'scientific')    // 1.23E4
FormatFloat(12345.12345, '#0.00')         // 12345.12
```

## Greater(`x`, `y`) {#greater}

Returns true if `x` is greater than `y` and false otherwise.

```javascript
Greater(5, 5)    // false
Greater(5, 4)    // true
```

## If(`condition`, `ifTrue`, `ifFalse`) {#if}

Returns `ifTrue`, if `condition` is true, or `ifFalse` otherwise.

```javascript
If(true, "a", "b")                 // "a"
If(false, "a", "b")                // "b"
If(true, If(true, "a", "b"), "c")  // "a"
If(Eq(10, 50), 1, 0)               // 0
```

## Ln(`x`) {#ln}

Returns the natural logarithm of `x`.

```javascript
Ln(1)    // 0
Ln(E)    // 1
```

## Log(`x`, `base`) {#log}

Returns the logarithm of `x` expressed in the base specified by `base`.

```javascript
Log(25, 5)    // 2
```

## Log10(`x`) {#log10}

Returns the 10-based logarithm of `x`.

```javascript
Log10(100)    // 2
```

## Max([`x1`, `x2`, `x3` ...]) {#max}

Returns the maximum value from the specified array of numbers.

```javascript
Max([15, 21])    // 21
```

## Median([`x1`, `x2`, `x3`...]) {#median}

Calculates the median of the numbers. Arguments are a set of numbers enclosed in square brackets.

```javascript
Median([0, 2, 5])       // 2
Median([0, 2, 5, 9])    // 2.5
```

## Min([`x1`, `x2`, `x3` ...]) {#min}

Returns the minimum value from the specified array of numbers.

```javascript
Min([15, 21])    // 15
```

## Mod(`x`, `y`) {#mod}

Returns the remainder of dividing `x` by `y`.

```javascript
Mod(8, 3)    // 2
```

## Mul(`x`, `y`) {#mul}

Returns the product of `x` and `y`.

```javascript
Mul(10, 1.5)    // 15
```

## Neg(`x`) {#neg}

Returns `x` with opposite sign.

```javascript
Neg(-5)    // 5
Neg(12)    // -12
```

## Not(`a`) {#not}

Returns logical negation of the `a`.

```javascript
Not(true)     // false
Not(5 < 1)    // true
```

## NotEq(`x`, `y`) {#noteq}

Returns false if `x` equal to `y` and true otherwise.

```javascript
NotEq(5, 5)           // false
NotEq(true, false)    // true
```

## NotGreater(`x`, `y`) {#notgreater}

Returns true if `x` is less than or equal to `y` and false otherwise.

```javascript
NotGreater(5, 5)    // true
NotGreater(6, 5)    // false
```

## NotSmaller(`x`, `y`) {#notsmaller}

Returns true if `x` is greater than or equal to `y`  and false otherwise.

```javascript
NotSmaller(5, 5)    // true
NotSmaller(5, 6)    // false
```

## Or(`a`, `b`) {#or}

Returns logical disjunction of boolean `a` and `b`.

```javascript
Or(true, false)        // true
Or(5 == 6, 20 < 10)    // false
```

## Percentile(`nums`, `percentage`) {#percentile}

Returns a value from `nums` below which a given `percentage` of values fall. `percentage` is in the range 0..1

```javascript
Percentile(([1, 2, 3, 4], 0.25))     // 2
Percentile(([1, 2, 3, 4], 0.75))     // 4
```

## Pow(`x`, `exponent`) {#pow}

Returns `x` to the power of `exponent`.

```javascript
Pow(2, 3)     // 8
Pow(2, -2)    // 0.25
```

## Qualifier(`x`) {#qualifier}

Extracts the qualifier from a qualified number `x`.

```javascript
Qualifier(Qnum(1.5, "="))   // =
Qualifier(Qnum(1.5, "<"))   // <
Qualifier(Qnum(1.5, ">"))   // >
```

## RandBetween(`m`, `n`) {#randbetween}

Returns a random integer number within the range from `m`, inclusive, to `n`, exclusive.

```javascript
RandBetween(5, 7)    // Randomly returns 5 or 6
```

## Rnd(`limit`) {#rnd}

Returns a random integer number within the range from 0, inclusive, to `n`, exclusive. The absolute value is taken if
the number is negative.

```javascript
Rnd(2)    // Randomly returns 0 or 1
```

## Round(`x`) {#round}

Returns the integer closest to `x`. Function rounds away from zero when there is no closest integer.

```javascript
Round(3.4)     // 3
Round(3.5)     // 4
Round(-3.5)    // -4
```

## RoundFloat(`x`, `decimalPlaces`) {#roundfloat}

Returns the number rounded up `x` to the number of decimal places specified by `decimalPlaces`.

`decimalPlaces` can be negative to round to even 10s, 100s, etc.

```javascript
RoundFloat(PI, 2)      // 3.14
RoundFloat(25, -1)     // 30
```

## Sin(`x`) {#sin}

Returns the sine of the `x`.

```javascript
Sin(0)         // 0
Sin(PI / 6)    // 0.5
```

## Smaller(`x`, `y`) {#smaller}

Returns true if `x` is less than `y` and false otherwise.

```javascript
Smaller(5, 5)    // false
Smaller(5, 6)    // true
```

## Sqrt(`x`) {#sqrt}

Returns the square root of the `x`.

```javascript
Sqrt(6.25)    // 2.5
```

## Sub(`x`, `y`) {#sub}

Returns the difference between `x` and `y`.

```javascript
Sub(10, 3)    // 7
```

## Tan(`x`) {#tan}

Returns the tangent of the `x`.

```javascript
Tan(0)         // 0
Tan(PI / 6)    // 0.5
```

## Xor(`a`, `b`) {#xor}

Returns logical exclusive disjunction of boolean `a` and `b`.

```javascript
Xor(true, false)        // true
Xor(5 == 5, 10 < 20)    // false
```
