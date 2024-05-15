---
title: "Statistical functions"
---

This type of function assumes working with whole columns or a set of numbers.

As parameters of the function, you can pass a whole column using the square brackets syntax `$[columnName]` or a list of
any numbers enclosed in square brackets using the syntax `[n1, n2, ....]`. Null arguments are ignored and do not affect
the results of statistical functions.

You can also use expressions such as `[${Width}, 270, ${Height}]` to pass row values as numbers in a list.

> By default, real numbers in a new column show only 2 digits after the integer part.
> You can change this behavior by setting the appropriate format for displaying the data in the column.
> To do this, click on the column heading and in the opened context menu select "Format" option.

*Function List:*

- [Avg](#avg)
- [Kurt](#kurt)
- [Max](#max)
- [Med](#med)
- [Min](#min)
- [MissingValueCount](#missingvaluecount)
- [Q1](#q1)
- [Q2](#q2)
- [Q3](#q3)
- [Skew](#skew)
- [StDev](#stdev)
- [Sum](#sum)
- [TotalCount](#totalcount)
- [ValueCount](#valuecount)
- [Variance](#variance)

## Avg(`data`) {#avg}

Returns the average (arithmetic mean) of the data.

```javascript
Avg([1, 2, 3, 4])     // 2.5
Avg($[Age])           // Average age
```

## Kurt(`data`) {#kurt}

Returns the kurtosis of the data. Kurtosis characterizes the relative peakedness or flatness of a distribution compared
with the normal distribution.

```javascript
Kurt([1, 2, 3])      // -1.5
Kurt($[Readings])
```

## Max(`data`) {#max}

Returns the maximum value of the data.

```javascript
Max([1.5, -2, 1.9])    // 1.9
Max($[Length])         // Maximum length
```

## Med(`data`) {#med}

Returns the median of the data.

```javascript
Med([1, 2, 3])         // 2
Med([$[Dev]])
```

## Min(`data`) {#min}

Returns the minimum value of the data.

```javascript
Min([1.5, -2, 1.9])    // -2
Min($[Length])         // Minimum length
```

## MissingValueCount(`data`) {#missingvaluecount}

Returns the number of empty (null) elements.

```javascript
MissingValueCount([10, null, 7])    // 1
MissingValueCount($[Weight])
```

## Q1(`data`) {#q1}

Returns the first quartile of the data.

```javascript
Q1([7, 2, -3, 4])    // 2
Q1($[Value])
```

## Q2(`data`) {#q2}

Returns the second quartile of the data.

```javascript
Q2([7, 2, -3, 4])    // 3
Q2($[Value])
```

## Q3(`data`) {#q3}

Returns the third quartile of the data.

```javascript
Q3([7, 2, -3, 4])   // 7
Q3($[Value])
```

## Skew(`data`) {#skew}

Returns the skewness of the data.

```javascript
Skew([1, 2, 3])         // 0
Skew($[Indications])
```

## StDev(`data`) {#stdev}

Returns the standard deviation of the data.

```javascript
StDev([7, 14, 21])    // 7
StDev($[Weight])
```

## Sum(`data`) {#sum}

Returns the sum of elements.

```javascript
Sum([-1, 4, 12, 5])    // 20
Sum($[Price])    //
```

## TotalCount(`data`) {#totalcount}

Returns the total number of elements.

```javascript
TotalCount([8, null, 1])    // 3
TotalCount($[Weight])
```

## ValueCount(`data`) {#valuecount}

Returns the number of non-null elements.

```javascript
ValueCount([8, null, 1])    // 2
ValueCount($[Weight])
```

## Variance(`data`) {#variance}

Returns the variance of elements.

```javascript
Variance([-5, 1])    // 18
Variance($[Size])
```
