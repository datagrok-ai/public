<!-- TITLE: Add New Column -->
<!-- SUBTITLE: -->

# Adding New Columns

Adds a column of the specified type to the current table, and initializes it using the specified
expression (mathematical function, constants, platform objects properties and functions).

To reference a column, specify its name in the curly brackets, preceded by the dollar sign, like
that: ${WIDTH}. While editing the formula, press '$' to opens up a column list popup; use arrows and
Enter to select it.

For formulas where row index is required, 'row' variable is available.

Examples:
```
1.57 * ${WEIGHT} / ${HEIGHT}
-log(${IC50} * LN10)
```

| Constants | Description                             |
|-----------|-----------------------------------------|
| E         | Euler's number (approx. 2.718)          |
| LN2       | Natural logarithm of 2 (approx. 0.693)  |
| LN10      | Natural logarithm of 10 (approx. 2.302) |
| LOG2E     | Base-2 logarithm of E (approx. 1.442)   |
| LOG10E    | Base-10 logarithm of E (approx. 0.434)  |
| PI        | PI (approx. 3.14)                       |
| SQRT1_2   | Square root of 1/2 (approx. 0.707)      |
| SQRT2     | Square root of 2 (approx. 1.414)        |


| Functions   | Description                                                       |
|-------------|-------------------------------------------------------------------|
| Abs(x)      | Absolute value of x                                               |
| Acos(x)     | Arccosine of x, in radians                                        |
| Asin(x)     | Arcsine of x, in radians                                          |
| Atan(x)     | Arctangent of x as a numeric value between -PI/2 and PI/2 radians |
| Atan2(y, x) | Arctangent of the quotient of its arguments                       |
| Ceil(x)     | x, rounded upwards to the nearest integer                         |
| Cos(x)      | Cosine of x (x is in radians)                                     |
| Exp(x)      | E^x                                                               |
| Floor(x)    | x, rounded downwards to the nearest integer                       |
| Log(x)      | Natural logarithm (base E) of x                                   |
| Max(x, y)   | Number with the highest value                                     |
| Min(x, y)   | Number with the lowest value                                      |
| Pow(x, y)   | Value of x to the power of y                                      |
| Random()    | A floating-point, pseudo-random number 0..1                       |
| Round(x)    | Rounds x to the nearest integer                                   |
| Sin(x)      | Returns the sine of x (x is in radians)                           |
| Sqrt(x)     | Square root of x                                                  |
| Tan(x)      | Tangent of an angle                                               |


To treat data as strings use quotes, for example: 
```
"a" + "b"
```
results in "ab".

Detection of datetime for most datetime formats works automatically. 
But in some cases it is hard to recognize is it formula or datetime.
For example to get "1/1/2019" date, to solve this set datetime type and use quotes:
```
"1/1/2019"
```

### Videos

<iframe width="560" height="315" src="https://www.youtube.com/embed/-yTTaS_WOU4" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

See also:

  * [Grok Scripting](../features/grok-script.md)
  * [Function](../entities/function.md)
  * [Column selectors](../viewers/column-selectors.md)
