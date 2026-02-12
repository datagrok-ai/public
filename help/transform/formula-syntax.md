---
title: Formula syntax
---

```mdx-code-block
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';
```

# Formula syntax

Datagrok uses a unified formula syntax across the platform for calculated columns, filters, and other expressions.  

## Referencing data in formulas

### Current row value

Use `${ColumnName}` to reference the value in the current row.  
This is the most common case, when a formula is evaluated independently for each row.  

The same syntax is used when passing a column to functions that process all values at once and return a result for each row, such as formula equalities.

```javascript
${AGE} > 18
${PRICE} * ${QUANTITY}
${Chemical Space Y} = ln(${Chemical Space X}) + 5  // formula line equality
```

### Entire column

Use `$[ColumnName]` to reference all values in a column at once. This is used when the formula needs column-level information, such as averages, totals, or other aggregated values.

```javascript
${IC50} / Avg($[IC50])   
${VALUE} > Median($[VALUE])
```

### Row index 

Use `row` when the formula explicitly requires the row index.

```javascript
row % 2 == 0
```

## Operators, constants, and literals

### Operators

Datagrok supports arithmetic, comparison, logical, and membership operators.  
You can use them with numbers, logical values (`true`, `false`), [constants](#constants) like `PI` and `E`, function results, or current row values `${ColumnName}`.  

You can combine operators and functions to build complex expressions.  
Use parentheses `()` to control the order of evaluation.

```javascript
Sin(PI / 6) * (17 - ${LENGTH}) < 9    // The result is a boolean value
```

Each operator also has a corresponding [function](functions/math-functions.md), so you can choose whichever form is more convenient.
<details> 
<summary> Reference: list of operators </summary> 
<div>

 Here `A` is the left operand of the operator and the `B` is the right operand.

| Operator | Description                                                      | Similar Function                              |
|----------|------------------------------------------------------------------|-----------------------------------------------|
| `/`      | The result of dividing `A` by `B`                                | [Div(A, B)](functions/math-functions.md#div)            |
| `*`      | The product of `A` and `B`                                       | [Mul(A, B)](functions/math-functions.md#mul)            |
| `%`      | The remainder of dividing `A` by `B`                             | [Mod(A, B)](functions/math-functions.md#mod)            |
| `^`      | Returns `A` to the power of `B`                                  | [Pow(A, B)](functions/math-functions.md#pow)            |
| `+`      | The sum of two numbers `A` and `B`                               | [Add(A, B)](functions/math-functions.md#add)            |
| `-`      | The difference between `A` and `B`                               | [Sub(A, B)](functions/math-functions.md#sub)            |
| `==`     | True if `A` equal to `B` and false otherwise                     | [Eq(A, B)](functions/math-functions.md#eq)              |
| `!=`     | False if `A` equal to `B` and true otherwise                     | [NotEq(A, B)](functions/math-functions.md#noteq)        |
| `>`      | True if `A` is greater than `B` and false otherwise              | [Greater(A, B)](functions/math-functions.md#greater)    |
| `<`      | True if `A` is less than `B` and false otherwise                 | [Smaller(A, B)](functions/math-functions.md#smaller)    |
| `>=`     | True if `A` is greater than or equal to `B`  and false otherwise | [NotSmaller(A, B)](functions/math-functions.md#notsmaller) |
| `\<=`     | True if `A` is less than or equal to `B` and false otherwise     | [NotGreater(A, B)](functions/math-functions.md#notgreater) |
| `and`    | Logical conjunction of boolean `A` and `B`                       | [And(A, B)](functions/math-functions.md#and)            |
| `&&`     | Logical conjunction of boolean `A` and `B`                       | [And(A, B)](functions/math-functions.md#and)            |
| `or`     | Logical disjunction of boolean `A` and `B`                       | [Or(A, B)](functions/math-functions.md#or)              |
| `xor`    | Logical exclusive disjunction of boolean `A` and `B`             | [Xor(A, B)](functions/math-functions.md#xor)            |
| `not`    | Logical negation of the `B`                                      | [Not(B)](functions/math-functions.md#not)               |
| `!`      | Logical negation of the `B`                                      | [Not(B)](functions/math-functions.md#not)               |
| `in`     | In operator, `A in [A,B]` returns `true`                         | In(A, B)             |

</div> 
</details>

### Constants

Datagrok provides predefined constants that you can use in formulas by simply writing their names.

```javascript
E * ln(${VALUE}) / SQRT2
```

<details> 
<summary> Reference: list of constants </summary> 

<Tabs>
<TabItem value="Common" label="Common" default>

| Name      | Description                                                                                        | Value              |
| --------- | -------------------------------------------------------------------------------------------------- | ------------------ |
| `E`       | [Euler's number](https://en.wikipedia.org/wiki/E_(mathematical_constant))                          | 2.718281828459045  |
| `LN2`     | [Natural logarithm of 2](https://en.wikipedia.org/wiki/Natural_logarithm_of_2)                     | 0.6931471805599453 |
| `LN10`    | [Natural logarithm of 10](https://en.wikipedia.org/wiki/Natural_logarithm#Natural_logarithm_of_10) | 2.302585092994046  |
| `LOG2E`   | Base-2 logarithm of E                                                                              | 1.4426950408889634 |
| `LOG10E`  | Base-10 logarithm of E                                                                             | 0.4342944819032518 |
| `PI`      | [Pi](https://en.wikipedia.org/wiki/Pi)                                                           | 3.1415926535897932 |
| `SQRT1_2` | Square root of 1/2                                                                                 | 0.7071067811865476 |
| `SQRT2`   | [Square root of 2](https://en.wikipedia.org/wiki/Square_root_of_2)                                 | 1.4142135623730951 |

</TabItem>
<TabItem value="mathematical" label="Other mathematical">

| Name      | Description                                                                                                      | Value              |
| --------- | ---------------------------------------------------------------------------------------------------------------- | ------------------ |
| `APERY`   | [Apéry's constant](https://en.wikipedia.org/wiki/Ap%C3%A9ry%27s_constant)                                        | 1.202056903159594  |
| `BACKH`   | [Backhouse's constant](https://en.wikipedia.org/wiki/Backhouse%27s_constant)                                     | 1.456074948582689  |
| `BERN`    | [Bernstein's constant](https://en.wikipedia.org/wiki/Bernstein%27s_constant)                                     | 0.280169499023869  |
| `CATALAN` | [Catalan's constant](https://en.wikipedia.org/wiki/Catalan%27s_constant)                                         | 0.915965594177219  |
| `CONWAY`  | [Conway's constant](https://en.wikipedia.org/wiki/Look-and-say_sequence#Growth_in_length)                        | 1.303577269034296  |
| `EUMA`    | [Euler–Mascheroni constant](https://en.wikipedia.org/wiki/Euler%E2%80%93Mascheroni_constant)                     | 0.577215664901532  |
| `ERDBOR`  | [Erdős–Borwein constant](https://en.wikipedia.org/wiki/Erd%C5%91s%E2%80%93Borwein_constant)                      | 1.606695152415291  |
| `FEIG1`   | [Feigenbaum constant 1](https://en.wikipedia.org/wiki/Feigenbaum_constants)                                      | 4.669201609102990  |
| `FEIG2`   | [Feigenbaum constant 2](https://en.wikipedia.org/wiki/Feigenbaum_constants)                                      | 2.502907875095892  |
| `FHL`     | [First Hardy–Littlewood conjecture](https://en.wikipedia.org/wiki/Twin_prime#Conjectures)                        | 0.660161815846869  |
| `FROB`    | [Fransén–Robinson constant](https://en.wikipedia.org/wiki/Frans%C3%A9n%E2%80%93Robinson_constant)                | 2.807770242028519  |
| `GAKUWI`  | [Gauss–Kuzmin–Wirsing operator](https://en.wikipedia.org/wiki/Gauss%E2%80%93Kuzmin%E2%80%93Wirsing_operator)     | 0.303663002898732  |
| `GODI`    | [Golomb–Dickman constant](https://en.wikipedia.org/wiki/Golomb%E2%80%93Dickman_constant)                         | 0.624329988543550  |
| `GRATIO`  | [Golden ratio](https://en.wikipedia.org/wiki/Golden_ratio)                                                       | 1.6180339887498948 |
| `HASAMC`  | [Hafner–Sarnak–McCurley constant](https://en.wikipedia.org/wiki/Hafner%E2%80%93Sarnak%E2%80%93McCurley_constant) | 0.353236371854995  |
| `KHIN`    | [Khinchin's constant](https://en.wikipedia.org/wiki/Khinchin%27s_constant)                                       | 2.685452001065306  |
| `LANRAM`  | [Landau–Ramanujan constant](https://en.wikipedia.org/wiki/Landau%E2%80%93Ramanujan_constant)                     | 0.764223653589220  |
| `LEVY`    | [Lévy's constant](https://en.wikipedia.org/wiki/L%C3%A9vy%27s_constant)                                          | 3.275822918721811  |
| `MEME`    | [Meissel–Mertens constant](https://en.wikipedia.org/wiki/Meissel%E2%80%93Mertens_constant)                       | 0.261497212847642  |
| `MILLS`   | [Mills' constant](https://en.wikipedia.org/wiki/Mills%27_constant)                                               | 1.306377883863080  |
| `NIVEN`   | [Niven's constant](https://en.wikipedia.org/wiki/Niven%27s_constant)                                             | 1.705211140105367  |
| `OMEGA`   | [Omega constant](https://en.wikipedia.org/wiki/Omega_constant)                                                   | 0.567143290409783  |
| `PLASTIC` | [Plastic constant](https://en.wikipedia.org/wiki/Plastic_number)                                                 | 1.324717957244746  |
| `RAMSOL`  | [Ramanujan–Soldner constant](https://en.wikipedia.org/wiki/Ramanujan%E2%80%93Soldner_constant)                   | 1.451369234883381  |
| `RECFIB`  | [Reciprocal Fibonacci constant](https://en.wikipedia.org/wiki/Reciprocal_Fibonacci_constant)                     | 3.359885666243177  |
| `SIERP`   | [Sierpiński's constant](https://en.wikipedia.org/wiki/Sierpi%C5%84ski%27s_constant)                              | 2.584981759579253  |
| `SQRT3`   | [Square root of 3](https://en.wikipedia.org/wiki/Square_root_of_3)                                               | 1.732050807568877  |
| `SQRT5`   | [Square root of 5](https://en.wikipedia.org/wiki/Square_root_of_5)                                               | 2.236067977499789  |
| `UNIPAR`  | [Universal parabolic constant](https://en.wikipedia.org/wiki/Universal_parabolic_constant)                       | 2.295587149392638  |

</TabItem>
<TabItem value="Chemical and physical" label="Chemical and physical">

| Name       | Description                                                             | Value                                                    |
| ---------- | ----------------------------------------------------------------------- | -------------------------------------------------------- |
| `AVOGADRO` | [Avogadro's Number](https://en.wikipedia.org/wiki/Avogadro_constant)    | 6.02214e+23                                              |
| `FARADAY`  | [Faraday Constant](https://en.wikipedia.org/wiki/Faraday_constant)      | 96485.33 (C·mol<sup>-1</sup>)                            |
| `ATOM`     | [Atomic Mass Constant](https://en.wikipedia.org/wiki/Dalton_(unit))     | 1.66053906660e-27 (kg)                                   |
| `GAS`      | [Molar Gas Constant](https://en.wikipedia.org/wiki/Gas_constant)        | 8.31446261815324 (m3⋅Pa⋅K<sup>−1</sup>⋅mol<sup>−1</sup>) |
| `COULOMB`  | [Coulomb constant](https://en.wikipedia.org/wiki/Coulomb_constant)      | 8.9875517923e+9 (N·m<sup>2</sup>/C<sup>2</sup>)          |
| `MAXSPEED` | [Speed of Light (Vacuum)](https://en.wikipedia.org/wiki/Speed_of_light) | 299792458 (m/s)                                          |
| `BOLTZ`    | [Boltzmann constant](https://en.wikipedia.org/wiki/Boltzmann_constant)  | 1.380649e-23 (J⋅K<sup>−1</sup>)                          |
| `ECHARGE`  | [Elementary charge](https://en.wikipedia.org/wiki/Elementary_charge)    | 1.602176634e-19 (C)                                      |
| `GRAVITY`  | [Standard gravity](https://en.wikipedia.org/wiki/Standard_gravity)      | 9.80665 (m/s<sup>2</sup>)                                |
| `PLANCK`   | [Planck constant](https://en.wikipedia.org/wiki/Planck_constant)        | 6.62607004e-34 (J⋅Hz<sup>−1</sup>)                       |
| `EMASS`    | [Mass of electron](https://en.wikipedia.org/wiki/Electron_rest_mass)    | 9.1093837015e-31 (kg)                                    |
| `NMASS`    | [Mass of neutron](https://en.wikipedia.org/wiki/Neutron)                | 1.67492749804e-27 (kg)                                   |
| `PMASS`    | [Mass of proton](https://en.wikipedia.org/wiki/Proton)                  | 1.67262192369e-27 (kg)                                   |

</TabItem>
</Tabs>
 
</details>


### String literals

String values can be written using single or double quotes. Both forms are supported and interchangeable.

```javascript
'Medium'
"Medium"
```

### Whitespace and line breaks

Spaces around operators are optional, and formulas can include line breaks. 
Whitespace and formatting do not affect how formulas are evaluated.

The following examples are valid and equivalent:

```javascript
${A}+${B}

${A} + ${B}

${A}
+
${B}
```

```javascript
in(20,60,40)
in(20, 60, 40)
```

## Functions

### Built-in functions

Datagrok provides built-in functions, including [calculations](functions/math-functions.md), [text manipulation](functions/text-functions.md), [statistical analysis](functions/stats-functions.md), [date/time](functions/datetime-functions.md) operations, [conversions](functions/conversion-functions.md), [binning](functions/binning-functions.md), and [timespan](functions/timespan-functions.md) operations. 

```javascript
RoundFloat(${IC50} / Median($[IC50]) * E, 3)
```

### Custom functions

Datagrok allows you to create custom functions in Python, R, Julia, JavaScript, C++, and other languages, as long as they are properly annotated.   
See [Function annotations](../datagrok/concepts/functions/func-params-annotation.md) for details.  

Functions can be written as standalone scripts or included in packages.  

#### Package function

Use `PackageName:FunctionName` to call package function. 

```javascript
Chem:getInchis(${Structure})
```

#### User-defined function

Use `UserLogin:FunctionName` to call user-defined function.

The following example shows a custom function written in Python:

```python
#name: Len
#language: python
#input: string s
#output: int length

length = len(s)
```

Once the script is saved, the function can be called in formulas like this:

```javascript
UserName:Len(${ColumnName})
```

## Calculated columns

Calculated columns generate new columns based on formulas. They can reference other calculated columns and update automatically when the source data changes.

The **Add New Column** dialog supports creating calculated columns with real-time preview, autocomplete, and formula validation, including checks for syntax errors, missing columns, and type mismatches.   
See [Add New Column](add-new-column.md) for details.

In-viewer filter expressions can be edited in the **Edit Formula** dialog, which provides the same enhanced functionality as **Add New Column**.

### Complex calculated columns (multiple outputs)

In some cases, a single formula can produce multiple columns at once. This is useful when a function returns several related results that should stay synchronized with the source data. When used from **Add New Column**, all resulting columns are created automatically and remain synchronized with the source data.   
For details, see [Сomplex calculated columns](../datagrok/concepts/functions/func-params-annotation.md#complex-calculated-columns).

<details>
<summary>See visual</summary>
<Tabs>
<TabItem value="Calculated column" label="Calculated column" default>

</TabItem>
<TabItem value="Complex" label="Complex calculated columns" default>

</TabItem>
<TabItem value="Filter" label="Filter expression" default>
![In-viewer filter formula example](in-viewer-filter.gif)
</TabItem>
</Tabs>
</details>

 
## See also

* [Functions](../datagrok/concepts/functions/functions.md)
* [Data wrangling](../transform/transform.md)
* [Macros](../datagrok/navigation/panels/panels.md#console)