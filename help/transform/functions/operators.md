<!-- TITLE: Operators -->
<!-- SUBTITLE: -->

# Operators

All operators have analogs among [functions](math-functions.md). In each case you can use the most appropriate option.

Template of using binary operators: `<Operand A> <operator> <Operand B>`

Template of using unary operators: `<operator> <Operand B>`

As operands of the operators, you can pass numeric scalars, numeric functions, [math constants](constants.md), boolean scalars, boolean functions, or a column name. To pass a column cell, you can use the syntax `${columnName}`. Other ways to use operands: `true`, `false`, `PI`, `E` etc.

Operators and functions can organize expressions as complex as you like. You can also use parentheses to change the standard sequence for evaluating operators. For example:

```javascript
Sin(PI / 6) * (17 - ${LENGTH}) < 9    // The result is a boolean value
```

*Operator List:*

| Operator | Description                                                              | Similar Function                                 |
| -------- | ------------------------------------------------------------------------ | ------------------------------------------------ |
| `/`      | Returns the result of dividing `A` by `B`                                | [Div(A, B)](math-functions.md#div)               |
| `*`      | Returns the product of `A` and `B`                                       | [Mul(A, B)](math-functions.md#mul)               |
| `'%'`    | Returns the remainder of dividing `A` by `B`                             | [Mod(A, B)](math-functions.md#mod)               |
| `^`      | Returns `A` to the power of `B`                                          | [Pow(A, B)](math-functions.md#pow)               |
| `+`      | Returns the sum of two numbers `A` and `B`                               | [Add(A, B)](math-functions.md#add)               |
| `-`      | Returns the difference between `A` and `B`                               | [Sub(A, B)](math-functions.md#sub)               |
| `==`     | Returns true if `A` equal to `B` and false otherwise                     | [Eq(A, B)](math-functions.md#eq)                 |
| `!=`     | Returns false if `A` equal to `B` and true otherwise                     | [NotEq(A, B)](math-functions.md#noteq)           |
| `>`      | Returns true if `A` is greater than `B` and false otherwise              | [Greater(A, B)](math-functions.md#greater)       |
| `<`      | Returns true if `A` is less than `B` and false otherwise                 | [Smaller(A, B)](math-functions.md#smaller)       |
| `>=`     | Returns true if `A` is greater than or equal to `B`  and false otherwise | [NotSmaller(A, B)](math-functions.md#notsmaller) |
| `<=`     | Returns true if `A` is less than or equal to `B` and false otherwise     | [NotGreater(A, B)](math-functions.md#notgreater) |
| `and`    | Returns logical conjunction of boolean `A` and `B`                       | [And(A, B)](math-functions.md#and)               |
| `&&`     | Returns logical conjunction of boolean `A` and `B`                       | [And(A, B)](math-functions.md#and)               |
| `or`     | Returns logical disjunction of boolean `A` and `B`                       | [Or(A, B)](math-functions.md#or)                 |
| `\|\|`   | Returns logical disjunction of boolean `A` and `B`                       | [Or(A, B)](math-functions.md#or)                 |
| `not`    | Returns logical negation of the `B`                                      | [Not(B)](math-functions.md#not)                  |
| `!`      | Returns logical negation of the `B`                                      | [Not(B)](math-functions.md#not)                  |
