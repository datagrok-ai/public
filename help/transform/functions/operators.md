<!-- TITLE: Operators -->
<!-- SUBTITLE: -->

# Operators

All operators have analogs among [functions](math-functions.md). In each case you can use the most appropriate option.

Template of using binary operators: `<Operand A> <operator> <Operand B>`

Template of using unary operators: `<operator> <Operand B>`

As operands of the operators, you can pass numeric scalars, numeric functions, [math constants](constants.md), boolean
scalars, boolean functions, or a column name. To pass a column cell, you can use the syntax `${columnName}`. Other ways
to use operands: `true`
, `false`, `PI`, `E` etc.

Operators and functions can organize expressions as complex as you like. You can also use parentheses to change the
standard sequence for evaluating operators. For example:

```javascript
Sin(PI / 6) * (17 - ${LENGTH}) < 9    // The result is a boolean value
```

*Operator List:*

> Here `A` is the left operand of the operator and the `B` is the right operand.

| Operator | Description                                                      | Similar Function                              |
|----------|------------------------------------------------------------------|-----------------------------------------------|
| `/`      | The result of dividing `A` by `B`                                | [Div(A, B)](math-functions.md#div)            |
| `*`      | The product of `A` and `B`                                       | [Mul(A, B)](math-functions.md#mul)            |
| `%`      | The remainder of dividing `A` by `B`                             | [Mod(A, B)](math-functions.md#mod)            |
| `^`      | Returns `A` to the power of `B`                                  | [Pow(A, B)](math-functions.md#pow)            |
| `+`      | The sum of two numbers `A` and `B`                               | [Add(A, B)](math-functions.md#add)            |
| `-`      | The difference between `A` and `B`                               | [Sub(A, B)](math-functions.md#sub)            |
| `==`     | True if `A` equal to `B` and false otherwise                     | [Eq(A, B)](math-functions.md#eq)              |
| `!=`     | False if `A` equal to `B` and true otherwise                     | [NotEq(A, B)](math-functions.md#noteq)        |
| `>`      | True if `A` is greater than `B` and false otherwise              | [Greater(A, B)](math-functions.md#greater)    |
| `<`      | True if `A` is less than `B` and false otherwise                 | [Smaller(A, B)](math-functions.md#smaller)    |
| `>=`     | True if `A` is greater than or equal to `B`  and false otherwise | [NotSmaller(A, B)](math-functions.md#notsmaller) |
| `<=`     | True if `A` is less than or equal to `B` and false otherwise     | [NotGreater(A, B)](math-functions.md#notgreater) |
| `and`    | Logical conjunction of boolean `A` and `B`                       | [And(A, B)](math-functions.md#and)            |
| `&&`     | Logical conjunction of boolean `A` and `B`                       | [And(A, B)](math-functions.md#and)            |
| `or`     | Logical disjunction of boolean `A` and `B`                       | [Or(A, B)](math-functions.md#or)              |
| `xor`    | Logical exclusive disjunction of boolean `A` and `B`             | [Xor(A, B)](math-functions.md#xor)            |
| `not`    | Logical negation of the `B`                                      | [Not(B)](math-functions.md#not)               |
| `!`      | Logical negation of the `B`                                      | [Not(B)](math-functions.md#not)               |
| `in`     | In operator, `A in [A,B]` returns `true`                         | In(A, B)             |
