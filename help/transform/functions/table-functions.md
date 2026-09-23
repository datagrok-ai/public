---
title: "Table and column functions"
description: Reference for formula functions that reference tables and columns by name, look values up in other tables, and compute running values over a column.
keywords:
  - lookup function
  - vlookup formula
  - running total
  - cumulative sum
  - moving average
  - row count
---

These functions work with tables and columns rather than with a single value. The table name is optional where shown
in brackets and defaults to the table the formula runs on. In the **Add New Column** dialog, the open tables and the
table's columns are offered as you type each argument.

*Function List:*

- [Column](#column)
- [CumSum](#cumsum)
- [IndexOf](#indexof)
- [Lookup](#lookup)
- [MovingAvg](#movingavg)
- [RowCount](#rowcount)
- [Table](#table)
- [Value](#value)

## Column(`columnName`, [`tableName`]) {#column}

Returns a column by name. Use it where a whole column is expected, including a column of another open table.

```javascript
Avg(Column("Width"))                  // Same as Avg($[Width])
Avg(Column("price", "products"))      // Average price in the "products" table
Column("Width").stats.max             // Column properties are available too
```

## CumSum(`column`, [`orderBy`], [`descending`], [`by`]) {#cumsum}

Returns the running total of a column. Empty values are skipped and stay empty. Reference the columns as `${name}`.

The optional arguments are passed by name:

| Argument     | Meaning                                                                                    |
|--------------|--------------------------------------------------------------------------------------------|
| `orderBy`    | A column to accumulate along, instead of the row order. Equal values keep their row order. |
| `descending` | `true` to walk `orderBy` from the largest value to the smallest.                           |
| `by`         | A column whose groups get their own running total.                                         |

```javascript
CumSum(${amount})                                   // 10, 25, (empty), 60 for 10, 15, (empty), 35
CumSum(${amount}) / Sum($[amount])                  // Cumulative share of the total
CumSum(${amount}, orderBy=${date})                  // Running total by date, whatever the row order
CumSum(${amount}, orderBy=${date}, by=${region})    // The same, separately for each region
CumSum(${amount}, orderBy=${date}, descending=true) // From the latest date back
```

## IndexOf(`column`, `value`) {#indexof}

Returns the zero-based index of the first row of the column that holds the value, or -1 when there is none.
It pairs with [Value](#value), which takes the same zero-based row.

```javascript
IndexOf(${id}, "A-17")                // Row of the first "A-17"
IndexOf(${id}, ${id}) == row - 1      // True for the first occurrence of each id
```

## Lookup(`tableName`, `keyColumn`, `key`, `valueColumn`) {#lookup}

Finds the first row of a table where the key column equals the key, and returns the value column of that row.
Returns an empty value when there is no match. This is the equivalent of a spreadsheet VLOOKUP.

```javascript
Lookup("products", "id", ${productId}, "price")              // Price of each row's product
${quantity} * Lookup("products", "id", ${productId}, "price")
```

## MovingAvg(`column`, `window`, [`orderBy`], [`descending`], [`by`], [`minPeriods`]) {#movingavg}

Returns the average of the current and the preceding rows, `window` rows in total. Empty values are ignored, and the
first rows average the rows available so far. Reference the columns as `${name}`.

`orderBy`, `descending`, and `by` work as in [CumSum](#cumsum). With `minPeriods`, the result stays empty until the
window holds that many values.

```javascript
MovingAvg(${price}, 7)                              // 7-row trailing average
MovingAvg(${price}, 7, orderBy=${date})             // Along the date column
MovingAvg(${price}, 7, orderBy=${date}, by=${ticker}, minPeriods=7)   // Per ticker, full windows only
```

## RowCount([`tableName`]) {#rowcount}

Returns the number of rows in a table.

```javascript
RowCount()                            // Rows in the current table
row / RowCount()                      // Relative position of the row
RowCount("products")
```

## Table([`tableName`]) {#table}

Returns a table by name, or the current table. Table properties are available on the result.

```javascript
Table("products").rowCount
Table().name
```

## Value(`columnName`, [`row`], [`tableName`]) {#value}

Returns the value of a column in the given row. The row is zero-based and defaults to the current row.

```javascript
Value("price", 0)                     // Price in the first row
${price} - Value("price", row - 2)    // Difference with the previous row (`row` is one-based)
Value("rate", 0, "settings")          // A value from another table
```
