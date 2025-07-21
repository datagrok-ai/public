---
title: "Add new column"
---

Adds a column of the specified type to the current table, and initializes it using the specified expression (
mathematical function, constants, platform objects properties and functions).

To add new columns, click on the `Add New Column` icon on the toolbar or go to top menu `Edit` -> `Add New Column`. Key features are:

* Support for functions implemented with Python, R, Julia, JavaScript, C++, and others. To add function to an editor, type it manually, drag and drop from
functions registry on the right or use *plus* icon. You can combine functions written in different languages in one formula.
![add function](./add_function_to_editor.gif)
* Auto-suggested functions based on the column type and semantic type. To use:
  * set functions sorting type in function registry to 'By relevance'
  * select column of interest
  * the functions in functions registry are sorted automatically. More relevant functions are on top of the list.
  * drag and drop function to the editor field or click *plus* icon. Corresponding parameter is prefilled automatically with selected column.
![functions suggestions](./add_new_column_functions_suggestions.gif)
* Interactive preview of results as you type
* Autocompletion for functions (including packages names) and columns. Suggestions appear as you type.
* Different highlights within formula for better readability. For instance, column names are highlighted with in bold blue font.
* Validation against various types of mistakes including syntax errors, missing columns detection, incorrect data types, unmatching brackets.
* Resulting column type autodetection
* Fast function and column search
* History, saving and reusing formulas

Adding columns to formulas:

* scalar functions
  * To reference each row of a column, specify its name in the curly brackets, preceded by the dollar sign: `${Width}`.
  For example you can use this expression in function like that: `Round(${Width})`.
  * To reference a whole column, specify its name in the square brackets, preceded by the dollar sign: `$[Width]`. For
  example you can use this expression in function like that: `Avg($[Width])`.

* vector function
  * To reference a whole column, specify its name in the curly brackets, preceded by the dollar sign: `${molecule}`. For example you can use this expression in function like that: `Chem:getInchis(${molecule})`.

To add a column to a formula, drag it to the editor. Alternatively, use the keyboard:

1. Open a column list popup by pressing '**$**'.
1. Select the column you want using the **up** and **down arrows**, then press **Enter**.

For formulas where row index is required, `row` variable is available.

Example:

```javascript
1.57 * RoundFloat(${Weight}, 2) / Avg($[Weight]) - log(${IC50} * PI)
```

To treat data as strings use quotes, for example:

```javascript
"Police" + "man"    // "Policeman"
```

The platform supports a large number of functions, constants and operators. You can find out about them in the
corresponding sections of the help system:

- [Binning functions](functions/binning-functions.md)
- [Constants](functions/constants.md)
- [Conversion functions](functions/conversion-functions.md)
- [DateTime functions](functions/datetime-functions.md)
- [Math functions](functions/math-functions.md)
- [Operators](functions/operators.md)
- [Stats functions](functions/stats-functions.md)
- [Text functions](functions/text-functions.md)
- [TimeSpan functions](functions/timespan-functions.md)

## Videos

[![Add New Columns](../uploads/youtube/add_new_columns.png "Open on Youtube")](https://www.youtube.com/watch?v=-yTTaS_WOU4)

See also:

- [Scripting](../develop/under-the-hood/grok-script.md)
- [Function](../datagrok/concepts/functions/functions.md)
- [Column selectors](../visualize/viewers/column-selectors.md)
