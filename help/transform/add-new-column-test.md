<!-- TITLE: Tests: Add new column -->
<!-- SUBTITLE: -->

# Tests: Add new column

[Adds a column](add-new-column.md) of the specified type to the current table, and initializes it using the specified
expression (mathematical function, constants, platform objects properties and functions).

## Testing scenarios

1. Open *"demog"* table

1. Open [*"Add new column"*](add-new-column.md) dialog by clicking corresponding icon on toolbar

1. Enter *"test_column"* to *"Name"* field

* Column name changed to *"test_column"* in preview

1. Enter `1` to *"Formula"* field

* Type of new column is automatically determined as *int*
* All cells of new column were filled with value "1" in preview

1. Enter `1.5` to *"Formula"* field

* Type of new column is automatically determined as *double*
* All cells of new column were filled with value "1.5" in preview

1. Enter `test` to *"Formula"* field

* Type of new column is automatically determined as *string*
* All cells of new column were filled with value "test" in preview

1. Enter `true` to *"Formula"* field

* Type of new column is automatically determined as *bool*
* All cells of new column were filled with value "true" in preview

1. Enter `1+1.5` to *"Formula"* field

* Type of new column is automatically determined as *double*
* All cells of new column were filled with value "2.5" in preview

1. Forcibly change type of new column from *double* to *int* in *"Type"* field. (`1+1.5` entered in *"Formula"*
   field)

* New column type changed to *int*
* All cells values of new column rounded to 3 in preview

1. Change type of new column to *string* to *"Type"* field. (`1+1.5` entered in *"Formula"*
   field)

* New column type changed to *string*
* All cells of new column were filled with value "2.5" (*string*) in preview

1. Mark  *"Treat as string"* as true

* All cells of new column were filled with value "1+1.5" (*string*) in preview

1. Re-open [*"Add new column"*](add-new-column.md) dialog

1. Add "Height" column to *"Formula"* field. Enter ```$``` to call drop-down list of columns or drag column from
   anywhere in platform ([Column Manager](../explore/column-manager.md)
   , [Table](../datagrok/table.md),
   [Properties panel](../datagrok/navigation.md#properties))

* Type of new column was automatically defined as *double* (*"Height"* column type)
* Values of new column are equal to values of *"Height"* column in preview
* *"Name"* field is automatically filled by entered in *"Formula"* field

1. Enter `$HEIGHT+${WEIGHT}` to *"Formula"* field

* Each row of new column is equal to sum of corresponding rows of *"Height"* and *"Weight"*
  columns

1. Add constant value to *"Formula"* field (PI, E, LN2, LN10, etc.)

* Values of new column are calculated depending on values of entered constants

1. Use math functions in  *"Formula"* field (Abs(), Acos(), Min(), Max(), Log(), etc.)

* Values of new column are calculated depending on entered funtions

1. Add columns of different types. Use different combinations of column types and entered values in formula.

1. Remove some values from each type of column. Use columns containing empty cells in *"Formula"*
   field

* For math operations with empty cells, they must remain empty in new column

1. Enter non-correct formula in appropriate fields

* Preview is empty
* `OK` not active for clicking

1. Register new function `norm2`(*) using JavaScript (**Tools | Scripting | JavaScript**)

1.

Enter `(norm2(Demo:DemoScripts:BMI($HEIGHT, $WEIGHT), $AGE)+norm2(Demo:DemoScripts:BSA($HEIGHT, $WEIGHT), $AGE))*100`
to *"Formula"* field and execute dialog

* New column added to table witch result of entered formula
* Added column has the name "New Column" (when you enter long formula, name is not auto complete)

1. Test non-functional modules (UI, help, navigation, console, history)

See also:

* [Add new column Auto Test](add-new-column-test.side)

(*):

```
grok.functions.register({
       signature: 'double norm2(double x, double y)',
       run: (x, y) => Math.sqrt(x / y)});
```
