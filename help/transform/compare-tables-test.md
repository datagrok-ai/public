<!-- TITLE: Tests: Compare tables -->
<!-- SUBTITLE: -->

# Tests: Compare tables

Instrument to compare the content of two tables

## Testing scenarios

1. Select two tables and key columns for comparison. Input parameters must be valid (key columns must have the same type
   and unique values)

* After clicking "OK" button comparison of tables is done. Created a new table with result. Check result for correct
  execution

1. Use "Value Columns" section to specify the columns to compare.

* Сolumns added for comparison

1. Use ["Сompare tables"](compare-tables.md) with non-valid input parameters (negative testing)

* "OK" button should be greyed out. Show a warning for the user with a reason.

1. Make sure different changes are picked up (changed values, deleted/inserted rows, deleted/inserted columns, etc)

1. Click on "History" button and return previous states

* "Compare Tables" dialog return to previous state

1. Test non-functional modules (popup menu, help, navigation, properties, etc.)

* Non-functional modules work correctly and are intuitive

See also:

* [Сompare tables](compare-tables.md)
* [Join tables test](../transform/join-tables-test.md)
