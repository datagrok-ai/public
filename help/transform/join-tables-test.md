<!-- TITLE: Tests: Join tables -->
<!-- SUBTITLE: -->

# Tests: Join tables

Joins two tables, using the specified key columns. Сan be used inner join, outer join, left and right joins

## Testing scenarios

1. Select two tables and key columns for join. Use different kinds of join.

    * Created a new table with the result after clicking the "OK" button
    * Result corresponds to the selected join type
    * Check the result for correct execution depending on the selected type of join

2. Use 'Value Columns' section to specify the columns to join.

3. Use ["Join tables"](join-tables.md) with non-valid input parameters (different types of key columns, not a unique key
   column)

    * "OK" button should be greyed out. Show a warning for the user with a reason

4. Click on "History" button and return previous states

    * ["Join tables"](join-tables.md) dialog return to previous state

5. Test non-functional modules (preview window, popup menu, help, navigation, etc.)

    * Non-functional modules work correctly and are intuitive

See also:

* [Join tables](join-tables.md)
* [Compare tables test](../explore/compare-tables-test.md)
