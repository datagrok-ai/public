<!-- TITLE: Tests: Split Columns -->
<!-- SUBTITLE: -->

# Tests: Split columns

Splits textual data from the specified column into other columns, based on the specified delimiter.

## Testing scenarios

1. Split the column using different types of delimiters (numbers, letters, special characters, space)
   * The result is shown in the dialog box "Results dialogs". After clicking the "OK" button new columns are created in the table

1. Change the prefix of new columns
   * The names of the new columns will be "prefix+n" (n={1,2,3…})

1. Enter more than one character to "Split by" field. Mark "Consecutive as one" 
   * Embedded characters are treated as one delimiter

1. Enter more than one character to "Split by" field. Mark "By chars" 
   * Splitting is carried out for each entered character 

1. Click on "History" button and return previous states
["Split column"](../dialogs/text-to-columns.md) dialog return to previous state

1. Test non-functional modules (help, navigation, UI, UX)
   * Non-functional modules work correctly and are intuitive

See also:
 * [Split column](../dialogs/text-to-columns.md)
 * [Split column by RegExp test](../tests/split-column-regexp-test.md)
 * [Split columns Auto Test](text-to-columns-test.side)
