<!-- TITLE: Tests: Split columns by RegExp -->
<!-- SUBTITLE: -->

# Tests: Split columns by regular expression

Matches the specified regular expression against the content of the specified column. Matched groups are added as new
columns.

## Testing scenarios

1. Split the column using different regular expressions (
   about [Regular expressions](http://www.regular-expressions.info))

* The result is shown in the dialog box "Results dialogs". After clicking the "OK" button new columns are created in the
  table

1. Change the prefix of new columns

* The names of the new columns will be "prefix+n" (n={1,2,3…})

1. Click on "History" button and return previous states

* ["Split column by RegExp"](extract-regexp.md) dialog return to previous state

1. Test non-functional modules (help, navigation, UI, UX)

* Non-functional modules work correctly and are intuitive

See also:

* [Split column by regular expression](extract-regexp.md)
* [About regular expressions](http://www.regular-expressions.info)
