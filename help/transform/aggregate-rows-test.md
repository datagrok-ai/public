<!-- TITLE: Tests: Aggregate rows -->
<!-- SUBTITLE: -->

# Test: Aggregate rows

["Aggregate rows"](aggregate-rows.md) tool allows to interactivity define aggregation logic, and immediately see results

## Testing scenarios

1. Add variables to "Columns", "Rows" fields by using button "+" or drag-and-drop columns. Select different aggregation
   metrics in "Measures" field. Try different types of variables (string, int, long, etc.)

* Results change interactively in the appropriate window

1. Open file "demog.csv". Add "Sex" column to field "Columns" and "Race" to "Rows" in "Aggregate rows" dialog. Select
   "avg(HEIGHT)" in "Measures" field. Use third-party software (e.g. MS Excel) to recount result

* Results of "Aggregate rows" module and in MS Excel are the same

1. Call context menu for "avg (HEIGHT)" in "Measures"

* You can change column and aggregation from context menu

1. Change column to "WEIGHT" and aggregation to "max". ("max (WEIGHT)")

* Result changed according to selected measure

1. Add one more column to field "Columns" (DIS_POP), field "Rows" (SEX) and field "Measures" (min (
   HEIGHT))

* Result changed according to added columns and measures

1. Use "Remove others" from context menu for one value from each field

* Values ​​besides one for which "Remove others" was called are deleted from dialog fields

1. Click on "History" button and return previous states

* "Aggregate rows" dialog return to previous state

1. Click on "Reset" button

* "Aggregate rows" return to start state. All fields became empty

1. Test non-functional modules (UI, popup menu, help, navigation, properties, etc.)

* Non-functional modules work correctly and are intuitive

See also:

* [Aggregate rows](aggregate-rows.md)
* [Aggregate rows Auto Test](aggregate-rows-test.side)
