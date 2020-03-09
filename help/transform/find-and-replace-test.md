<!-- TITLE: Tests: Find and Replace -->
<!-- SUBTITLE: -->

# Tests: Find and Replace

This typical ["Find and Replace"](find-and-replace-test.md) dialog that you see in every text editor, except that
it supports specifying the columns to run against, and \[search patterns\] for matching non-textual columns.

## Testing scenario

1. Open *"demog"* dataset

1. Select 10 first rows in *"demog"* table

1. Filter rows in *"demog"* table for ```Age = [35, 45]```

1. Select column "SITEID" in current focus

1. Open [Find and Replace](find-and-replace-test.md) dialog from **Edit** menu. (Or use ```Ctrl + H```)
   * [Find and Replace](find-and-replace-test.md) dialog is open
   * Help switched to [Find and Replace](find-and-replace-test.md)

1. Enter value ```1``` to "Find what" field

1. Enter value ```1000``` to "Replace with" field

1. Select ```Current``` value in "In columns" field

1. Set values of "Selected row only" and "Filtered rows only" as ```true``` (remaining bool fields should be ```false```)

1. Click on ```Ok``` button for execute dialog
   * "SITEID" column values which were equal to 1 replaced by 1000
   * Only filtered values ​​are changed for which ```Age = [35, 45]```
   * Only values ​​for first ten rows are replaced which were selected
   
1. Remove all rows filters and  selections. 
   * Make sure that values ​​that were not selected and not filtered are not changed. ("SITEID" column values remained 1)   

1. Open [Find and Replace](find-and-replace-test.md) dialog 
   
1. Select column "RACE" in current focus  

1. Enter value ```asian``` to "Find what" field

1. Enter value ```test``` to "Replace with" field

1. Select ```Current``` value in "In columns" field

1. Set value of "Match case" as ```true``` (remaining bool fields should be ```false```)

1. Click on ```Ok``` button for execute dialog
   * Values ​​in "RACE" column have not replaced

1. Repeat 12-18 steps. In step 14, enter ```Asian``` to "Find what" field 
    * Values ```Asian``` ​​in "RACE" column have replaced by values ```test```

1. Open [Find and Replace](find-and-replace-test.md) dialog 

1. Select column "DIS_POP" in current focus  

1. Enter value ```A``` to "Find what" field

1. Enter value ```test``` to "Replace with" field

1. Select ```Current``` value in "In columns" field

1. Set value of "Match whole word" as ```true``` (remaining bool fields should be ```false```)

1. Click on ```Ok``` button for execute dialog
   * Values ​​in "DIS_POP" column have not replaced

1. Repeat 20-26 steps. In step 22, enter ```RA``` to "Find what" field 
    * Values ```RA``` ​​in "DIS_POP" column have replaced by values ```test```

1. Open [Find and Replace](find-and-replace-test.md) dialog 

1. Select column "SEX" in current focus

1. Set value of "Filter matching rows" as ```true``` (remaining bool fields should be ```false```)   

1. Enter value ```F``` to "Find what" field

1. Enter value ```test``` to "Replace with" field

1. Click on ```Ok``` button for execute dialog
   * Values ​​in ```F```  "SEX" column have replaced by values ```test```
   * Filter has been applied to rows with replaced values.

1. Open [Find and Replace](find-and-replace-test.md) dialog 
   
1. Use search [patterns](../explore/data-search-patterns.md)  for numerical values (*">"*, *"<"*. *"<="*, *">="*, *"range"*) in "Find what" field
   * Numerical values ​​found correspond to the entered [pattern](../explore/data-search-patterns.md)
   
1. Use search [patterns](../explore/data-search-patterns.md)  for string values (*"starts with"*, *"ends with"*. *"contains"*, *"regex"*)
   * String values ​​found correspond to the entered [pattern](../explore/data-search-patterns.md)
   
1. Use search [patterns](../explore/data-search-patterns.md)  for datetime values (eg, *"1990-1991"*, *"June 1990"*. *"Oct 17, 1990"*, *"before 10/17/1990"*, *"after 10/17/1990"*, *"today"*, "this week", *"regex"*, *"yesterday"* and etc. )
    * Datetime values ​​found correspond to the entered [pattern](../explore/data-search-patterns.md)

1. Test non-functional modules (help, navigation, UI, UX)
   * Non-functional modules work correctly and are intuitive


See also:
 * [Find and Replace](find-and-replace.md)
