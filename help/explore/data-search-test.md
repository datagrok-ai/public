<!-- TITLE: Tests: Search -->
<!-- SUBTITLE: -->

# Tests: Search

Data [search](data-search.md) allows you to find data in a dataset. 
Here supported all patterns for numeric, string and datetime types. It is possible to select and filter the values found


## Testing scenarios


1. Open "demog" dataset
  
1. Open "[Search](data-search.md)" tab on toolbox 

1. Close "[Search](data-search.md)" tab on toolbox 

1. Press ```Ctrl + F``` from keyboard to open [Search](data-search.md).
Dataset should be in focus
   
1. Enter alternately *"50"*, *"Caucasian"*, *"1990"*, *"true"* in the [search](data-search.md) field
   * When entering data, the values ​​in the columns of the corresponding types are highlighted

1. Enter *"50"* in the [search](data-search.md) field. Press ```↓``` and ```↑``` keys on keyboard
   * After pressing ```↓``` key go to next *"50"* value 
   * After pressing ```↑``` key go to previous *"50"* value 
   
1. Enter *"50"* in the [search](data-search.md) field. Select "Select matching" item from 
   [search](data-search.md) menu (˅) or press ```Alt + S``` on keyboard
   * All rows in which the values found are selected

1. Enter *"50"* in the [search](data-search.md) field. Choose "Filter matching" item from search 
   menu (˅) or press ```Alt + F``` on keyboard
   * All rows in which the values found are selected are added to [filter](../viewers/filters.md)

1. Clear the [search](../features/data-search) field. Choose "Auto-select" item from [search](data-search.md) 
   menu (˅) and then enter *"50"* in the [search](data-search.md) field
   * After entering a value in the [search](data-search.md) field, rows with the found cells are selected automatically

1. Clear the [search](data-search.md) field. Choose "Auto-filter" item from 
   [search](data-search.md) menu (˅) and then enter *"50"* in the [search](data-search.md) field
   * After entering a value in the [search](data-search.md) field, rows with the found cells are added to 
     [filter](../viewers/filters.md) automatically

1. Use search [patterns](data-search-patterns.md)  for numerical values (*">"*, *"<"*. *"<="*, *">="*, *"range"*)
   * Numerical values ​​found correspond to the entered [pattern](data-search-patterns.md)
   
1. Use search [patterns](data-search-patterns.md)  for string values (*"starts with"*, *"ends with"*. *"contains"*, *"regex"*)
   * String values ​​found correspond to the entered [pattern](data-search-patterns.md)
   
1. Use search [patterns](../features/data-search-patterns)  for datetime values (eg, *"1990-1991"*, 
   *"June 1990"*. *"Oct 17, 1990"*, *"before 10/17/1990"*, *"after 10/17/1990"*, *"today"*, "this week", *"regex"*, 
   *"yesterday"* and etc. )
    * Datetime values ​​found correspond to the entered [pattern](data-search-patterns.md)

1.  Enter alternately *"Age = 50"*, *"Race contains as"*, *"Starred after 8/1/1990"*, *"Control true"* in the [search](data-search.md) field
    * Values ​​are searched only for the specified column in the [search](data-search.md) field

See also:
  * [Search](data-search.md)
  * [Search patterns](data-search-patterns.md)
  * [Smart search](../overview/smart-search.md)
  * [Search Auto Test](../selenium/data-search-test.side)
