<!-- TITLE: Tests: Sorting -->
<!-- SUBTITLE: -->

# Tests: Sorting

Platform provides ability to sort data by values ​​of columns. 
Supports sorting in descending, ascending, and default order (original order from data source).

## Testing scenario

1. Open *demog* table

1. Remove multiple values ​​from *RACE*, *STARTED* and *CONTROL* columns
   * Columns of each type contains ```null``` values

1. Double click on *AGE* column
   * Order of table rows has changed by decreasing values ​​of *AGE* column
   * In heading of *AGE* column shows down arrow
   * Rows with ```AGE ​​= null``` moving to bottom of table   

1. Double click on *AGE* column
   * Order of table rows has changed by ascending values ​​of *AGE* column
   * In heading of *AGE* column shows up arrow
   * Rows with ```AGE ​​= null``` moving to bottom of table   
      
1. Double click on *AGE* column
   * Order of table rows has changed to default order (original order from data source).

1. Double click on *RACE* column
   * Order of table rows has changed by decreasing values ​​of *RACE* column (by reverse alphabet)
   * In heading of *RACE* column shows down arrow
   * Rows with ```RACE ​​= null``` moving to bottom of table   

1. Double click on *RACE* column
   * Order of table rows has changed by ascending values ​​of *RACE* column (by alphabet)
   * In heading of *RACE* column shows up arrow
   * Rows with ```RACE ​​= null``` moving to bottom of table   
      
1. Double click on *RACE* column
   * Order of table rows has changed to default order (original order from data source).   
   
1. Double click on *HEIGHT* column
   * Order of table rows has changed by decreasing values ​​of *HEIGHT* column 
   * In heading of *HEIGHT* column shows down arrow
   * Rows with ```HEIGHT ​​= null``` moving to bottom of table   

1. Double click on *HEIGHT* column
   * Order of table rows has changed by ascending values ​​of *HEIGHT* column
   * In heading of *HEIGHT* column shows up arrow
   * Rows with ```HEIGHT ​​= null``` moving to bottom of table   
      
1. Double click on *HEIGHT* column
   * Order of table rows has changed to default order (original order from data source).   
   
1. Double click on *CONTROL* column
   * Order of table rows has changed by decreasing values ​​of *CONTROL* column (first values ​​```true```)
   * In heading of *CONTROL* column shows down arrow
   * Rows with ```CONTROL ​​= null``` moving to bottom of table   

1. Double click on *CONTROL* column
   * Order of table rows has changed by ascending values ​​of *CONTROL* column (first values ​​```false```)
   * In heading of *CONTROL* column shows up arrow
   * Rows with ```CONTROL ​​= null``` moving to bottom of table   
      
1. Double click on *CONTROL* column
   * Order of table rows has changed to default order (original order from data source).         
   
1. Double click on *STARTED* column
   * Order of table rows has changed by decreasing values ​​of *CONTROL* column (by reverse chronology)
   * In heading of *CONTROL* column shows down arrow
   * Rows with ```CONTROL ​​= null``` moving to bottom of table   

1. Double click on *STARTED* column
   * Order of table rows has changed by ascending values ​​of *STARTED* column (by chronology)
   * In heading of *STARTED* column shows up arrow
   * Rows with ```STARTED ​​= null``` moving to bottom of table   
      
1. Double click on *STARTED* column
   * Order of table rows has changed to default order (original order from data source).         
      
Also use context menu of columns to sorting. Submenu **Sort**. 
Repeat steps 3-17 using the column context menu.   
   
See also:
  * [Table](../entities/table.md)
  * [Viewers Test](../../visualize/viewers/viewers-test.md)
