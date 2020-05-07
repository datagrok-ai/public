<!-- TITLE: Tests: Column Manager -->
<!-- SUBTITLE: -->

# Tests: Column Manager

[Column Manager](column-manager.md), available via **View | Columns**, contains a list of columns in the currently 
open tables. Use it to navigate between columns, select them, or perform batch actions. 
It also allows to view metadata on multiple columns in a tabular format, as well as other column
properties, such as statistics.

The implementation is based on the [grid](../visualize/viewers/grid.md), so many of the grid's features apply.

## Testing scenario

1. Open demog.csv

1. Open "Columns manager" via **View | Tables** (or ```Alt + C```)

1. Alternately click on columns in ["Column Manager"](column-manager.md)
   * Selected column is highlighted in [table](../overview/table.md)
   * Selected column is shown on [Property Panel](../overview/property-panel.md)
   
1. Call context menu for column in ["Column Manager"](column-manager.md)
   * Submenu is displayed that refers to selected column (tittle like *'column name' + "column"*, eg. *'Age' column*), as well as "Add stats" submenu for all columns
   * There are items "Remove", "Rename" "Extract", "Change type", "Rename" in submenu for selected column, which work correctly
        
1. With ```Shift``` held down, select all columns from ["Column Manager"](column-manager.md)
   * In context menu, submenu appears whose actions apply to all selected columns (tittle like *"n selected columns"*)
   * [Property Panel](../overview/property-panel.md) displays actions related to all selected columns
   
1. Expand "Add stats" submenu from the context menu
   * Submenu "Add stats" shows tree blocks of columns attributes to be added to ["Column Manager"](column-manager.md)
   * First block contains items "All" (for add all attributes), "Selected" and "Filtered", which work correctly 
   * Second block with "Histogram" item for add the histogram for numerical columns (for non-numeric columns cells are empty)
   * Third block contains standart columns stats (min, max, avg, etc.)
   
1. Click on "All" item from "Add stats" submenu from the context menu
   * All colimns attributes are added to ["Column Manager"](column-manager.md)
   * "All" item and and all attributes in list marked with icon "✓"
   * If attribute is not available for a specific column, the corresponding cell is empty
 
1. Click on "All" item from "Add stats" submenu from the context menu again
   * Attributes are no longer displayed in ["Column Manager"](column-manager.md)
   * "All" item and all attributes are no longer marked with icon "✓"
   
1. Alternately add individual column attributes from all blocks
   * When you add column attribute, corresponding column is added to ["Column Manager"](column-manager.md)
   * Added column attribute is marked with icon "✓" in the list
   * When you click on attribute again - corresponding column is deleted from ["Column Manager"](column-manager.md) and icon "✓" for it no longer displayed

1. Open ["Aggregate Rows"](../transform/aggregate-rows.md) dialog. Drag columns from ["Column Manager"](column-manager.md) 
   to ["Aggregate Rows"](../transform/aggregate-rows.md) dialog fields
   * ["Aggregate Rows"](../transform/aggregate-rows.md) dialog fields are filled with dragged columns
   * In the same way, test the ability to drag and drop into all dialogs and [viewers](../visualize/viewers.md) that it supports

1. Test non-functional features (help, navigation, working with window, etc.)
   * Non-functional modules work correctly and are intuitive   


See also: 
  * ["Column Manager"](column-manager.md)
  * [Table Manager Test](../overview/table-manager-test.md)
