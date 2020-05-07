<!-- TITLE: Tests: Table Manager -->
<!-- SUBTITLE: -->

# Tests: Table Manager

["Table Manager"](table-manager.md), available via **View | Tables**, contains a list of currently open tables. Use it to 
navigate between tables, select them, or perform batch actions. It also allows to view metadata on
multiple tables in a tabular format.

The implementation is based on the [grid](../visualize/viewers/grid.md), so many of the grid's features apply.

## Testing scenarios

1. Open two [tables](table.md)

1. Open "Table" dialog via **View | Tables** (or ```Alt + T```)

1. Click on first [table](table.md) in ["Table Manager"](table-manager.md) then click on second table
   * When selecting a [table](table.md) from the ["Table Manager"](table-manager.md), tab with the selected table opens
   * Selected [table](table.md) in ["Table Manager"](table-manager.md) is displayed on [Property Panel](../overview/property-panel.md)
   
1. Call context menu for [table](table.md) in ["Table Manager"](table-manager.md)
   * Submenu is displayed that refers to selected [table](table.md) (tittle like table name), as well 
     as common submenus for all tables ("Show", "Save as table")
   * There are items "Add view", "Remove" "Open in Jupyter", "Save as CSV", "Rename", "Clone", "Save as project" in 
     submenu for selected [table](table.md), which work correctly
        
1. Click on "Save as table" item in context menu
   * New tab has been created with [table](table.md) in which information from 
     ["Table Manager"](table-manager.md)
   * New row has been created in ["Table Manager"](table-manager.md) corresponding to new table
   
1. With ```Shift``` held down, select all three tables from ["Table Manager"](table-manager.md)
   * In context menu, submenu appears whose actions apply to all selected tables (tittle is "3 tables")
   * [Property Panel](../overview/property-panel.md) displays actions related to all selected tables
   
1. Expand "Show" submenu from the context menu
   * Submenu "Show" shows tables attributes to be added to ["Table Manager"](table-manager.md)
   * There is item for adding all attributes ("All")
   
1. Click on "All" item from "Show" submenu from the context menu
   * All tables attributes are added to ["Table Manager"](table-manager.md)
   * "All" item and and all attributes in list marked with icon "✓"
   * If attribute is not available for a specific [table](table.md), the corresponding cell is empty
   
1. Click on "All" item from "Show" submenu from the context menu again
   * Tables attributes are no longer displayed in ["Table Manager"](table-manager.md)
   * "All" item and all attributes are no longer marked with icon "✓"
   
1. Alternately add individual [table](table.md) attributes 
   * When you add [table](table.md) attribute, corresponding column is added to ["Table Manager"](table-manager.md)
   * Added [table](table.md) attribute is marked with icon "✓" in the list
   * When you click on [table](table.md) attribute again - corresponding column is deleted from 
     ["Table Manager"](table-manager.md) and icon "✓" for it no longer displayed

1. Test non-functional features (help, navigation, working with window, etc.)
   * Non-functional modules work correctly and are intuitive   


See also: 
  * [Table Manager](table-manager.md)
  * [Column Manager Test](../explore/column-manager-test.md)
