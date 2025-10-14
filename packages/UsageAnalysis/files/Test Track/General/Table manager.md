### Table Manager

**Table Manager**, available via **View | Tables**, contains a list of currently open tables. Use it to navigate between tables, select them, or perform batch actions. It also allows to view metadata on multiple tables in
a tabular format. 

The implementation is based on the **grid**, so many of the grid's features apply.

- Precondition: load 3-4 radom datasets.

1. Open "Table" dialog via **View | Tables** (or ```Alt + T```): 
* A list of opened datasets opens
* When pressing on a dataset in the Table Manager, tab with the selected table opens.
* Call context panel for dataet in Table Manager. Selected dataset in Table Manager is displayed on Context Panel. 

2. Click on "Save as table" item in context menu:
* New tab has been created with table constaining information from Table Manager`s datasets.
* New row has been created in Table Manager.

3. With ```Shift``` held down, select all added tables from Table Manager:
* In context menu, submenu appears whose actions apply to all selected tables (tittle is "3 tables")
* Context Panel displays actions related to all selected tables

4. Expand "Show" submenu pressing RMB on the Table Manager view. Submenu "Show" shows tables attributes to be added to Table Manager. There is item for adding all attributes ("All").
- Click on "All" item from "Show" submenu from the context menu
  - All tables attributes are added to Table Manager
  - "All" item and and all attributes in list marked with icon "✓". (If attribute is not available for a specific table, the corresponding cell is empty)
- Click on "All" item from "Show" submenu from the context menu again:
  - Tables attributes are no longer displayed in Table Manager
  - "All" item and all attributes are no longer marked with icon "✓"

---
{
  "order": 6
}