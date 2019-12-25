<!-- TITLE: Tests: Hide Menus -->
<!-- SUBTITLE: -->

# Test: Hide Menus

Platform flexibility allows to hide unused items of various menus (top menu, context menu for objects, etc.)

## Testing scenario

1. Call context menu for "Chem" item of top menu
   * Called context menu to hide selected item
   * Context menu displays item name to be hidden

1. Click on "Hide" from from called context menu for "Chem"
   * "Chem" from top menu is hidden
   * Balun about item hiding is shown

1.  Call context menu for "Admin | Users" (or any other menu item) item of top menu
   * Called context menu to hide selected item
   * Context menu displays item name to be hidden

1. Click on "Hide" from from called context menu for "Admin | Users"
   * "Admin | Users" from top menu is hidden
   * Balun about item hiding is shown
   
1. Call context menu for "demog" (or any other project)

1. Right-click on item "Share..."
   * Called menu to hide selected item
   * Displays item name to be hidden
   
1. Click on "Hide" from from called menu for "Share..."
   * "Share..." from context menu for projects is hidden
   * Balun about item hiding is shown
   
1. Call context menu for other various projects
   * "Share..." is not displayed in context menu for all projects
   
1. As described in previous steps, hide one context menu item for all platform objects (Data source, Data connection, Data Query, Data Job, User, Group, Script, Notebook, Model, Function)

1. Refresh the platform page
   * After refreshing all previously hidden menu items are not displayed

1. Open **Tools | Settings | Menu**
   * All hidden menu items are displayed in format: ```top | Chem```  

1. Remove item "top | Chem" from **Tools | Settings | Menu** and left other items

1. Refresh the platform page
    * "Chem" item again displayed in top menu, while other hidden menus are not displayed
    
1. Clear all from  **Tools | Settings | Menu**   

1. Refresh the platform page
    * All previously hidden items now displayed