---
title: "Tests"
---

## Testing viewers

Each viewer should at least satisfy a number of well-known criteria. Some can (and should) be tested
automatically, others should be tested manually. Each viewer:

* Should open against any dataset. If the data structure is not applicable, a message is shown inside the viewer.
* Changing settings should not produce exceptions
* Changing settings should be immediately reflected on a viewer
* Should save and load from a layout
* Should save and load from a project
* Should not freeze browser. If you try to visualize too many rows, limit the number 
  (and show a message inside), or show a message (inside the viewer) that it can't be visualized.
... and a lot more.

## Testing functions

Functions are no different

* Function browser should not produce errors when building the UIs for any function
* Most important functions should have test cases specified as a function signature
* Depending on a function role, a function should behave appropriately:
  * Function role: cell renderer
    * Should have the right signature
    * Should execute and return an instance of GridCellRenderer