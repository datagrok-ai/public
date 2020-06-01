# Docking

Datagrok provides a flexible window management system, where windows
could be either dragged out and positioned manually, or set up automatically.

## Top-level window docking

Use `grok.shell.dockManager` to dock, undock, or reposition windows. 
See also `View.dockNode`   

## Nested docking

Some of the views contain a nested docking manager, which allows to manage
windows within that particular view (they cannot be undocked and docked on
a top level, or within a different view). See `dockManager` property of the 
View class.
