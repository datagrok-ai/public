<!-- TITLE: Navigation -->
<!-- SUBTITLE: -->

# Navigation

## Views, toolbox, toolbar, property panel

Typically, a view resides in the center and occupies all available area in the screen. A toolbox 
located on the left contains controls that are specific to that view.

Toolbar, located on the top right under the menu, has some static controls (such as 'import file' icon),
as well as controls that are specific for the active view. 

All panels can be turned on or off on/off by either clicking on the 'x' sign on the top right, 
using shortcuts, or via the menu:View | Toolbox    - Alt+XView | Properties - F4Help | Context Help - F1

Shortcuts

|                     |                             |                        |
|---------------------|-----------------------------|------------------------|
| View / Toolbox      | Alt+X                       | Show / hide toolbox    |
| View / Properties   | F4                          | Show / hide property panel  | 
| Help / Context Help | F1                          | Show / hide help   | 
|                     | Double-click on view header | Open view in full-screen | 


## Current object

Whenever user clicks on most objects in the platform, such as a viewer, view, command, or any
of the [40+ types of objects](../entities/entities.md) supported by the platform, this object becomes
what we call "current object". Its properties are shown in the [property panel](../features/property-panel.md),
and interactive help is displayed as well.  

To access current object from the console, use 'o' variable.

## Drag-and-Drop

Datagrok lets almost all objects be drag-and-dropped, which makes a lot of interactions more intuitive.
As soon as a user starts dragging an object, all possible drop zones get highlighed. Here
are some of commonly used drag-and drop operations: 

* [Variables view](../overview/variables-view.md)
* Columns 
  * Rearrange grid column order (drag by the grid column header)
  * Drag into [Add New Column](../transform/add-new-column.md) dialog

## Drag-and-Dropping Columns

Drag-and-dropping columns is an intuitive way to choose columns on viewers, or perform operations on
them. A number of controls can be used as a source of drag-and-drop: 
[Grid](../viewers/grid.md), [Column Manager](../explore/column-manager.md), [Table Manager](../overview/table-manager.md)  
Most of the viewers, and a lot of dialogs accept columns that are being dragged. 

## Help

Typically, when you click on an object it becomes a "current object", and its help appears in the help pane. 
Press `F1` to toggle its visibility. Use `<`, `>` icons on top of the help pane to move back/forward. 
Press `◰` to open that page in a new tab.

See also:

* [Property panel](property-panel.md)
* [Toolbox](../overview/toolbox.md)
* [Info panels](../discover/info-panels.md)
