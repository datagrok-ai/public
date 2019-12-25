<!-- TITLE: Tests: Workspace -->
<!-- SUBTITLE: -->

# Tests: Workspace

[Workspace](../features/workspace.md) is a collection of currently opened [projects](../entities/project.md).

[Workspace](../features/workspace) Tree is a hierarchical viewer that is used for navigating and managing [projects](../entities/project),
[tables](../entities/table.md), [connections](../entities/data-connection.md), and other entities.

## Testing scenarios

1. Open *"demog"* [project](../entities/project.md)

1. Open "[Workspace](../features/workspace.md)" via **View | Workspace** (or ```Alt + W```)
   * [Workspace](../features/workspace.md) is opened, in which empty [Skrachpad](../entities.md) and *"demog"* [project](../entities/project.md)

1. Expand *"demog"* in [Workspace](../features/workspace.md)
   * In *"demog"* branch there are subelements tables and view next to which number of internal elements is displayed (1)

1. Expand *"Tables"*
   * There is one table *"demog"*

1. Call context menu for table *"demog"*
   * "Add view" - add new view for table
   * "Remove" - remove table from [project](../entities/project.md)
   * "Open in Jupyter" - create "Jupyter Notebook" for table
   * "Save as CSV" - save table to .csv file
   * "Rename" - rename table
   * "Clone" - create copy of table
   * "Save as [project](../entities/project.md)"
   * "Exclude from [project](../entities/project.md)" - exclude table from [project](../entities/project.md). Table will be moved to [Scratchpad](../entities/scratchpad.md)

1. Expand table *"demog"*
   * Here are tabs *"Tags"* and *"Columns"*
   * In "Columns" you can see columns tags and columns stats
   * From context menu for columns, you can perform the operations "Rename", "Remove", "Extract", "Change type", "Extract numbers" (for sting type columns)

1. Expand "View" tab for [project](../entities/project.md) *"demog"*
   * Here are Table Views for all tables in [project](../entities/project.md)
   * You can remove view from its context menu

1. To add entities (queries, connections, data jobs, models, scripts, notebooks)to [project](../entities/project.md), 
   you can drag them from any location directly to [project](../entities/project.md) branch in [Workspace](../features/workspace.md)
   * After adding entity to [project](../entities/project.md), all available operations on them can be performed from 
     their context menu in [project](../entities/project.md) tree
   * Any entity can be excluded from [project](../entities/project.md). After that, it will be transferred to 
     [Scratchpad](../entities/scratchpad.md)

1. Open *"cars.csv"* file from local storage
   * File "cars.csv" is open
   * Table and view for it are added to [Scratchpad](../entities/scratchpad.md)

1. Drag *"cars"* table from [Scratchpad](../entities/scratchpad.md) to *"demog"* [project](../entities/project.md)
   * Table *"cars"* added to [project](../entities/project.md) *"demog"* and no longer exists in [Scratchpad](../entities/scratchpad.md)
   * Corresponding view was moved along with table in *"demog"* [project](../entities/project.md)

1. Choose *"Exclude"* from context menu of *"demog"* table in *"demog"* [project](../entities/project.md)
   * Table *"demog"* added to [Scratchpad](../entities/scratchpad.md) and no longer exists in *"demog"* [project](../entities/project.md)
   * Corresponding view was moved along with table in [Scratchpad](../entities/scratchpad.md)

1. Click on [Scratchpad](../entities/scratchpad.md) or on its internal elements. Then click on *"Upload"* button (or use context menu for upload)
   * "Create a project" dialog is open
   * New [project](../entities/project.md) open after uploading
   * All elements of [Scratchpad](../entities/scratchpad.md) have moved to new [project](../entities/project.md)

1. Click on *"demog"* [project](../entities/project.md) or on its internal elements. Then click on *"Save"* button (or use context menu for save)
   * [project](../entities/project.md) *"demog"* is reuploaded with changes

1. Click on *"Close others"* item in context menu for [Scratchpad](../entities/scratchpad.md)
   * All [projects](../entities/project.md) are closed


See also:
* [Workspace](../features/workspace.md)
* [Project](../entities/project.md)
* [Scratchpad](../entities/scratchpad.md)
* [Table](../entities/table.md)
