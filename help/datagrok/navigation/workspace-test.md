<!-- TITLE: Tests: Workspace -->
<!-- SUBTITLE: -->

# Tests: workspace

Workspace is a collection of currently opened [projects](../concepts/project/project.md).

[Browse](../../datagrok/navigation/views/browse.md) Tree is a hierarchical viewer that is used for navigating and
managing [projects](../concepts/project/project.md),
[tables](../concepts/table.md), [connections](../access/access.md#data-connection), and other entities.

## Testing scenarios

1. Open *"demog"* [project](../concepts/project/project.md)

1. Open "Workspace" via **View | Workspace** (or ```Alt + W```)

* Workspace is opened, in which empty Scratchpad and *"
  demog"* [project](../concepts/project/project.md)

1. Expand *"demog"* in Workspace

* In *"demog"* branch there are subelements tables and view next to which number of internal elements is displayed (1)

1. Expand *"Tables"*

* There is one table *"demog"*

1. Call context menu for table *"demog"*

* "Add view" - add new view for table
* "Remove" - remove table from [project](../concepts/project/project.md)
* "Open in Jupyter" - create "Jupyter Notebook" for table
* "Save as CSV" - save table to .csv file
* "Rename" - rename table
* "Clone" - create copy of table
* "Save as [project](../concepts/project/project.md)"
* "Exclude from [project](../concepts/project/project.md)" - exclude table from [project](../concepts/project/project.md). Table will be moved
  to Scratchpad

1. Expand table *"demog"*

* Here are tabs *"Tags"* and *"Columns"*
* In "Columns" you can see columns tags and columns stats
* From context menu for columns, you can perform the operations "Rename", "Remove", "Extract", "
  Change type", "Extract numbers" (for sting type columns)

1. Expand "View" tab for [project](../concepts/project/project.md) *"demog"*

* Here are Table Views for all tables in [project](../concepts/project/project.md)
* You can remove view from its context menu

1. To add entities (queries, connections, data jobs, models, scripts, notebooks)
   to [project](../concepts/project/project.md), you can drag them from any location directly to [project](../concepts/project/project.md)
   branch in Workspace

* After adding entity to [project](../concepts/project/project.md), all available operations on them can be performed from their context
  menu in [project](../concepts/project/project.md) tree
* Any entity can be excluded from [project](../concepts/project/project.md). After that, it will be transferred to
  Scratchpad

1. Open *"cars.csv"* file from local storage

* File "cars.csv" is open
* Table and view for it are added to Scratchpad

1. Drag *"cars"* table from Scratchpad to *"demog"* [project](../concepts/project/project.md)

* Table *"cars"* added to [project](../concepts/project/project.md) *"demog"* and no longer exists in Scratchpad
* Corresponding view was moved along with table in *"demog"* [project](../concepts/project/project.md)

1. Choose *"Exclude"* from context menu of *"demog"* table in *"demog"* [project](../concepts/project/project.md)

* Table *"demog"* added to Scratchpad and no longer exists in *"
  demog"* [project](../concepts/project/project.md)
* Corresponding view was moved along with table in Scratchpad

1. Click on Scratchpad or on its internal elements. Then click on *"Upload"*
   button (or use context menu for upload)

* "Create a project" dialog is open
* New [project](../concepts/project/project.md) open after uploading
* All elements of Scratchpad have moved to new [project](../concepts/project/project.md)

1. Click on *"demog"* [project](../concepts/project/project.md) or on its internal elements. Then click on *"Save"*
   button (or use context menu for save)

* [project](../concepts/project/project.md) *"demog"* is reuploaded with changes

1. Click on *"Close others"* item in context menu for Scratchpad

* All [projects](../concepts/project/project.md) are closed

See also:

* [Browse](../../datagrok/navigation/views/browse.md)
* [Project](../concepts/project/project.md)
* [Table](../concepts/table.md)
