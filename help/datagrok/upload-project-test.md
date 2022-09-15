<!-- TITLE: Tests: Upload project -->
<!-- SUBTITLE: -->

# Tests: Upload project

A [project](project.md) is a collection of different objects that you want to use together.
[Project](project.md) can contain [tables](table.md), [queries](../access/data-query.md),
[connections](../access/data-connection.md), [data jobs](../access/data-job.md),
[Jupyter notebooks](../compute/jupyter-notebook.md)
, [predictive models](../learn/predictive-modeling.md), and other.

## Testing scenario

1. Open a local file. (for example, test_tables.xlsx)

* File is opened
* Corresponding views are created.
* Tables and their views added to the "Scratchpad" project. (for check: **View | Workspace**)

1. Run query "Test" from "northwind" connection of MySQL source

* New table (result of query) and its view added to the "Scratchpad" project.

1. Add all available viewers to views of opened [tables](table.md).

1. Click on "Upload" button

* "Create project" dialog is open

1. Enter the name and description for new [project](project.md). (for example, name="test"). Add connection, query and
   data job to project.

* Name and description have been changed successfully. [query](../access/data-query.md),
  [connections](../access/data-connection.md) and [data jobs](../access/data-job.md) added.

1. Click on "Ok" button of "Publish project" dialog

* "Share to" dialog is open

1. Click on "Copy to Clipboard" icon

* Project link copied
* Balun "Link copied to clipboard" is shown

1. Select "All users" for allow permissions to view the project

* Field "Send notification" appeared

1. Mark "Send notification" as true

* Label with number of users who receive message appeared

1. Share project with all users

* All users now have access to view project and received emails

1. Save project locally as ZIP-file

1. Open link copied in step 7 in browser

* Uploaded project open in presentation mode

1. Return to datasets browser

1. Delete "test" project

* "test" project deleted successfully

1. Open downloaded ZIP-file from step 11

* Project from ZIP-file open successfully

1. Upload opened project to server

See also:

* [Prject](project.md)
* [Projects tutorial](../_internal/tutorials/projects.md)
