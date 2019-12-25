<!-- TITLE: Tests: Upload project -->
<!-- SUBTITLE: -->

# Tests: Upload project

A [project](../entities/project.md) is a collection of different objects that you want to use together. 
[Project](../entities/project.md) can contain [tables](../entities/table.md), [queries](../entities/data-query.md), 
[connections](../entities/data-connection.md), [data jobs](../entities/data-job.md), 
[Jupyter notebooks](../plugins/jupyter-notebook.md), [predictive models](../plugins/predictive-modeling.md), and other.

## Testing scenario


1. Open a local file. (for example, test_tables.xlsx) 
   * File is opened 
   * Corresponding views are created. 
   * Tables and their views added to the "Scratchpad" project. (for check: **View | Workspace**)

1. Run query "Test" from "northwind" connection of MySQL source
   * New table (result of query) and its view added to the "Scratchpad" project. 

1. Add all available viewers to views of opened [tables](../entities/table.md). 

1. Click on "Upload" button
   * "Create project" dialog is open

1. Enter the name and description for new [project](../entities/project.md). (for example, name="test"). 
   Add connection, query and data job to project. 
   * Name and description have been changed successfully. [query](../entities/data-query.md), 
     [connections](../entities/data-connection.md) and [data jobs](../entities/data-job.md) added.

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
   * All users now have access to view project and received e-mails

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
 * [Prject](../entities/project.md)
 * [Projects tutorial](../tutorials/projects.md)
 * [Project browser test](../tests/project-browser-test.md)
 