<!-- TITLE: Project -->
<!-- SUBTITLE: -->

# Project

Project is a collection of entities along with the applied visualizations.
Projects are used to group and share data and other assets with other users. One of the most
common applications of projects are dashboards that consist of tables (with either
static or dynamic data), and visualizations applied to them.

## Uploading a project

Creating a project is easy. After getting the data of interest in the scratchpad project in [workspace](workspace.md),
click on the `UPLOAD` button. After the project gets uploaded to the server, a separate 
window pops us asking you whom to share the project with. By default, it is only accessible by you,
you have to share it in order for others to use it.

Or, if you are editing an existing project, click `SAVE` to save your changes.

Use `Share` context action to edit access permissions. Sharing a project will 
automatically share all entities and data inside.

### Dynamic data

Whenever a table is created by executing a [function](../overview/functions) 
(such as a [database query](../access/data-query.md)), this information gets stored 
with the table as a "generation script". This serves multiple purposes:
* Provides lineage regarding
* On-demand data refreshing (Table toolbox, "Query" panel, `REFRESH` button)
* Enables publishing dashboards with the dynamic data
  
In the "Upload project" dialog, a "Data sync" option appears next to the tables
that have a generation script defined. This option determines whether the data  
should be stored as a static snapshot, or as a generation script.
In the latter case, the function will be re-executed whenever the project is opened.     

![](project-upload-data-sync.png)

## Filtering

The following fields could be used to filter projects with [smart search](smart-search.md):

| Field        | Description                                 |
|--------------|---------------------------------------------|
| name         |                                             |
| description  |                                             |
| createdOn    |                                             |
| updatedOn    |                                             |
| author       | [User](../govern/user.md) object            |
| starredBy    | [User](../govern/user.md) object            |
| commentedBy  | [User](../govern/user.md) object            |
| usedBy       | [User](../govern/user.md) object            |

See also:

  * [Create Project](create-project.md)
  * [Project Gallery](project-gallery.md)
  * [Data Pipeline](../access/data-pipeline.md)
  * [Data Connection](../access/data-connection.md)
  * [Data Query](../access/data-query.md)
  * [Data Job](../access/data-job.md)
  * [Function Call](functions/function-call.md)
