<!-- TITLE: Tests: Project Browser -->
<!-- SUBTITLE: -->

# Tests: Project Browser

Project browser allows you to view previously created [projects](project.md), 
view their properties and manage them.

6.2 Browse and Open the project

## Testing scenarios

1. Open a "Datasets" view for "Welcome" tab
   * List of saved projects appeared on the "Welcome" tab

1. Use [search](smart-search.md) to find the "Test" [project](project.md)

1. Open "Share" dialog for the "Test" [project](project.md) from its context menu
   * "Share" dialog open
   * You can change users or groups permissions for project here.

1. Open "Details" dialog for the "Test" [project](project.md) from its context menu
   * "Details" dialog open
   * Here is the description of the project, [project](project.md) elements and activity.  

1. Open the "General" tab in [Property Panel](../features/property-panel.md) for [project](project.md)
   * "General" tab is open
   * The correct and actually information for all fields is displayed (Created, Created by)

1. Open the "Tables" tab in [Property Panel](../features/property-panel.md) for [project](project.md)
   * "Tables" tab is open
   * Display all [tables](table.md) which are included in [project](project.md)

1. Open the "Connections" tab in [Property Panel](../features/property-panel.md) for [project](project.md)
   * "Connections" tab is open
   * Display all [connections](../access/data-connection.md) which are included in [project](project.md)

1. Open the "Queries" tab in [Property Panel](../features/property-panel.md) for [project](project.md)
   * "Queries" tab is open
   * Display all [queries](../access/data-query.md) which are included in [project](project.md)

1. Open the "Jobs" tab in [Property Panel](../features/property-panel.md) for [project](project.md)
   * "Jobs" tab is open
   * Display all [data jobs](../access/data-job.md) which are included in [project](project.md)

1. Open the "Activity" tab in [Property Panel](../features/property-panel.md) for [project](project.md)
   * "Activity" tab is open
   * Display information about actual actions with [project](project.md)

1. Open the "Shared with" tab in [Property Panel](../features/property-panel.md)
   * "Shared with" tab is open
   * Display [users](../govern/user.md) and users [groups](../govern/group.md) which this [project](project.md) is available for view or edit

1. Open the "Resources" tab in [Property Panel](../features/property-panel.md)
  
1. Open the "Notebooks" tab in [Property Panel](../features/property-panel.md) 
   * "Notebooks" tab is open
   * Display all [notebooks](../compute/jupyter-notebook.md) which are included in [project](project.md)

1. Open the "Models" tab in [Property Panel](../features/property-panel.md) 
   * "Models" tab is open
   * Display all [models](../learn/predictive-modeling.md) which are included in project

1. Open [project](project.md) "Test" (by double clicking or from context menu)
   * [Project](project.md) "Test" is open
   * All [project](project.md) views were opened, all created [viewers](../viewers/viewers.md) on them opened correctly and in full

1. Open **View | Workspace**  ( ``` Alt + W ```)
   * In it there is a tree of the [project](project.md) "Test" and the default project "[Scratchpad](scratchpad.md)"

1. Test all elements of the context menu of the "Test" [project](project.md) in the  "Workspace" window

1. Expand tree of "Test" [project](project.md) in "Workspace" window
   * Here are all the elements created in the project ([tables](table.md), [views](../visualize/view-layout.md), 
     [connections](../access/data-connection.md), [queries](../access/data-query.md), etc.)

1. Test the work of project elements. (open views, run queries, run jobs, etc.)

See also:
 * [Project upload test](../dialogs/upload-project-test.md)
 * [Project](project.md)
 * [Projects tutorial](../_internal/tutorials/projects.md)
