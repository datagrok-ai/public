<!-- TITLE: Tests: Collections -->
<!-- SUBTITLE: -->

# Tests: Collections

Collections are sets of entities for which you can perform same actions.

## Testing scenario

1. Open "Datasets..." page of "Welcome" view

1. Click on "Demog" project
   * "Demog" project selected on [Property Panel](../overview/property-panel.md)
     
1. Click with ```CTRL``` on "Demo" project
   * Projects "Demog" and "Demo" added to collection
   * [Property Panel](../overview/property-panel.md) switched to "2 projects"
   * Actions and commands on [Property Panel](../overview/property-panel.md) concern both projects from collection
   
1. Expand "Items" tab
   * All projects from collection are presented in "Items" tab
   
1. Expand "Commands" tab
   * Tab contains available actions for all collection projects
   
1. Click on "Add to favorites" action
   * Both projects added to "Favorites"

1. Click on "Tag" action
   * Tag adding dialog is open
      
1. Write \#test_collection tag and click OK
   * \#test_collection added to both projects
   
1. Click on "Share" action
   * "Share" dialog open for each project
   
1. Give access to view for "All users" for each project
   * "All users" received rights to view both projects

1. Click on "Open" action
   * Both projects have opened
   
1. Click on "Details" action
   * "Details" tab is open for each project
   
1. Click on "Delete" action
   * "Are you sure?" dialogs opened for each project
   
1. Cancel project deletion    

1. Add one more project to collection by clicking with ```CTRL```

1. Repeat 4-14 steps for collection with three projects

1. Make sure 3 projects are open (if not, open 3 projects)

1. Select into collection two of three opened projects

1. Click on "Close Others" action
   * Project that was not selected in collection is closed
   
1. Open **Admin | Data connections**

1. Select "northwind" (PostgreSQL) and "AirNow" (Web) for collection by clicking with ```CTRL```
   * Connections "northwind" (PostgreSQL) and "AirNow" (Web) added to collection
   * [Property Panel](../overview/property-panel.md) switched to "2 connections"
   * Actions and commands on [Property Panel](../overview/property-panel.md) concern both connections from collection
   
1. Expand "Items" tab
   * All connections from collection are presented in "Items" tab
   
1. Expand "Commands" tab
   * Tab contains available actions for all collection projects
   
1. Run all actions in turn from "Commands" tab for connections collection 
   * All actions runs to all connections in collection
   
1. Change browser view to "Card" and add connections to collection 
   * Connections are added to collection in "Card" view

1. Change browser view to "Table" and add connections to collection 
   * Connections are added to collection in "Table" view

1. Repeat 21-26 steps for all entities browsers (Scripts, Models, Layouts, Queries, Query runs, Jobs, Job runs, Users, Groups)
