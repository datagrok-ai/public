@realizes:views.projects
Feature: Projects regressions: a pivot saved after its source table is closed
  core#531: orders.csv from Browse > Files > Demo > northwind gets a pivot table (Group by
  ShipCountry, Pivot ShipVia, count(OrderID)) published into the workspace with ADD, and the orders
  view is then closed. The Save dialog still offers Data sync for the published table, whose creation
  script aggregates the closed "orders". A project saved with it on cannot be opened again: the
  reopen stops on "Data loading error" with "Could not resolve table "orders"". The scenario is a known
  failure until the dialog disables Data sync for a table whose source is gone.

  The pivot is configured through its own rows and chips, as in projects-uploading. The project is
  named with the run's time and removed with its table and view when the scenario starts and when the
  feature ends.

  Background:
    Given user is logged in
    And simple mode is off
    And the browse panel is open

  # core#531: Data sync is offered for "orders aggregation" after "orders" is closed
  @known-failure:49 @realizes:core#531
  Scenario: Data sync is not offered for a pivot whose source table is closed
    Given no project named "BDDRegClosedSrc{time}" is on the server
    And user clears the saved pivot table parameters
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    And Files---Demo---northwind tree node inside browse tree is expanded
    When user double-clicks Files---Demo---northwind---orders.csv tree node inside browse tree
    Then the "orders" view should be current
    And table "orders" should have 830 rows
    Given the toolbox pane is shown
    When user clicks on "pivot table" icon in toolbox
    Then pivot table viewer should be visible
    When user adds "ShipCountry" to the "group by" row of pivot table viewer
    And user clicks on the "remove group by chip ShipName" area of pivot table viewer
    And user adds "ShipVia" to the "pivot" row of pivot table viewer
    And user clicks on the "remove pivot chip CustomerID" area of pivot table viewer
    And user picks "Aggregation > count" from the context menu of the "aggregate chip avg(OrderID)" area of pivot table viewer
    And user closes the context menu
    Then the "group by" reading of pivot table viewer should be "ShipCountry"
    And the "pivot" reading of pivot table viewer should be "ShipVia"
    And the "aggregate" reading of pivot table viewer should be "count(OrderID)"
    When user clicks on the "add to workspace" area of pivot table viewer
    Then the "orders aggregation" view should be current
    And the table should have 40 rows
    When user closes "orders" view tab
    Then the "orders aggregation" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "orders aggregation" project table in "Save project" dialog should be disabled
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
