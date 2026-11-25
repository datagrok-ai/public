@realizes:views.projects
Feature: Projects regressions: a derived table saved after its source table is closed
  core#531: orders.csv from Browse > Files > Demo > northwind gets a pivot table (Group by
  ShipCountry, Pivot ShipVia, count(OrderID)) published into the workspace with ADD, and the orders
  view is then closed. The Save dialog offered Data sync for the published table, whose creation
  script aggregates the closed "orders", and a project saved with it on could not be opened again: the
  reopen stopped on "Data loading error" with "Could not resolve table "orders"". The Save dialog now
  disables Data sync for a table whose creation script reads a table that is not saved with it, and
  saving it as a snapshot reopens with its rows.

  Closing the source asks first (Close orders?, core#456): the scenarios answer CLOSE, since
  the subject here is the Save dialog that follows. The same rule holds for a join, which reads two
  tables, and along a chain: a table derived from a derived table keeps Data sync when only the first
  source is closed, because the middle table is saved as a snapshot and loads before it.
  "Ask to save on close" is off for the feature (its default), since the pivot added to the source view would ask too.

  The pivot is configured through its own rows and chips, as in projects-uploading. The project is
  named with the run's time and removed with its table and view when the scenario starts and when the
  feature ends.

  Background:
    Given user is logged in
    And the "Ask to save on close" setting is off until the feature ends
    And simple mode is off
    And the browse panel is open

  @realizes:core#531
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
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "orders aggregation" project table in "Save project" dialog should be enabled
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
    When user closes "orders" view tab
    Then "Close orders?" dialog should be visible
    And "Close orders?" dialog should contain text 'Made from "orders"'
    When user clicks on CLOSE button in "Close orders?" dialog
    Then the "Close orders?" dialog should close
    And the "orders aggregation" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "orders aggregation" project table in "Save project" dialog should be disabled
    And Data sync switch in "orders aggregation" project table in "Save project" dialog should be unchecked
    When user enters "BDDRegClosedSrc{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegClosedSrc{time}" dialog should be visible
    And the "orders aggregation" table of the "BDDRegClosedSrc{time}" project should be saved as a snapshot
    When user clicks on CANCEL button in "Share BDDRegClosedSrc{time}" dialog
    Then the "Share BDDRegClosedSrc{time}" dialog should close
    When user closes all views
    And user opens the "BDDRegClosedSrc{time}" project and waits for its table
    Then the "orders aggregation" table view should open with 40 rows
    And "Data loading error" dialog should be absent

  Scenario: Data sync is not offered for a join one of whose tables is closed
    Given user opens cars dataset
    And user opens cars dataset keeping the first 10 rows as "cars10"
    When user picks "Data > Join Tables..." from the top menu
    Then "Join Tables" dialog should be visible
    When user selects "cars" in Tables input in "Join Tables" dialog
    And user selects "cars10" in second choice input in "Join Tables" dialog
    And user clicks on OK button in "Join Tables" dialog
    Then the "Join Tables" dialog should close
    And the "result" view should be current
    When user closes "cars10" view tab
    Then "Close cars10?" dialog should be visible
    And "Close cars10?" dialog should contain text "result"
    And "Close cars10?" dialog should contain text 'Made from "cars10"'
    When user clicks on CLOSE button in "Close cars10?" dialog
    Then the "Close cars10?" dialog should close
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "result" project table in "Save project" dialog should be disabled
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close

  Scenario: Along a chain of derived tables only the table whose source is closed loses Data sync
    Given user clears the saved pivot table parameters
    And user opens cars dataset
    And the toolbox pane is shown
    When user clicks on "pivot table" icon in toolbox
    And user clicks on the "add to workspace" area of pivot table viewer
    Then the "cars aggregation" view should be current
    When user clicks on "pivot table" icon in toolbox
    And user clicks on the "add to workspace" area of pivot table viewer
    Then the "cars aggregation aggregation" view should be current
    When user closes "cars" view tab
    Then "Close cars?" dialog should be visible
    When user clicks on CLOSE button in "Close cars?" dialog
    Then the "Close cars?" dialog should close
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "cars aggregation" project table in "Save project" dialog should be disabled
    And Data sync switch in "cars aggregation aggregation" project table in "Save project" dialog should be enabled
    When user hovers over the text "Data sync" in "cars aggregation" project table in "Save project" dialog
    Then tooltip should contain text 'Data sync is off: this table is made from "cars", which this project does not save'
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
