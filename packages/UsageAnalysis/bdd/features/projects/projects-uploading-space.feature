@journey @serial @realizes:views.projects @realizes:views.space @realizes:data.menu.link-tables @realizes:file.menu.save.tables-as-project
Feature: Projects uploaded from files in a space, with Data sync on and off
  The setup of the md: a space is created in the Browse tree and customers.csv and orders.csv are
  copied into it from Browse > Files > Demo > northwind (drag onto the space, Copy). TestCase4 opens
  both files from the space, TestCase5 customers.csv from the space and orders.csv from Files >
  Demo > northwind. As in TestCase1 (projects-uploading), the two tables are linked
  selection-to-filter through Data > Link Tables..., the link is checked in the status bar of
  orders, and both are saved in one project through the ribbon's Save dialog with Data sync on and
  off; the project reopens from its Dashboards card with both tables and a link that still filters,
  and the Save dialog of the reopened project shows each table's Creation script only with Data
  sync. Translated from the TestTrack case Projects/uploading.

  Fixture numbers: customers.csv has 91 rows and orders.csv 830; the first two customers have 6 and
  4 orders, the next two 7 and 13. The link is checked with the first pair before the save and the
  second pair after the reopen (the md table says rows 1 and 2 again; a filter the project merely
  restored would pass those whatever the link does). Parked (see the request document): the proof that a reopened table
  is not a frame that stayed open. The console errors are not claimed (the Save dialog's preview
  noise, GROK-18606).

  Names are letters and digits only (the Dashboards search misses "-" and "_") and carry the run's
  time. The space (with its files) and the four projects (with their tables and views) are removed
  when the feature starts and ends. It is serial: the Dashboards search and the uploads are shared
  with every feature that saves a project.

  Background:
    Given user is logged in
    And simple mode is off
    And the browse panel is open
    And no space named "BDDUpSpace{time}" is on the server
    And no project named "BDDUpload4Sync{time}" is on the server
    And no project named "BDDUpload4NoSync{time}" is on the server
    And no project named "BDDUpload5Sync{time}" is on the server
    And no project named "BDDUpload5NoSync{time}" is on the server

  Scenario: A space holds copies of customers.csv and orders.csv
    When user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    Then Create Space dialog should be visible
    When user enters "BDDUpSpace{time}" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then the "Create Space" dialog should close
    And 1 space named "BDDUpSpace{time}" should be on the server
    Given Spaces tree node inside browse tree is expanded
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user clicks on "Files > Demo > northwind" tree node inside browse tree
    Then the "Demo/northwind" view should be current
    When user drags customers.csv link in gallery to BDDUpSpace{time} tree node inside browse tree
    Then Move entity dialog should be visible
    When user selects "Copy" in Move entity dialog
    And user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    When user drags orders.csv link in gallery to BDDUpSpace{time} tree node inside browse tree
    Then Move entity dialog should be visible
    When user selects "Copy" in Move entity dialog
    And user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    When user double-clicks on BDDUpSpace{time} tree node inside browse tree
    Then the "BDDUpSpace{time}" view should be current
    And customers.csv link in gallery should be visible
    And orders.csv link in gallery should be visible

  Scenario Outline: customers.csv from the space and orders.csv from <Source>, linked and saved with Data sync <Sync>
    When user picks "Close All" from the context menu of browse tab
    Given the browse panel is open
    And Spaces tree node inside browse tree is expanded
    When user double-clicks on BDDUpSpace{time} tree node inside browse tree
    Then the "BDDUpSpace{time}" view should be current
    When user double-clicks on customers.csv link in gallery
    Then the "customers" view should be current
    And table "customers" should have 91 rows
    When user clicks on browse tab
    And user <open folder>
    Then the "<Folder view>" view should be current
    When user double-clicks on orders.csv link in gallery
    Then the "orders" view should be current
    And table "orders" should have 830 rows
    When user picks "Data > Link Tables..." from the top menu
    Then "Link Tables" dialog should be visible
    When user sets the tables of the Link Tables dialog to "customers" and "orders"
    And user sets key columns 1 of the Link Tables dialog to "CustomerID" and "CustomerID"
    And user selects "selection to filter" in "Link Type" input in "Link Tables" dialog
    And user clicks on LINK button in "Link Tables" dialog
    Then "Link Tables" dialog should contain text "customers -> orders"
    When user clicks on CLOSE button in "Link Tables" dialog
    Then the "Link Tables" dialog should close
    When user clicks on customers view
    And user clicks on the "row header 1" area of grid
    And user clicks on the "row header 2" area of grid holding Shift
    Then 2 rows of table "customers" should be selected
    When user clicks on orders view
    Then 10 rows of table "orders" should pass the filter
    And status bar should contain text "Filtered: 10"
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "<Project>" into Name text input in "Save project" dialog
    And user <switch> Data sync switch in "customers" project table in "Save project" dialog
    And user <switch> Data sync switch in "orders" project table in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "<Project>" uploaded' should have been shown
    And the "customers" table of the "<Project>" project should be saved <saved>
    And the "orders" table of the "<Project>" project should be saved <saved>
    And "Share <Project>" dialog should be visible
    When user clicks on CANCEL button in "Share <Project>" dialog
    Then the "Share <Project>" dialog should close
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "<Project>" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on <Project> gallery card
    Then the "customers" table view should open with 91 rows
    And the "orders" table view should open with 830 rows
    Given user switches to the "customers" table view
    Then the table should have been <how>
    Given user switches to the "orders" table view
    Then the table should have been <how>
    And no error or warning balloon should have been shown
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "customers" project table in "Save project" dialog should be switched <state>
    And "Creation script" button in "customers" project table in "Save project" dialog should be <script>
    And Data sync switch in "orders" project table in "Save project" dialog should be switched <state>
    And "Creation script" button in "orders" project table in "Save project" dialog should be <script>
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
    When user clicks on customers view
    And user clicks on the "row header 3" area of grid
    And user clicks on the "row header 4" area of grid holding Shift
    Then 2 rows of table "customers" should be selected
    When user clicks on orders view
    Then 20 rows of table "orders" should pass the filter
    And status bar should contain text "Filtered: 20"

    Examples:
      | Source                    | open folder                                                             | Folder view      | Sync | Project                | switch   | saved          | state | script  | how                   |
      | the space                 | double-clicks on BDDUpSpace{time} tree node inside browse tree          | BDDUpSpace{time} | ON   | BDDUpload4Sync{time}   | checks   | with data sync | on    | visible | reloaded by data sync |
      | the space                 | double-clicks on BDDUpSpace{time} tree node inside browse tree          | BDDUpSpace{time} | OFF  | BDDUpload4NoSync{time} | unchecks | as a snapshot  | off   | hidden  | loaded as a snapshot  |
      | Files > Demo > northwind  | clicks on "Files > Demo > northwind" tree node inside browse tree       | Demo/northwind   | ON   | BDDUpload5Sync{time}   | checks   | with data sync | on    | visible | reloaded by data sync |
      | Files > Demo > northwind  | clicks on "Files > Demo > northwind" tree node inside browse tree       | Demo/northwind   | OFF  | BDDUpload5NoSync{time} | unchecks | as a snapshot  | off   | hidden  | loaded as a snapshot  |
