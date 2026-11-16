@realizes:views.projects @realizes:core#456
Feature: Closing a table that other open tables depend on
  core#456 and core#531: closing a table that something else open still reads asks first. Four kinds
  of dependency are named: a table derived from it (a pivot published to the workspace, a join; its
  creation script reads the table by name), a formula column of another table that reads it
  (RowCount("cars")), a link between it and another table, and a viewer of another table's view
  bound to it. The "Close table" dialog lists what is affected; CANCEL leaves the tab, the table and
  every dependent as they were, CLOSE ANYWAY closes the table and leaves the dependents to what follows:
  the published pivot keeps its rows, the formula column empties, the link is gone, the viewer keeps
  the rows it showed. The last open view of a saved project closes the project with all of its tables,
  so it asks about what another open project reads from them (a table opened while a project is open
  joins that project), though a view of a saved project otherwise leaves its table open.

  It asks only where something is lost: never for a table nothing reads, for the dependent itself, for
  one of two views of a table (the table stays open), for a dependent closed before its source, or for
  Close All, which closes everything at once. Close others asks once for everything it closes, and not
  at all when it closes a source together with all its dependents. The table's Close command and the
  simple mode's close icon ask the same way as the tab's close button.

  Every table is opened through the API (cars, and its first 10 rows as cars10); the pivot is
  published through the toolbox. Only the two-project scenario saves anything: its projects are
  named with the run's time and removed when it starts and when the feature ends.
  "Warn on unsaved changes" is off for the feature: a pivot added to the source view is unsaved work,
  whose own question (projects-close-unsaved) would come first.

  Background:
    Given user is logged in
    And the "Warn on unsaved changes" setting is off until the feature ends

  Scenario: A table a published pivot was made from asks before it closes
    Given simple mode is off
    And user clears the saved pivot table parameters
    And user opens cars dataset
    And the toolbox pane is shown
    When user clicks on "pivot table" icon in toolbox
    And user clicks on the "add to workspace" area of pivot table viewer
    Then the "cars aggregation" view should be current
    When user closes "cars" view tab
    Then "Close table" dialog should be visible
    And "Close table" dialog should contain text "cars aggregation"
    When user clicks on CANCEL button in "Close table" dialog
    Then the "Close table" dialog should close
    And "cars" view tab should be visible
    And table "cars" should be open
    When user closes "cars" view tab
    And user clicks on CLOSE ANYWAY button in "Close table" dialog
    Then the "Close table" dialog should close
    And "cars" view tab should be absent
    And the "cars aggregation" view should be current
    And table "cars aggregation" should be open

  Scenario: A table a formula column of another table reads asks before it closes
    Given simple mode is off
    And user opens cars dataset
    And user opens cars dataset keeping the first 10 rows as "cars10"
    And user adds a calculated column "n" with formula 'RowCount("cars")'
    Then table "cars10" should have no missing values in "n" column
    When user closes "cars" view tab
    Then "Close table" dialog should be visible
    And "Close table" dialog should contain text "cars10"
    When user clicks on CANCEL button in "Close table" dialog
    Then the "Close table" dialog should close
    And "cars" view tab should be visible
    And table "cars10" should have no missing values in "n" column
    When user closes "cars" view tab
    And user clicks on CLOSE ANYWAY button in "Close table" dialog
    Then the "Close table" dialog should close
    And "cars" view tab should be absent
    And table "cars10" should have missing values in "n" column

  Scenario: A table linked to another table asks before it closes
    Given simple mode is off
    And user opens cars dataset
    And user opens cars dataset keeping the first 10 rows as "cars10"
    And tables "cars10" and "cars" are linked on "model"
    When user makes row 3 of table "cars10" current
    Then 1 rows of table "cars" should pass the filter
    When user closes "cars" view tab
    Then "Close table" dialog should be visible
    And "Close table" dialog should contain text "cars10"
    When user clicks on CANCEL button in "Close table" dialog
    Then the "Close table" dialog should close
    And "cars" view tab should be visible
    When user makes row 5 of table "cars10" current
    Then 1 rows of table "cars" should pass the filter
    When user closes "cars" view tab
    And user clicks on CLOSE ANYWAY button in "Close table" dialog
    Then the "Close table" dialog should close
    And "cars" view tab should be absent

  Scenario: A table a viewer of another view is bound to asks before it closes
    Given simple mode is off
    And user opens cars dataset
    And user opens cars dataset keeping the first 10 rows as "cars10"
    And user adds a scatter plot viewer with:
      | Table | cars |
    Then "Table" property of scatter plot viewer should be "cars"
    When user closes "cars" view tab
    Then "Close table" dialog should be visible
    And "Close table" dialog should contain text "cars10"
    When user clicks on CANCEL button in "Close table" dialog
    Then the "Close table" dialog should close
    And "cars" view tab should be visible
    And "Table" property of scatter plot viewer should be "cars"
    When user closes "cars" view tab
    And user clicks on CLOSE ANYWAY button in "Close table" dialog
    Then the "Close table" dialog should close
    And "cars" view tab should be absent
    And "Table" property of scatter plot viewer should be "cars"

  Scenario: A table that several things depend on asks once, naming all of them
    Given simple mode is off
    And user clears the saved pivot table parameters
    And user opens cars dataset
    And the toolbox pane is shown
    When user clicks on "pivot table" icon in toolbox
    And user clicks on the "add to workspace" area of pivot table viewer
    Then the "cars aggregation" view should be current
    Given user opens cars dataset keeping the first 10 rows as "cars10"
    And user adds a calculated column "n" with formula 'RowCount("cars")'
    And tables "cars10" and "cars" are linked on "model"
    When user closes "cars" view tab
    Then "Close table" dialog should be visible
    And "Close table" dialog should contain text "cars aggregation"
    And "Close table" dialog should contain text "cars10"
    When user clicks on CANCEL button in "Close table" dialog
    Then the "Close table" dialog should close
    And "cars" view tab should be visible

  Scenario: Closing the last view of a saved project asks about what another project reads from it
    Given simple mode is off
    And no project named "BDDCloseDepSource{time}" is on the server
    And no project named "BDDCloseDepReader{time}" is on the server
    And user opens cars dataset
    When user saves the current view as project "BDDCloseDepSource{time}"
    And user closes all views
    Given user opens cars dataset keeping the first 10 rows as "cars10"
    And user adds a calculated column "n" with formula 'RowCount("cars")'
    When user saves the current view as project "BDDCloseDepReader{time}"
    And user closes all views
    And user opens the "BDDCloseDepReader{time}" project and waits for its table
    And user opens the "BDDCloseDepSource{time}" project and waits for its table
    Then table "cars10" should have no missing values in "n" column
    When user closes "cars" view tab
    Then "Close table" dialog should be visible
    And "Close table" dialog should contain text "cars10"
    When user clicks on CANCEL button in "Close table" dialog
    Then the "Close table" dialog should close
    And "cars" view tab should be visible
    And the project "BDDCloseDepSource{time}" should be open

  Scenario: Nothing asks when nothing depends on what closes
    Given simple mode is off
    And user clears the saved pivot table parameters
    And user opens cars dataset
    When user closes "cars" view tab
    Then "cars" view tab should be absent
    And "Close table" dialog should be absent
    Given user opens cars dataset
    And the toolbox pane is shown
    When user clicks on "pivot table" icon in toolbox
    And user clicks on the "add to workspace" area of pivot table viewer
    Then the "cars aggregation" view should be current
    When user closes "cars aggregation" view tab
    Then "cars aggregation" view tab should be absent
    And "Close table" dialog should be absent
    When user closes "cars" view tab
    Then no table should be left in the workspace
    And "Close table" dialog should be absent

  Scenario: Closing one of two views of a table does not ask
    Given simple mode is off
    And user opens cars dataset
    And user opens a second view of table "cars" named "cars copy"
    And user opens cars dataset keeping the first 10 rows as "cars10"
    And user adds a calculated column "n" with formula 'RowCount("cars")'
    When user closes "cars copy" view tab
    Then "cars copy" view tab should be absent
    And "Close table" dialog should be absent
    And table "cars" should be open
    And table "cars10" should have no missing values in "n" column

  Scenario: Close All closes dependent tables without asking
    Given simple mode is off
    And user opens cars dataset
    And user opens cars dataset keeping the first 10 rows as "cars10"
    And user adds a calculated column "n" with formula 'RowCount("cars")'
    When user picks "Close All" from the context menu of the sidebar
    Then no table should be left in the workspace
    And "Close table" dialog should be absent

  Scenario: Close others asks once for the tables it closes
    Given simple mode is off
    And user clears the saved pivot table parameters
    And user opens cars dataset
    And the toolbox pane is shown
    When user clicks on "pivot table" icon in toolbox
    And user clicks on the "add to workspace" area of pivot table viewer
    Then the "cars aggregation" view should be current
    When user picks "Close others" from the context menu of the current view tab
    Then "Close table" dialog should be visible
    And "Close table" dialog should contain text "cars aggregation"
    When user clicks on CANCEL button in "Close table" dialog
    Then the "Close table" dialog should close
    And "cars" view tab should be visible
    When user picks "Close others" from the context menu of the current view tab
    And user clicks on CLOSE ANYWAY button in "Close table" dialog
    Then the "Close table" dialog should close
    And "cars" view tab should be absent
    And the "cars aggregation" view should be current

  Scenario: Close others does not ask when a table closes together with what depends on it
    Given simple mode is off
    And user clears the saved pivot table parameters
    And user opens cars dataset
    And the toolbox pane is shown
    When user clicks on "pivot table" icon in toolbox
    And user clicks on the "add to workspace" area of pivot table viewer
    Then the "cars aggregation" view should be current
    Given user opens cars dataset keeping the first 10 rows as "cars10"
    When user picks "Close others" from the context menu of the current view tab
    Then the open table views should be exactly "cars10"
    And "Close table" dialog should be absent

  Scenario: The table's Close command asks like the tab's close button
    Given simple mode is off
    And user opens cars dataset
    And user opens cars dataset keeping the first 10 rows as "cars10"
    And user adds a calculated column "n" with formula 'RowCount("cars")'
    When user clicks on the tab of the "cars" view
    And user picks "Close" from the context menu of the table name in the status bar
    Then "Close table" dialog should be visible
    And "Close table" dialog should contain text "cars10"
    When user clicks on CANCEL button in "Close table" dialog
    Then the "Close table" dialog should close
    And table "cars" should be open

  Scenario: The simple mode's close icon asks like the tab's close button
    Given user opens cars dataset
    And user opens cars dataset keeping the first 10 rows as "cars10"
    And user adds a calculated column "n" with formula 'RowCount("cars")'
    When user closes the "cars" view from the view selector
    Then "Close table" dialog should be visible
    And "Close table" dialog should contain text "cars10"
    When user clicks on CANCEL button in "Close table" dialog
    Then the "Close table" dialog should close
    And table "cars" should be open
