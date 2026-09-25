@journey @realizes:powerpack.cp.add-new-column-persists @realizes:powerpack.int.add-new-column-datasync-roundtrip @realizes:GROK-17109
Feature: Calculated columns over a file in My files follow a rename and an edit, and survive a project round trip
  A copy of demog in the user's own My files share, opened from Browse > Files > My files. Weight2 =
  ${WEIGHT} + 100 and Weight3 = ${Weight2} + 100 are added through the Add New Column dialog. WEIGHT
  is renamed through its header's Column Properties dialog and one of its cells edited in the grid:
  Weight2's formula follows the new name and both columns recalculate. The view is saved as a project
  with Data sync on, closed and reopened: the table is read from the file again, so the edited cell
  holds the file's value, while the rename and the two columns come back with their formulas and
  values that follow the source; a second rename and a second edit on the reopened table are followed the same way
  (GROK-17109). Translated from TestTrack PowerPack/add-new-column-advanced.md, the Home dir source.

  The views are closed through the platform's API (the workspace then holds no table), and the
  project is reopened from Browse > Dashboards.

  Background:
    Given user is logged in
    And no project named "bdd-anc-home-{run}" is on the server
    And a copy of the "System:DemoFiles/demog.csv" file is in the home folder as "bdd-anc-{run}.csv"
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---My-files tree node inside browse tree is expanded
    When user double-clicks on Files---My-files---bdd-anc-{run}.csv tree node inside browse tree
    Then the "bdd-anc-{run}" view should be current
    And the table should have 5850 rows

  Scenario: Two chained columns, then a rename and an edit of their source
    When user clicks on "Add New Column..." icon
    And user types "Weight2" into column name input
    And user types "${WEIGHT} + 100" into formula editor
    And user clicks on OK button in "Add New Column" dialog
    Then "Add New Column" dialog should be hidden
    When user clicks on "Add New Column..." icon
    And user types "Weight3" into column name input
    And user types "${Weight2} + 100" into formula editor
    And user clicks on OK button in "Add New Column" dialog
    Then "Add New Column" dialog should be hidden
    And every value of "Weight2" column should equal "WEIGHT" column plus 100
    And every value of "Weight3" column should equal "Weight2" column plus 100
    When user picks "Column Properties..." from the context menu of the "header WEIGHT" area of grid
    And user types "BaseWeight" into "New name:" input in "WEIGHT" dialog
    And user clicks on OK button in "WEIGHT" dialog
    Then the table should have a column "BaseWeight"
    And "Weight2" column should have tag "formula" equal to "${BaseWeight} + 100"
    When user double-clicks on the "cell 1 of BaseWeight" area of grid
    And user presses Control+A in cell editor
    And user types "500" into cell editor
    And user presses Enter
    Then the value of "BaseWeight" column in row 1 should be "500"
    And the value of "Weight2" column in row 1 should be "600"
    And every value of "Weight2" column should equal "BaseWeight" column plus 100
    And every value of "Weight3" column should equal "Weight2" column plus 100
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Saved with Data sync, closed and reopened, the columns keep their formulas
    When user clicks on Save button in toolbar
    Then "Save project" dialog should be visible
    When user enters "bdd-anc-home-{run}" into Name text input in "Save project" dialog
    And user switches on Data sync input in "Save project" dialog
    Then Data sync input in "Save project" dialog should be switched on
    When user clicks on OK button in "Save project" dialog
    Then "Save project" dialog should be hidden
    And no error or warning balloon should have been shown
    And 1 project named "bdd-anc-home-{run}" should be on the server
    When user presses Escape
    And user closes all views
    Then no table should be open
    When user clicks on Dashboards tree node inside browse tree
    Then the "Projects" view should be current
    When user types "bdd-anc-home-{run}" into gallery search
    Then "bdd-anc-home-{run}" project card should become visible within 60 seconds
    When user double-clicks on "bdd-anc-home-{run}" project card
    Then the "bdd-anc-{run}" view should be current
    And the table should have 5850 rows
    And the table should have a column "Weight2"
    And the table should have a column "Weight3"
    And the table should have a column "BaseWeight"
    And "Weight2" column should have tag "formula" equal to "${BaseWeight} + 100"
    And "Weight3" column should have tag "formula" equal to "${Weight2} + 100"
    # Data sync reads the file again: the edited cell holds what the file holds (73.2), the rename stays
    And the value of "BaseWeight" column in row 1 should be "73.19999694824219"
    And every value of "Weight2" column should equal "BaseWeight" column plus 100
    And every value of "Weight3" column should equal "Weight2" column plus 100
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: On the reopened table a second rename and edit are followed too
    When user picks "Column Properties..." from the context menu of the "header BaseWeight" area of grid
    And user types "BaseWeight2" into "New name:" input in "BaseWeight" dialog
    And user clicks on OK button in "BaseWeight" dialog
    Then the table should have a column "BaseWeight2"
    And "Weight2" column should have tag "formula" equal to "${BaseWeight2} + 100"
    And "Weight3" column should have tag "formula" equal to "${Weight2} + 100"
    When user double-clicks on the "cell 2 of BaseWeight2" area of grid
    And user presses Control+A in cell editor
    And user types "400" into cell editor
    And user presses Enter
    Then the value of "BaseWeight2" column in row 2 should be "400"
    And the value of "Weight3" column in row 2 should be "600"
    And every value of "Weight2" column should equal "BaseWeight2" column plus 100
    And every value of "Weight3" column should equal "Weight2" column plus 100
    And no error or warning balloon should have been shown
    And no errors should have been logged
