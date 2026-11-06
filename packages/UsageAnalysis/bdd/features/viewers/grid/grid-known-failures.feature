@viewers @realizes:viewers.grid
Feature: Grid known failures
  Temporary: the known failures of the grid while their tickets are open. Each scenario is expected
  to fail on its pinned step; when its ticket is fixed, remove the tag and move the scenario back
  into its topical feature as an ordinary @edge-cases scenario: core#1198 into
  grid-appearance.feature, core#1195 into grid-column-groups.feature, core#1193 into
  grid-editing.feature, core#1194 into grid-rows.feature. The known failures of grid-pinning.feature
  and grid-forms-column.feature are not here: those features are not journeys.

  Each scenario is a test of its own on demog-1000.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    Then grid should show 1000 rows

  @edge-cases @known-failure:22 @realizes:core#1198
  Scenario: A malformed number format marks the format field invalid
    When user picks "Format > Custom..." from the context menu of the "header AGE" area of grid
    Then Format AGE dialog should be visible
    When user enters "0.0.0" into Custom input in Format AGE dialog
    Then Custom input in Format AGE dialog should be invalid

  @edge-cases @known-failure:30 @realizes:core#1198
  Scenario: A malformed number format is not taken and the old format stays
    When user picks "Format > Custom..." from the context menu of the "header AGE" area of grid
    Then Format AGE dialog should be visible
    When user enters "0.0.0" into Custom input in Format AGE dialog
    And user clicks on OK button in Format AGE dialog
    Then "AGE" column should have tag "format" equal to ""

  @edge-cases @known-failure:42 @realizes:core#1195
  Scenario: A column group named with spaces only is not made
    When user clicks on the "header AGE" area of grid holding Control
    And user clicks on the "header SEX" area of grid holding Control
    Given the context panel is open
    And Actions accordion header in context panel is expanded
    When user clicks on "Group columns..." text in context panel
    Then Group input in dialog should be visible
    When user types "   " into Group input in dialog
    And user clicks on OK button in dialog
    Then "AGE" column should have tag "group" equal to ""

  @edge-cases @known-failure:50 @realizes:core#1193
  Scenario: A read-only grid keeps its rows on Shift+Delete
    When user sets "Allow Edit" property of grid to "false"
    And user clicks on the "row header 5" area of grid holding Control
    Then 1 row should be selected
    When user presses Shift+Delete
    Then the table should have 1000 rows

  @edge-cases @known-failure:61 @realizes:core#1194
  Scenario: Shift+Enter selects nothing while Allow Row Selection is off
    When user clicks on the "cell 4 of SEX" area of grid
    And user presses Shift+Enter
    Then 447 rows should be selected
    When user clears the row selection
    And user sets "Allow Row Selection" property of grid to "false"
    And user clicks on the "cell 1 of SEX" area of grid
    And user presses Shift+Enter
    Then no rows should be selected

  @edge-cases @known-failure:70 @realizes:core#1194
  Scenario: Shift+Enter selects only the same rows that pass the filter
    When user filters rows where "AGE" is between 20 and 30
    And user clicks on the "cell 2 of SEX" area of grid
    And user clicks on the "cell 1 of SEX" area of grid
    And user presses Shift+Enter
    Then some rows should be selected
    And 82 rows should be selected
