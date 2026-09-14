@journey @viewers @realizes:viewers.filters
Feature: Filter panel of a cloned view
  View > Layout > Clone View opens a second view over the same table, and its filter panel comes up
  with the original's cards, criteria and switches; a category clicked on a card of the original is
  followed by the clone's card on the same column, a card switched off in one view is switched off
  in the other while a card of another type on the same column is left alone, and removing a card
  in the clone leaves the original's card as it was. A layout saved from the clone brings its cards
  and switches back. A missing-values choice made from a card's own menu survives the clone, and
  the clone refuses a second card on that column (github-1984). A second table built on the very
  column object of the first is not filtered by the first table's card (GROK-13582). One journey on
  demog-1000: RACE is Caucasian 896, Black 27, Other 62 (985 together); 18 of the Black rows are F;
  HEIGHT is blank in 128 rows.
  Not translated: the md's AGE for the missing-values choice — demog-1000 has no blank AGE, so the
  choice is made on HEIGHT, the one numeric column with blanks, from the card's menu as D8 requires;
  the view switches are the shell's (the view tabs are hidden in the simple mode a bdd page runs
  in), the gestures in each view are made on that view's own panel; the second card is offered by
  dragging the column's header (the picker, refusing it, leaves its own grid pending), after a
  drag that does add a card; the "not filtered" reads are single reads, not a held window; the
  row count after the clone's layout comes back is not claimed — it is the two panels' together
  over one table, and the original view's panel is not reset by the clone's layout (measured on dev:
  11 rows after the reload against 380 at the save).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user opens an empty filter panel
    And user adds a card for "RACE" to the filter panel
    And user adds a card for "SEX" to the filter panel
    And user adds a categorical filter on "RACE" keeping "Caucasian, Black, Other"
    Then 985 rows should pass the filter
    And counter of filter panel should have text "1"

  Scenario: A cloned view comes up with the same cards, criteria and switches
    When user picks "View > Layout > Clone View" from the top menu
    Then the "demog-1000 copy" view should be current
    And filter panel should be visible
    And the "cards" reading of filter panel should be "SEX, RACE"
    And the "selected categories of RACE" reading of filter panel should be "Black, Caucasian, Other"
    And the "filtering of RACE" reading of filter panel should be "true"
    And the "filtering of SEX" reading of filter panel should be "false"
    And "RACE" filter card should be enabled
    And "SEX" filter card should be enabled
    And 985 rows should pass the filter
    And no errors should have been logged

  Scenario: A category clicked on the original's card is followed by the clone's card
    When user switches to the "demog-1000" view
    And user clicks on the "category Black of RACE" area of filter panel
    Then 27 rows should pass the filter
    And the "selected categories of RACE" reading of filter panel should be "Black"
    When user switches to the "demog-1000 copy" view
    Then the "selected categories of RACE" reading of filter panel should be "Black"
    And the "filtering of RACE" reading of filter panel should be "true"
    And the "filtering of SEX" reading of filter panel should be "false"
    And no errors should have been logged

  Scenario: A card switched off in one view is switched off in the other
    When user switches to the "demog-1000" view
    And user clicks on the "category F of SEX" area of filter panel
    Then 18 rows should pass the filter
    When user hovers over "SEX" filter card
    And user unchecks checkbox of "SEX" filter card
    Then 27 rows should pass the filter
    When user switches to the "demog-1000 copy" view
    Then "SEX" filter card should be disabled
    And the "enabled of SEX" reading of filter panel should be "false"
    And "RACE" filter card should be enabled
    When user switches to the "demog-1000" view
    And user switches the "SEX" filter card back on
    Then 18 rows should pass the filter
    When user switches to the "demog-1000 copy" view
    Then "SEX" filter card should be enabled
    And no errors should have been logged

  Scenario: Switching off a card of another type in the clone leaves the original's card alone
    When user switches to the "demog-1000" view
    And user adds a range filter on "AGE" from 30 to 60
    Then fewer than 18 rows should pass the filter
    And the "rows shown" reading of filter panel should be at least 1
    When user remembers the "rows shown" reading of filter panel
    And user switches to the "demog-1000 copy" view
    And user adds a card for "AGE" to the filter panel
    And user hovers over "AGE" filter card
    And user clicks on "Switch to categorical filter" icon in "AGE" filter card
    Then the "type of AGE" reading of filter panel should be "categorical"
    When user hovers over "AGE" filter card
    And user unchecks checkbox of "AGE" filter card
    Then "AGE" filter card should be disabled
    And the "rows shown" reading of filter panel should be as remembered
    When user switches to the "demog-1000" view
    Then "AGE" filter card should be enabled
    And the "type of AGE" reading of filter panel should be "histogram"
    And the "rows shown" reading of filter panel should be as remembered
    And no errors should have been logged

  Scenario: Removing a card in the clone leaves the original's card as it was
    When user switches to the "demog-1000 copy" view
    And user hovers over "SEX" filter card
    And user clicks on close of "SEX" filter card
    Then "SEX" filter card should be absent
    And the "rows shown" reading of filter panel should be as remembered
    When user switches to the "demog-1000" view
    Then "SEX" filter card should be visible
    And "SEX" filter card should be enabled
    And the "selected categories of SEX" reading of filter panel should be "F"
    And the "rows shown" reading of filter panel should be as remembered
    And no errors should have been logged

  Scenario: A layout saved from the clone brings its cards and switches back
    When user switches to the "demog-1000 copy" view
    And user hovers over "RACE" filter card
    And user unchecks checkbox of "RACE" filter card
    Then "RACE" filter card should be disabled
    When user remembers the "rows shown" reading of filter panel
    And user remembers the "cards" reading of filter panel
    And user saves the layout of the current table view to the server
    And user switches the "RACE" filter card back on
    Then the "rows shown" reading of filter panel should not be as remembered
    When user clicks on close icon of filters viewer
    Then filter panel should be hidden
    When user loads the saved layout
    Then filter panel should be visible
    And "RACE" filter card should be disabled
    And the "enabled of RACE" reading of filter panel should be "false"
    And the "enabled of AGE" reading of filter panel should be "false"
    And the "cards" reading of filter panel should be as remembered
    And no errors should have been logged

  Scenario: A missing-values choice made from the card's menu survives the clone
    When user opens demog-1000 dataset keeping the first 1000 rows as "demog-missing"
    And user clicks on filter icon in toolbar
    Then "HEIGHT" filter card should be visible
    And all rows should pass the filter
    When user opens the indicator menu of the "HEIGHT" filter card
    Then the open menu should list "Missing values > Keep missing values"
    And "Keep missing values" menu item should be selected
    And "Filter out missing values" menu item should not be selected
    And "Show only missing value" menu item should not be selected
    When user closes the context menu
    And user picks "Missing values | Filter out missing values" from the indicator menu of the "HEIGHT" filter card
    And user closes the context menu
    Then 872 rows should pass the filter
    And the "filtering of HEIGHT" reading of filter panel should be "true"
    When user picks "View > Layout > Clone View" from the top menu
    Then the "demog-missing copy" view should be current
    And 872 rows should pass the filter
    When user opens the indicator menu of the "HEIGHT" filter card
    Then the open menu should list "Missing values > Filter out missing values"
    And "Filter out missing values" menu item should be selected
    And "Keep missing values" menu item should not be selected
    And "Show only missing value" menu item should not be selected
    When user closes the context menu
    And user remembers the "cards" reading of filter panel
    And user drags the "header USUBJID" area of grid onto the "view" area of filter panel
    Then the "cards" reading of filter panel should not be as remembered
    And "USUBJID" filter card should be visible
    When user remembers the "cards" reading of filter panel
    And user drags the "header HEIGHT" area of grid onto the "view" area of filter panel
    Then the "cards" reading of filter panel should be as remembered
    When user switches to the "demog-missing" view
    Then 872 rows should pass the filter
    When user opens the indicator menu of the "HEIGHT" filter card
    Then the open menu should list "Missing values > Filter out missing values"
    And "Filter out missing values" menu item should be selected
    And "Keep missing values" menu item should not be selected
    And "Show only missing value" menu item should not be selected
    When user closes the context menu
    Then 872 rows should pass the filter
    And no errors should have been logged

  Scenario: A table built on the same column object is not filtered by the first table's card
    When user opens a table "shared SEX" that shares the "SEX" column of the current table
    Then 1000 rows of table "shared SEX" should pass the filter
    When user switches to the "demog-missing" view
    And user clicks on the "category F of SEX" area of filter panel
    Then fewer than 553 rows should pass the filter
    And the "rows shown" reading of filter panel should be at least 1
    And no rows where "SEX" is "M" should pass the filter
    And 1000 rows of table "shared SEX" should pass the filter
    And no errors should have been logged
