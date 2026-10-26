@viewers @realizes:viewers.filters
Feature: Filter panel regressions across layouts and cloned views
  Closed filter panel issues about what a saved layout or a cloned view brings back, each a claim
  that fails if its bug comes back: a card search reset before the layout was saved comes back as no
  filtering (public#3059); the master toggle switched off in one view is switched back on from the
  clone, and the cards filter again in both (public#1837, public#1460); a layout saved with the
  master toggle off comes back off, and switched on again every card has its checked box and its
  criterion (public#1103); a missing-values choice that came back with a layout can be changed back
  (public#791); a layout with unchecked categories brings the card's indicator and its Select all
  back (public#377); the order the cards were dragged into survives the layout (public#1464); the
  card's sort menu puts RA (434 rows) above Psoriasis, and the layout keeps that order (public#1889); a view
  cloned after Reset filters comes up unfiltered (public#319). Each scenario on its own demog-1000,
  where RACE is Black in 27 rows, DIS_POP is Psoriasis or PsA in 242 and HEIGHT is blank in 128.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user opens an empty filter panel
    Then the filter panel should have 0 filters
    And all rows should pass the filter

  Scenario: A layout saved after a card search was reset brings back no filtering (public#3059)
    When user adds a card for "DIS_POP" to the filter panel
    And user hovers over "DIS_POP" filter card
    And user clicks on search icon of "DIS_POP" filter card
    And user types "Ps" into the search of the "DIS_POP" filter card
    Then 242 rows should pass the filter
    And the "filtering of DIS_POP" reading of filter panel should be "true"
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    When user saves the layout of the current table view to the server
    And user clicks on the "category RA of DIS_POP" area of filter panel
    Then 434 rows should pass the filter
    When user loads the saved layout
    Then "DIS_POP" filter card should be visible
    And all rows should pass the filter
    And the "filtering of DIS_POP" reading of filter panel should be "false"
    And counter of filter panel should be hidden
    And no errors should have been logged

  Scenario: The master toggle switched off in one view is switched back on from the clone (public#1837, public#1460)
    When user adds a card for "RACE" to the filter panel
    And user clicks on the "category Black of RACE" area of filter panel
    Then 27 rows should pass the filter
    When user picks "View > Layout > Clone View" from the top menu
    Then the "demog-1000 copy" view should be current
    When user switches to the "demog-1000" view
    And user hovers over filter panel
    And user unchecks master of filter panel
    Then all rows should pass the filter
    And the "active" reading of filter panel should be "false"
    When user switches to the "demog-1000 copy" view
    And user hovers over filter panel
    And user checks master of filter panel
    Then 27 rows should pass the filter
    And the "active" reading of filter panel should be "true"
    When user switches to the "demog-1000" view
    Then 27 rows should pass the filter
    And the "active" reading of filter panel should be "true"
    When user switches to the "demog-1000 copy" view
    And user hovers over filter panel
    And user unchecks master of filter panel
    Then all rows should pass the filter
    And no errors should have been logged

  Scenario: A layout saved with the master toggle off comes back off, and switched on every card filters (public#1103)
    When user adds a card for "RACE" to the filter panel
    And user clicks on the "category Black of RACE" area of filter panel
    Then 27 rows should pass the filter
    When user hovers over filter panel
    And user unchecks master of filter panel
    Then all rows should pass the filter
    When user saves the layout of the current table view to the server
    And user loads the saved layout
    Then "RACE" filter card should be visible
    And the "active" reading of filter panel should be "false"
    And all rows should pass the filter
    When user hovers over filter panel
    And user checks master of filter panel
    Then 27 rows should pass the filter
    And the "active" reading of filter panel should be "true"
    And the "enabled of RACE" reading of filter panel should be "true"
    When user hovers over "RACE" filter card
    Then checkbox of "RACE" filter card should be checked
    And no errors should have been logged

  Scenario: A missing-values choice that came back with a layout can be undone (public#791)
    When user adds a card for "HEIGHT" to the filter panel
    And user picks "Missing values | Filter out missing values" from the indicator menu of the "HEIGHT" filter card
    And user closes the context menu
    Then 872 rows should pass the filter
    When user saves the layout of the current table view to the server
    And user loads the saved layout
    Then "HEIGHT" filter card should be visible
    And 872 rows should pass the filter
    When user picks "Missing values | Keep missing values" from the indicator menu of the "HEIGHT" filter card
    And user closes the context menu
    Then all rows should pass the filter
    When user opens the indicator menu of the "HEIGHT" filter card
    Then the open menu should list "Missing values > Keep missing values"
    And "Keep missing values" menu item should be selected
    And "Filter out missing values" menu item should not be selected
    When user closes the context menu
    Then no errors should have been logged

  Scenario: A layout with unchecked categories brings the card's indicator back (public#377)
    When user adds a card for "RACE" to the filter panel
    And user clicks on the "category Black of RACE" area of filter panel
    Then 27 rows should pass the filter
    When user saves the layout of the current table view to the server
    And user picks "Remove All" from the filter panel menu
    Then the filter panel should have 0 filters
    And all rows should pass the filter
    When user loads the saved layout
    Then "RACE" filter card should be visible
    And 27 rows should pass the filter
    When user hovers over caption of "RACE" filter card
    Then indicator of "RACE" filter card should be visible
    When user picks "Select all" from the indicator menu of the "RACE" filter card
    Then all rows should pass the filter
    And the "selected categories of RACE" reading of filter panel should be "Asian, Black, Caucasian, Other"
    And no errors should have been logged

  Scenario: The order the cards were dragged into comes back with the layout (public#1464)
    When user adds a card for "RACE" to the filter panel
    And user adds a card for "SEX" to the filter panel
    And user adds a card for "AGE" to the filter panel
    Then the "cards" reading of filter panel should be "AGE, SEX, RACE"
    When user drags the "RACE" filter card above the "SEX" filter card
    Then the "cards" reading of filter panel should be "AGE, RACE, SEX"
    When user saves the layout of the current table view to the server
    And user picks "Remove All" from the filter panel menu
    Then the filter panel should have 0 filters
    When user loads the saved layout
    Then the "cards" reading of filter panel should be "AGE, RACE, SEX"
    And no errors should have been logged

  Scenario: The card's sort menu orders the categories by count, and the layout keeps the order (public#1889)
    When user adds a card for "DIS_POP" to the filter panel
    Then the "category RA of DIS_POP" area of filter panel should lie below the "category Psoriasis of DIS_POP" area
    When user hovers over "DIS_POP" filter card
    And user clicks on "sort-amount-down" icon in "DIS_POP" filter card
    And user picks "Count (high to low)" from the open menu
    Then the "category RA of DIS_POP" area of filter panel should lie above the "category Psoriasis of DIS_POP" area
    When user saves the layout of the current table view to the server
    And user picks "Remove All" from the filter panel menu
    Then the filter panel should have 0 filters
    When user loads the saved layout
    Then "DIS_POP" filter card should be visible
    And the "category RA of DIS_POP" area of filter panel should lie above the "category Psoriasis of DIS_POP" area
    And no errors should have been logged

  Scenario: A view cloned after Reset filters comes up unfiltered (public#319)
    When user adds a card for "RACE" to the filter panel
    And user clicks on the "category Black of RACE" area of filter panel
    Then 27 rows should pass the filter
    When user clicks on close icon of filters viewer
    Then filter panel should be hidden
    And all rows should pass the filter
    When user clicks on filter icon in toolbar
    Then filter panel should be visible
    And 27 rows should pass the filter
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    When user picks "View > Layout > Clone View" from the top menu
    Then the "demog-1000 copy" view should be current
    And filter panel should be visible
    And "RACE" filter card should be visible
    And all rows should pass the filter
    And the "filtering of RACE" reading of filter panel should be "false"
    And no errors should have been logged
