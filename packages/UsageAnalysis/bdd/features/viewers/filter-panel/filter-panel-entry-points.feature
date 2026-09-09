@journey @viewers @realizes:viewers.filters
Feature: Filter panel entry points
  How cards get into the panel and out of it: the header's column picker and the panel's own
  Add Filter menu both insert at the top, a card's close icon removes it and releases whatever it
  was keeping, "Remove others" drops the cards that restrict nothing while "Remove All" empties the
  panel, and closing the panel releases its filtering while reopening it brings the card back with
  its criterion. One journey on demog-1000 (RACE: Caucasian 896, Black 27).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user opens an empty filter panel
    Then the filter panel should have 0 filters
    And all rows should pass the filter
    And counter of filter panel should be hidden

  Scenario: Cards added from the header picker stack at the top
    When user adds a card for "RACE" to the filter panel
    And user adds a card for "SEX" to the filter panel
    Then the filter panel should have 2 filters
    And the "cards" reading of filter panel should be "SEX, RACE"
    And "RACE" filter card should be visible
    And "SEX" filter card should be visible
    And all rows should pass the filter
    And counter of filter panel should be hidden
    And no errors should have been logged

  Scenario: A card added from the panel's own menu goes to the top too
    When user picks "Add Filter | Combined Boolean" from the filter panel menu
    Then "Flags" filter card should be visible
    And the filter panel should have 3 filters
    And the "cards" reading of filter panel should be "Flags, SEX, RACE"
    And all rows should pass the filter
    And counter of filter panel should be hidden
    And no errors should have been logged

  Scenario: A card's close icon removes it and releases its criterion
    When user clicks on the "category Caucasian of RACE" area of filter panel
    Then 896 rows should pass the filter
    And counter of filter panel should have text "1"
    When user hovers over "RACE" filter card
    And user clicks on close of "RACE" filter card
    Then "RACE" filter card should be absent
    And the filter panel should have 2 filters
    And all rows should pass the filter
    And counter of filter panel should be hidden
    And no errors should have been logged

  Scenario: Remove others keeps the cards that restrict rows, Remove All empties the panel
    When user adds a card for "RACE" to the filter panel
    And user clicks on the "category Caucasian of RACE" area of filter panel
    Then 896 rows should pass the filter
    And the filter panel should have 3 filters
    When user picks "Remove others" from the filter counter menu
    Then "RACE" filter card should be visible
    And "SEX" filter card should be absent
    And "Flags" filter card should be absent
    And the filter panel should have 1 filter
    And 896 rows should pass the filter
    When user picks "Remove All" from the filter panel menu
    Then the filter panel should have 0 filters
    And all rows should pass the filter
    And counter of filter panel should be hidden
    And no errors should have been logged

  Scenario: Closing the panel releases its filtering and reopening restores the card
    When user adds a card for "RACE" to the filter panel
    And user clicks on the "category Black of RACE" area of filter panel
    Then 27 rows should pass the filter
    When user clicks on close icon of filters viewer
    Then filter panel should be hidden
    And all rows should pass the filter
    When user opens the filter panel
    Then filter panel should be visible
    And "RACE" filter card should be visible
    And the "selected categories of RACE" reading of filter panel should be "Black"
    And 27 rows should pass the filter
    When user picks "Remove All" from the filter panel menu
    Then all rows should pass the filter
    And no errors should have been logged
