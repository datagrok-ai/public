@journey @viewers @realizes:viewers.filters
Feature: Combined boolean filter card
  One card for every boolean column of the table: a boolean column added later joins it, a click on
  a flag's name keeps that flag's true rows and drops whatever the card kept before, the state
  survives a layout round-trip, and removing the card releases the rows while Add Filter brings it
  back exactly once (GROK-16488). One journey on demog-1000 with its own CONTROL (true in 6 rows)
  and a calculated SEX_bool (true in the 553 F rows).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user opens an empty filter panel
    And user picks "Add Filter | Combined Boolean" from the filter panel menu
    Then "Flags" filter card should be visible
    And the filter panel should have 1 filter
    And the "categories of Flags" reading of filter panel should be "CONTROL"
    And all rows should pass the filter
    And counter of filter panel should be hidden

  Scenario: A boolean column added later joins the card
    When user adds a calculated column "SEX_bool" with formula "${SEX} == \"F\""
    Then "SEX_bool" column should have type "bool"
    And the "categories of Flags" reading of filter panel should be "CONTROL, SEX_bool"
    And the "cards" reading of filter panel should be "Flags"
    And the filter panel should have 1 filter
    And all rows should pass the filter
    And no errors should have been logged

  Scenario: A click on a flag keeps its true rows
    When user clicks on the "category CONTROL of Flags" area of filter panel
    Then 6 rows should pass the filter
    And the filter should pass exactly the rows where "CONTROL" is "true"
    And the "summary of Flags" reading of filter panel should be "CONTROL: True"
    And the "filters" reading of filter panel should be 1
    And counter of filter panel should have text "1"
    And no errors should have been logged

  Scenario: A click on another flag replaces the one the card kept
    When user clicks on the "category SEX_bool of Flags" area of filter panel
    Then 553 rows should pass the filter
    And the filter should pass exactly the rows where "SEX" is "F"
    And the "summary of Flags" reading of filter panel should be "SEX_bool: True"
    And counter of filter panel should have text "1"
    And no errors should have been logged

  Scenario: A saved layout brings the card and its flag back
    When user saves the layout of the current table view to the server
    And user picks "Remove All" from the filter panel menu
    Then the filter panel should have 0 filters
    And all rows should pass the filter
    When user loads the saved layout
    Then "Flags" filter card should be visible
    And the filter panel should have 1 filter
    And 553 rows should pass the filter
    And the "summary of Flags" reading of filter panel should be "SEX_bool: True"
    And counter of filter panel should have text "1"
    And no errors should have been logged

  Scenario: Removing the card releases the rows and Add Filter brings it back once
    When user hovers over "Flags" filter card
    And user clicks on close of "Flags" filter card
    Then "Flags" filter card should be absent
    And the filter panel should have 0 filters
    And all rows should pass the filter
    And counter of filter panel should be hidden
    When user picks "Add Filter | Combined Boolean" from the filter panel menu
    Then "Flags" filter card should be visible
    And the "cards" reading of filter panel should be "Flags"
    And the filter panel should have 1 filter
    And all rows should pass the filter
    When user picks "Remove All" from the filter panel menu
    And user removes "SEX_bool" column
    Then the table should have 11 columns
    And all rows should pass the filter
    And no errors should have been logged
