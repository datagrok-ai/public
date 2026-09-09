@journey @viewers @realizes:viewers.filters
Feature: Categorical filter card
  What the categorical card does with its categories: a click on a name keeps one while a click on
  a checkbox adds one, the indicator menu's batch operations select, deselect and invert them all,
  Radio mode keeps exactly one and drops the batch items from the menu, the in-card search narrows
  the table to the categories that match, and a numeric column's card switches to categorical and
  back without filtering anything. One journey on demog-1000, whose DIS_POP is RA 434, Psoriasis
  204, Indigestion 152, UC 123, AS 49 and PsA 38.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user opens an empty filter panel
    And user adds a card for "DIS_POP" to the filter panel
    Then the "type of DIS_POP" reading of filter panel should be "categorical"
    And the "categories of DIS_POP" reading of filter panel should be "AS, Indigestion, PsA, Psoriasis, RA, UC"
    And all rows should pass the filter
    And counter of filter panel should be hidden

  Scenario: A name click keeps one category, a checkbox click adds another
    When user clicks on the "category RA of DIS_POP" area of filter panel
    Then 434 rows should pass the filter
    And the filter should pass exactly the rows where "DIS_POP" is "RA"
    And the "selected categories of DIS_POP" reading of filter panel should be "RA"
    And counter of filter panel should have text "1"
    When user clicks on the "checkbox UC of DIS_POP" area of filter panel
    Then 557 rows should pass the filter
    And the "selected categories of DIS_POP" reading of filter panel should be "RA, UC"
    And counter of filter panel should have text "1"
    And no errors should have been logged

  Scenario: The indicator menu selects, deselects and inverts every category
    When user picks "Deselect all" from the indicator menu of the "DIS_POP" filter card
    Then 0 rows should pass the filter
    And the "selected categories of DIS_POP" reading of filter panel should be ""
    And counter of filter panel should have text "1"
    When user picks "Invert all" from the indicator menu of the "DIS_POP" filter card
    Then all rows should pass the filter
    And the "selected categories of DIS_POP" reading of filter panel should be "AS, Indigestion, PsA, Psoriasis, RA, UC"
    And counter of filter panel should be hidden
    When user picks "Deselect all" from the indicator menu of the "DIS_POP" filter card
    Then 0 rows should pass the filter
    When user picks "Select all" from the indicator menu of the "DIS_POP" filter card
    Then all rows should pass the filter
    And counter of filter panel should be hidden
    And no errors should have been logged

  Scenario: Radio mode keeps exactly one category and offers no batch operations
    When user clicks on the "category RA of DIS_POP" area of filter panel
    Then 434 rows should pass the filter
    When user picks "Mode | Radio" from the indicator menu of the "DIS_POP" filter card
    And user closes the context menu
    Then 434 rows should pass the filter
    And the "selected categories of DIS_POP" reading of filter panel should be "RA"
    When user opens the indicator menu of the "DIS_POP" filter card
    Then context menu should contain the text "Mode"
    And context menu should not contain the text "Select all"
    And context menu should not contain the text "Invert all"
    When user closes the context menu
    And user clicks on the "category UC of DIS_POP" area of filter panel
    Then 123 rows should pass the filter
    And the "selected categories of DIS_POP" reading of filter panel should be "UC"
    When user picks "Mode | Multi-Select" from the indicator menu of the "DIS_POP" filter card
    And user closes the context menu
    Then 123 rows should pass the filter
    And no errors should have been logged

  Scenario: The in-card search narrows the table to the categories that match
    When user picks "Select all" from the indicator menu of the "DIS_POP" filter card
    Then all rows should pass the filter
    When user clicks on search icon of "DIS_POP" filter card
    And user types "Ps" into the search of the "DIS_POP" filter card
    Then 242 rows should pass the filter
    And the filter should pass exactly the rows where "DIS_POP" contains "Ps"
    And counter of filter panel should have text "1"
    When user clears the search of the "DIS_POP" filter card
    Then all rows should pass the filter
    And counter of filter panel should be hidden
    And no errors should have been logged

  Scenario: A numeric card switches to categorical and back without filtering
    When user adds a card for "AGE" to the filter panel
    Then the "type of AGE" reading of filter panel should be "histogram"
    When user hovers over "AGE" filter card
    And user clicks on "Switch to categorical filter" icon in "AGE" filter card
    Then the "type of AGE" reading of filter panel should be "categorical"
    And all rows should pass the filter
    And counter of filter panel should be hidden
    When user hovers over "AGE" filter card
    And user clicks on "Switch to histogram filter" icon in "AGE" filter card
    Then the "type of AGE" reading of filter panel should be "histogram"
    And all rows should pass the filter
    And counter of filter panel should be hidden
    And no errors should have been logged
