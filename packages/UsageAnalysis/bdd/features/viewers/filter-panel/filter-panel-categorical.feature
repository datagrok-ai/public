@journey @viewers @realizes:viewers.filters
Feature: Categorical filter card
  What the categorical card does with its categories: a click on a name keeps one while a click on
  a checkbox adds one, the indicator menu's batch operations select, deselect and invert them all,
  Radio mode keeps exactly one of the two kept and drops the batch items from the menu while
  Multi-Select brings them back, the in-card search narrows the table to the categories that match
  and a pasted list of lines keeps exactly the listed ones (GROK-20242), rows deleted and brought
  back by Ctrl+Z bring their category back to the card (GROK-19537), a numeric column's card
  switches to categorical with its menu at hand (GROK-19915) and back to a histogram that keeps
  nothing it kept as categories, and a numeric column holding one value gets a histogram card that
  takes a click without an error (GROK-19897). One journey on demog-1000, whose DIS_POP is RA 434,
  Psoriasis 204, Indigestion 152, UC 123, AS 49 and PsA 38; its AGE takes 72 values from 18 to 89.
  Not translated: the row deletion is "Edit > Remove > Selected Rows" over a selection made through
  the table's API — the md's own grid-selection gesture is not the subject; Ctrl+Z is pressed for real.
  The panel is opened, and the one-value column made and removed, through the API; the card's own
  count of matching categories under the search is not read, the rows the search keeps are.

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
    When user clicks on the "checkbox UC of DIS_POP" area of filter panel
    Then 557 rows should pass the filter
    And the "selected categories of DIS_POP" reading of filter panel should be "RA, UC"
    When user picks "Mode | Radio" from the indicator menu of the "DIS_POP" filter card
    And user closes the context menu
    Then 434 rows should pass the filter
    And the "selected categories of DIS_POP" reading of filter panel should be "RA"
    When user opens the indicator menu of the "DIS_POP" filter card
    Then the open menu should list "Mode"
    And the open menu should not list "Select all"
    And the open menu should not list "Deselect all"
    And the open menu should not list "Invert all"
    When user closes the context menu
    And user clicks on the "checkbox UC of DIS_POP" area of filter panel
    Then 123 rows should pass the filter
    And the "selected categories of DIS_POP" reading of filter panel should be "UC"
    And counter of filter panel should be visible
    And counter of filter panel should have text "1"
    When user picks "Mode | Multi-Select" from the indicator menu of the "DIS_POP" filter card
    And user closes the context menu
    Then 123 rows should pass the filter
    When user opens the indicator menu of the "DIS_POP" filter card
    Then the open menu should list "Select all"
    And the open menu should list "Deselect all"
    And the open menu should list "Invert all"
    When user closes the context menu
    Then no errors should have been logged

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
    When user types "ori" into the search of the "DIS_POP" filter card
    Then 204 rows should pass the filter
    And the filter should pass exactly the rows where "DIS_POP" contains "ori"
    When user clears the search of the "DIS_POP" filter card
    Then all rows should pass the filter
    And counter of filter panel should be hidden
    And no errors should have been logged

  Scenario: A pasted list of lines keeps exactly the categories it lists
    When user pastes "RA\nUC\n" into the search of the "DIS_POP" filter card
    Then 557 rows should pass the filter
    And no rows where "DIS_POP" is "Psoriasis" should pass the filter
    And no rows where "DIS_POP" is "PsA" should pass the filter
    And all rows where "DIS_POP" is "UC" should pass the filter
    And counter of filter panel should have text "1"
    When user clears the search of the "DIS_POP" filter card
    Then all rows should pass the filter
    When user pastes "RA\nUC" into the search of the "DIS_POP" filter card
    Then 557 rows should pass the filter
    And no rows where "DIS_POP" is "Psoriasis" should pass the filter
    When user clears the search of the "DIS_POP" filter card
    Then all rows should pass the filter
    And counter of filter panel should be hidden
    And no errors should have been logged

  Scenario: Rows deleted and brought back by Ctrl+Z bring their category back to the card
    When user selects rows where "DIS_POP" is "PsA"
    And user picks "Edit > Remove > Selected Rows" from the top menu
    Then the table should have 962 rows
    And the "categories of DIS_POP" reading of filter panel should be "AS, Indigestion, Psoriasis, RA, UC"
    And filter panel should not have a "category PsA of DIS_POP" area
    When user clicks on grid
    And user presses Control+Z
    Then the table should have 1000 rows
    And the "categories of DIS_POP" reading of filter panel should be "AS, Indigestion, PsA, Psoriasis, RA, UC"
    And filter panel should have a "category PsA of DIS_POP" area
    And all rows should pass the filter
    When user clears the row selection
    Then no errors should have been logged

  Scenario: A numeric card switches to categorical and back without filtering
    When user adds a card for "AGE" to the filter panel
    Then the "type of AGE" reading of filter panel should be "histogram"
    When user hovers over "AGE" filter card
    And user clicks on "Switch to categorical filter" icon in "AGE" filter card
    Then the "type of AGE" reading of filter panel should be "categorical"
    And all rows should pass the filter
    And counter of filter panel should be hidden
    When user hovers over caption of "AGE" filter card
    Then indicator of "AGE" filter card should be visible
    When user opens the indicator menu of the "AGE" filter card
    Then the open menu should list "Select all"
    When user closes the context menu
    And user clicks on the "category 20 of AGE" area of filter panel
    Then the filter should pass exactly the rows where "AGE" is "20"
    And the "selected categories of AGE" reading of filter panel should be "20"
    And counter of filter panel should have text "1"
    When user hovers over "AGE" filter card
    And user clicks on "Switch to histogram filter" icon in "AGE" filter card
    Then the "type of AGE" reading of filter panel should be "histogram"
    And all rows should pass the filter
    And counter of filter panel should be hidden
    And no errors should have been logged

  Scenario: A column holding one value gets a card that takes a click without an error
    When user adds a calculated column "probe_constant" with formula "1"
    And user adds a card for "probe_constant" to the filter panel
    Then "probe_constant" filter card should be visible
    And the "type of probe_constant" reading of filter panel should be "histogram"
    When user clicks on body of "probe_constant" filter card
    Then all rows should pass the filter
    And no errors should have been logged
    When user removes "probe_constant" column
    Then "probe_constant" filter card should be absent
    And all rows should pass the filter
    And no errors should have been logged
