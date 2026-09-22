@journey @viewers @realizes:viewers.filters
Feature: Filter panel entry points
  How cards get into the panel, move in it and get out of it: the header's column picker, the
  panel's own Add Filter menu, a pinned grid header dropped on the panel (GROK-19516) and the Add
  filter link of a column's context panel all insert at the top, and a card from that link is
  suspended by its own checkbox like any other (GROK-18765) — every filter change makes the
  filtered rows the current object, so the column is clicked again before its Filter section is
  collapsed, which must happen: a section left expanded keeps a filter of its own on the column
  and the last scenario's "all rows pass" would then fail; a card dragged by its caption moves
  above another and nothing filters differently; a card's close icon removes it and
  releases whatever it was keeping, "Remove others" drops the cards that restrict nothing while
  "Remove All" empties the panel; a panel that works its cards out again leaves out a column hidden
  in the grid and gives it back once it is shown; Select Columns puts a card on every column or on
  none; a multi-value card built from its dialog lists each value its cells hold and keeps the rows
  that hold the one clicked; Filter to Column writes the panel's filtering into a boolean column;
  closing the panel releases its filtering while reopening it brings every card back as it was —
  criterion, switched-off card and removed card alike. One journey on demog-1000 (RACE: Caucasian
  896, Black 27; DIS_POP RA 434; SEX F 553). TOKENS is a column the feature makes,
  "<SEX>;<DIS_POP>", to have one whose cells hold several values.
  Not translated: the pin is proven by the grid's column order only (the grid reports no pinned
  columns); the unpinned header drag is in filter-panel-ladder.feature; the panel is first opened
  empty through its API; two values ticked together on the multi-value card and its AND/OR switch over
  them — a click on a value's name replaces the value kept, and the card reports no checkbox areas
  for its include column (only "category <value>"), so the second tick needs that area in the core;
  the GROK-12955 freeze is claimed only as "no errors" after the drags, the page staying responsive
  being what every later step already needs.

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

  Scenario: A pinned grid header dropped on the panel still makes a card
    When user picks "Pin | Pin Column" from the context menu of the "header AGE" area of grid
    Then the "column order" reading of grid should be "AGE, USUBJID, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"
    When user adds a card for "SEX" to the filter panel
    And user drags the "header AGE" area of grid onto the "view" area of filter panel
    Then "AGE" filter card should be visible
    And the "cards" reading of filter panel should be "AGE, SEX"
    And the "type of AGE" reading of filter panel should be "histogram"
    And all rows should pass the filter
    And counter of filter panel should be hidden
    And no errors should have been logged

  Scenario: A card from the column's Add filter link goes on top and its checkbox suspends it
    Given the context panel is open
    When user clicks on the "header DIS_POP" area of grid
    Then the context panel should show "DIS_POP"
    When user expands "Filter" section in context panel
    And user clicks on "Add filter" text in context panel
    Then "DIS_POP" filter card should be visible
    And the "cards" reading of filter panel should be "DIS_POP, AGE, SEX"
    And all rows should pass the filter
    When user adds a range filter on "AGE" from 30 to 60
    Then 708 rows should pass the filter
    When user clicks on the "category RA of DIS_POP" area of filter panel
    Then 299 rows should pass the filter
    And the "filters" reading of filter panel should be 2
    And counter of filter panel should have text "2"
    When user hovers over "DIS_POP" filter card
    And user unchecks checkbox of "DIS_POP" filter card
    Then 708 rows should pass the filter
    And "DIS_POP" filter card should be disabled
    And the "filtering of DIS_POP" reading of filter panel should be "false"
    And the "filters" reading of filter panel should be 1
    And counter of filter panel should have text "1"
    When user hovers over "DIS_POP" filter card
    And user checks checkbox of "DIS_POP" filter card
    Then 299 rows should pass the filter
    And the "selected categories of DIS_POP" reading of filter panel should be "RA"
    And the "filters" reading of filter panel should be 2
    And counter of filter panel should have text "2"
    When user clicks on the "header DIS_POP" area of grid
    Then the context panel should show "DIS_POP"
    When user collapses "Filter" section in context panel
    Then no errors should have been logged

  Scenario: An expression card from the panel's menu goes on top and filters nothing yet
    When user picks "Add Filter | Expression" from the filter panel menu
    Then "Expression" filter card should be visible
    And the "cards" reading of filter panel should be "Expression, DIS_POP, AGE, SEX"
    And the "type of Expression" reading of filter panel should be "expression"
    And the "filtering of Expression" reading of filter panel should be "false"
    And 299 rows should pass the filter
    And counter of filter panel should have text "2"
    And no errors should have been logged

  Scenario: A card dragged by its caption moves and nothing filters differently
    When user drags the "AGE" filter card above the "Expression" filter card
    Then the "cards" reading of filter panel should be "AGE, Expression, DIS_POP, SEX"
    And 299 rows should pass the filter
    And counter of filter panel should have text "2"
    When user drags the "Expression" filter card above the "AGE" filter card
    And user drags the "DIS_POP" filter card above the "AGE" filter card
    Then the "cards" reading of filter panel should be "Expression, DIS_POP, AGE, SEX"
    And 299 rows should pass the filter
    When user picks "Remove All" from the filter panel menu
    Then the filter panel should have 0 filters
    And all rows should pass the filter
    And no errors should have been logged

  Scenario: Filter to Column writes the panel's filtering into a boolean column
    When user adds a card for "SEX" to the filter panel
    And user clicks on the "category F of SEX" area of filter panel
    Then 553 rows should pass the filter
    When user picks "Filter to Column..." from the filter panel menu
    Then Name input in "Filter to Column" dialog should not have the value ""
    When user clicks on OK button in "Filter to Column" dialog
    Then the table should have 12 columns
    And "SEX: F" column should have type "bool"
    And the filter should pass exactly the rows where "SEX: F" is "true"
    And 553 rows should pass the filter
    When user picks "Remove All" from the filter panel menu
    Then all rows should pass the filter
    And the filter panel should have 0 filters
    When user removes "SEX: F" column
    Then the table should have 11 columns
    And all rows should pass the filter
    And no errors should have been logged

  Scenario: A column hidden in the grid gets no card when the panel works its cards out again
    When user picks "Hide" from the context menu of the "header RACE" area of grid
    Then the "column order" reading of grid should not contain "RACE"
    When user clicks on close icon of filters viewer
    Then filter panel should be hidden
    When user clicks on filter icon in toolbar
    Then filter panel should be visible
    And the "cards" reading of filter panel should not contain "RACE"
    And the "cards" reading of filter panel should contain "DIS_POP"
    And all rows should pass the filter
    When user picks "Order or Hide Columns..." from the context menu of the "header SEX" area of grid
    Then the "text of cell 4 of __name" reading of grid in "Order or Hide Columns" dialog should be "RACE"
    When user clicks on the "cell 4 of x" area of grid in "Order or Hide Columns" dialog
    Then the "column order" reading of grid should contain "RACE"
    When user clicks on CLOSE button in "Order or Hide Columns" dialog
    And user picks "Remove All" from the viewer menu of filter panel
    And user clicks on close icon of filters viewer
    And user clicks on filter icon in toolbar
    Then the "cards" reading of filter panel should contain "RACE"
    And all rows should pass the filter
    When user picks "Remove All" from the viewer menu of filter panel
    Then the filter panel should have 0 filters
    And no errors should have been logged

  Scenario: Select Columns puts a card on every column, and on none
    When user picks "Select Columns..." from the filter panel menu
    Then "0 checked" text in "Select columns..." dialog should be visible
    When user clicks on "All" link in "Select columns..." dialog
    Then "11 checked" text in "Select columns..." dialog should be visible
    When user clicks on OK button in "Select columns..." dialog
    Then the filter panel should have 11 filters
    And the "cards" reading of filter panel should be "SEVERITY, STARTED, CONTROL, DEMOG, WEIGHT, HEIGHT, DIS_POP, RACE, SEX, AGE, USUBJID"
    And all rows should pass the filter
    When user picks "Select Columns..." from the viewer menu of filter panel
    Then "11 checked" text in "Select columns..." dialog should be visible
    When user clicks on "None" link in "Select columns..." dialog
    Then "0 checked" text in "Select columns..." dialog should be visible
    When user clicks on OK button in "Select columns..." dialog
    Then the filter panel should have 0 filters
    And all rows should pass the filter
    And no errors should have been logged

  Scenario: A multi-value card built from its dialog filters on each value its cells hold
    When user adds a calculated column "TOKENS" with formula "${SEX} + \";\" + ${DIS_POP}"
    And user picks "Add Filter | Multi Value..." from the filter panel menu
    Then OK button in "Select column and separator for multi-value filter" dialog should be disabled
    And Separator input in "Select column and separator for multi-value filter" dialog should be invalid
    When user picks "TOKENS" in the column selector of the "Select column and separator for multi-value filter" dialog
    And user types ";" into Separator input in "Select column and separator for multi-value filter" dialog
    Then OK button in "Select column and separator for multi-value filter" dialog should be enabled
    When user clicks on OK button in "Select column and separator for multi-value filter" dialog
    Then "TOKENS" filter card should be visible
    And the filter panel should have 1 filter
    And the "type of TOKENS" reading of filter panel should be "multi-value"
    And the "categories of TOKENS" reading of filter panel should be "AS, F, Indigestion, M, PsA, Psoriasis, RA, UC"
    And all rows should pass the filter
    And counter of filter panel should be hidden
    When user clicks on the "category F of TOKENS" area of filter panel
    Then 553 rows should pass the filter
    And the filter should pass exactly the rows where "TOKENS" contains "F"
    And the "summary of TOKENS" reading of filter panel should be "Includes F"
    And counter of filter panel should have text "1"
    When user clicks on the "category RA of TOKENS" area of filter panel
    Then 434 rows should pass the filter
    And the filter should pass exactly the rows where "TOKENS" contains "RA"
    And the "summary of TOKENS" reading of filter panel should be "Includes RA"
    When user picks "Remove All" from the filter panel menu
    And user removes "TOKENS" column
    Then the table should have 11 columns
    And no errors should have been logged

  Scenario: Closing the panel releases its filtering and reopening restores every card as it was
    When user adds a card for "RACE" to the filter panel
    And user adds a card for "AGE" to the filter panel
    And user adds a card for "SEX" to the filter panel
    And user clicks on the "category Black of RACE" area of filter panel
    Then 27 rows should pass the filter
    When user hovers over "AGE" filter card
    And user unchecks checkbox of "AGE" filter card
    And user hovers over "SEX" filter card
    And user clicks on close of "SEX" filter card
    Then "SEX" filter card should be absent
    And 27 rows should pass the filter
    And counter of filter panel should have text "1"
    When user clicks on close icon of filters viewer
    Then filter panel should be hidden
    And all rows should pass the filter
    When user clicks on filter icon in toolbar
    Then filter panel should be visible
    And "RACE" filter card should be visible
    And the "selected categories of RACE" reading of filter panel should be "Black"
    And 27 rows should pass the filter
    And "AGE" filter card should be disabled
    And the "enabled of AGE" reading of filter panel should be "false"
    And "SEX" filter card should be absent
    And the "cards" reading of filter panel should be "AGE, RACE"
    When user picks "Remove All" from the filter panel menu
    Then all rows should pass the filter
    And no errors should have been logged
