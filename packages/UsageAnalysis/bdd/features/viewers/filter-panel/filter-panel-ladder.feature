@journey @viewers @realizes:viewers.filters
Feature: Filter panel core ladder
  The ladder every other filter panel feature stands on: a grid column header dragged onto the panel
  becomes a card that filters nothing, a click on a category name keeps that category alone while a
  click on its checkbox adds one, a second criterion composes with the first, the min and max fields
  of a histogram card still answer after an end out of the column's range (github-2307), the header
  counter counts the cards that restrict rows and its tooltip names them, the master toggle stashes
  and gives back every card's state and keeps a card switched off on its own switched off, Escape on
  the focused panel does the same and keeps every criterion, the header search hides cards and
  nothing else, the reset icon clears the criteria and switches every card back on (github-1103),
  and a card's own checkbox suspends its criterion without losing it. One journey on demog-1000:
  RACE is Black 27, Other 62, Caucasian 896, Asian 15, and 708 of the 1000 rows are aged 30 to 60 —
  19 of them Black (18 up to 55, 25 from 30 on).
  Not translated: the tooltip's freshness against the one the previous hover left — the library's
  tooltip is the one shown now; the layout and project round-trips of the md are in
  filter-panel-persistence.feature; the panel is emptied through its API rather than Remove All,
  and AGE's window is set through the card's state before its min and max fields are typed into.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user opens an empty filter panel
    Then the filter panel should have 0 filters
    And all rows should pass the filter
    And counter of filter panel should be hidden

  Scenario: A card dragged from a grid header filters nothing
    When user drags the "header RACE" area of grid onto the "view" area of filter panel
    Then "RACE" filter card should be visible
    And the "cards" reading of filter panel should be "RACE"
    And the "type of RACE" reading of filter panel should be "categorical"
    And the "categories of RACE" reading of filter panel should be "Asian, Black, Caucasian, Other"
    And all rows should pass the filter
    And the "filters" reading of filter panel should be 0
    And the "filtering of RACE" reading of filter panel should be "false"
    And counter of filter panel should be hidden
    And no errors should have been logged

  Scenario: A click on a category name keeps that category alone
    When user clicks on the "category Black of RACE" area of filter panel
    Then 27 rows should pass the filter
    And the filter should pass exactly the rows where "RACE" is "Black"
    And the "selected categories of RACE" reading of filter panel should be "Black"
    And the "rows shown" reading of filter panel should be 27
    And the "filters" reading of filter panel should be 1
    And counter of filter panel should be visible
    And counter of filter panel should have text "1"
    And no errors should have been logged

  Scenario: A click on a checkbox adds a category to the ones kept
    When user clicks on the "checkbox Other of RACE" area of filter panel
    Then 89 rows should pass the filter
    And the "selected categories of RACE" reading of filter panel should be "Black, Other"
    And counter of filter panel should have text "1"
    When user clicks on the "category Black of RACE" area of filter panel
    Then 27 rows should pass the filter
    And the "selected categories of RACE" reading of filter panel should be "Black"
    And no errors should have been logged

  Scenario: A second criterion composes with the first
    When user adds a range filter on "AGE" from 30 to 60
    Then "AGE" filter card should be visible
    And 19 rows should pass the filter
    And the "min of AGE" reading of filter panel should be 30
    And the "max of AGE" reading of filter panel should be 60
    And the "filters" reading of filter panel should be 2
    And counter of filter panel should have text "2"
    And no errors should have been logged

  Scenario: The max field still answers after an end out of the column's range
    When user picks "Min / max" from the indicator menu of the "AGE" filter card
    And user closes the context menu
    And user enters "55" into the max field of the "AGE" filter card
    Then 18 rows should pass the filter
    And the "max of AGE" reading of filter panel should be 55
    When user enters "999" into the max field of the "AGE" filter card
    Then 25 rows should pass the filter
    And the "max of AGE" reading of filter panel should be 89
    When user enters "60" into the max field of the "AGE" filter card
    Then 19 rows should pass the filter
    And the "max of AGE" reading of filter panel should be 60
    And the "min of AGE" reading of filter panel should be 30
    And counter of filter panel should have text "2"
    And no errors should have been logged

  Scenario: The counter's tooltip names the cards that restrict rows
    When user adds a card for "SEX" to the filter panel
    Then the "filters" reading of filter panel should be 2
    When user hovers over counter of filter panel
    Then tooltip should contain the text "RACE"
    And tooltip should contain the text "Black"
    And tooltip should contain the text "AGE"
    And tooltip should contain the text "[30,60]"
    And tooltip should not contain the text "Caucasian"
    And tooltip should not contain the text "Other"
    And tooltip should not contain the text "SEX"
    And no errors should have been logged

  Scenario: The master toggle stashes every card's state and gives it back
    When user hovers over filter panel
    And user unchecks master of filter panel
    Then all rows should pass the filter
    And the "active" reading of filter panel should be "false"
    And the "cards" reading of filter panel should be "SEX, AGE, RACE"
    And the "selected categories of RACE" reading of filter panel should be "Black"
    And "AGE" filter card should be disabled
    And counter of filter panel should be hidden
    When user checks master of filter panel
    Then 19 rows should pass the filter
    And counter of filter panel should be visible
    And the "active" reading of filter panel should be "true"
    And counter of filter panel should have text "2"
    And no errors should have been logged

  Scenario: The master toggle keeps a card switched off on its own switched off
    When user hovers over "RACE" filter card
    And user unchecks checkbox of "RACE" filter card
    Then 708 rows should pass the filter
    And "RACE" filter card should be disabled
    And counter of filter panel should have text "1"
    When user hovers over filter panel
    And user unchecks master of filter panel
    Then all rows should pass the filter
    And the "active" reading of filter panel should be "false"
    When user checks master of filter panel
    Then 708 rows should pass the filter
    And the "active" reading of filter panel should be "true"
    And "RACE" filter card should be disabled
    And the "enabled of RACE" reading of filter panel should be "false"
    And counter of filter panel should have text "1"
    When user hovers over "RACE" filter card
    Then checkbox of "RACE" filter card should be unchecked
    When user switches the "RACE" filter card back on
    Then 19 rows should pass the filter
    And "RACE" filter card should be enabled
    And the "selected categories of RACE" reading of filter panel should be "Black"
    And counter of filter panel should have text "2"
    And no errors should have been logged

  Scenario: Escape on the focused panel switches the cards off and gives them back unchanged
    When user clicks on empty plot space of filter panel
    And user presses Escape
    Then all rows should pass the filter
    And the "active" reading of filter panel should be "false"
    And "AGE" filter card should be disabled
    And the "selected categories of RACE" reading of filter panel should be "Black"
    And the "min of AGE" reading of filter panel should be 30
    And the "max of AGE" reading of filter panel should be 60
    And counter of filter panel should be hidden
    When user presses Escape
    Then 19 rows should pass the filter
    And the "active" reading of filter panel should be "true"
    And counter of filter panel should have text "2"
    And no errors should have been logged

  Scenario: The header search hides cards and leaves the rows alone
    When user hovers over filter panel
    And user clicks on search icon of filter panel
    And user types "RACE" into search of filter panel
    Then "AGE" filter card should be hidden
    And "SEX" filter card should be hidden
    And "RACE" filter card should be visible
    And 19 rows should pass the filter
    And counter of filter panel should have text "2"
    When user clears search of filter panel
    Then "AGE" filter card should be visible
    And "SEX" filter card should be visible
    And 19 rows should pass the filter
    And no errors should have been logged

  Scenario: The reset icon clears the criteria, keeps the cards and switches every one back on
    When user hovers over "AGE" filter card
    And user unchecks checkbox of "AGE" filter card
    Then 27 rows should pass the filter
    And "AGE" filter card should be disabled
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    And the "filters" reading of filter panel should be 0
    And counter of filter panel should be hidden
    And the "cards" reading of filter panel should be "SEX, AGE, RACE"
    And "AGE" filter card should be enabled
    And "RACE" filter card should be enabled
    And "SEX" filter card should be enabled
    And the "enabled of AGE" reading of filter panel should be "true"
    When user hovers over "AGE" filter card
    Then checkbox of "AGE" filter card should be checked
    When user hovers over "RACE" filter card
    Then checkbox of "RACE" filter card should be checked
    When user hovers over "SEX" filter card
    Then checkbox of "SEX" filter card should be checked
    And the "filtering of RACE" reading of filter panel should be "false"
    And no errors should have been logged

  Scenario: A card's checkbox suspends its criterion and keeps it
    When user clicks on the "category Black of RACE" area of filter panel
    Then 27 rows should pass the filter
    When user hovers over "RACE" filter card
    And user unchecks checkbox of "RACE" filter card
    Then all rows should pass the filter
    And "RACE" filter card should be disabled
    And the "enabled of RACE" reading of filter panel should be "false"
    And the "selected categories of RACE" reading of filter panel should be "Black"
    And counter of filter panel should be hidden
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then "RACE" filter card should be enabled
    And the "enabled of RACE" reading of filter panel should be "true"
    And all rows should pass the filter
    And no errors should have been logged
