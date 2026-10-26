@journey @viewers @realizes:viewers.filters
Feature: Filter panel regressions
  Closed filter panel issues, each a claim that fails if its bug comes back: Reset filters puts a
  numeric card's missing-values choice back to "Keep missing values" (public#790); a min/max range
  that holds none of the column's values keeps no row (public#2364); the range handles dragged back
  past the ends of the slider keep every row (public#1343); with one card's search open, a second
  card's search icon takes the typing (public#2359); a grid header dragged onto the panel leaves the
  grid's column order alone (public#1986); Ctrl+C in a card's search copies the text selected in
  it, not a picture of the panel (public#3027); a string column holding one value is listed as a
  compact row after the cards rather than as a card of its own (core#479); adding viewers does not
  shrink the panel (public#2501);
  the column's header popup still filters after the column's card left the panel (public#2943,
  public#1242). One journey on demog-1000: HEIGHT is blank in 128 rows, AGE holds whole years from
  18 to 89, DIS_POP is RA in 434 rows and Psoriasis or PsA in 242; on the AGE card at its default bins the min handle on bin 10
  leaves 306 rows and the max handle on bin 15 then 260 (the column popup's histogram, binned
  otherwise, leaves 326 from bin 10). Every scenario but the last two empties the panel again; the
  calculated column of the single-value scenario and the viewers the last ones add stay.
  Not translated: the range of public#2364 is 88.2 to 88.8, inside the column's span, since the min
  and max fields clamp a typed end to it (filter-panel-ladder.feature: 999 becomes 89) and 95 to 120
  becomes 89 to 89, which keeps the one row aged 89; the categorical half of public#2943 — the column
  popup reports no `category <label>` areas (grid/MISSING.md); public#344 (an invalid viewer Filter
  formula) logs an error on the current master and is in MISSING-regressions.md; the close and menu icons of
  the single-value row (core#479) are not claimed — the row reports no parts.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user opens an empty filter panel
    Then the filter panel should have 0 filters
    And all rows should pass the filter

  @realizes:public#790
  Scenario: Reset filters gives back the rows a missing-values choice dropped (public#790)
    When user adds a card for "HEIGHT" to the filter panel
    And user picks "Missing values | Filter out missing values" from the indicator menu of the "HEIGHT" filter card
    And user closes the context menu
    Then 872 rows should pass the filter
    And the "filtering of HEIGHT" reading of filter panel should be "true"
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    And the "filtering of HEIGHT" reading of filter panel should be "false"
    When user opens the indicator menu of the "HEIGHT" filter card
    Then the open menu should list "Missing values > Keep missing values"
    And "Keep missing values" menu item should be selected
    And "Filter out missing values" menu item should not be selected
    When user closes the context menu
    And user picks "Remove All" from the filter panel menu
    Then the filter panel should have 0 filters
    And no errors should have been logged

  @realizes:public#2364
  Scenario: A min/max range wholly outside the column's values keeps no row (public#2364)
    When user adds a card for "AGE" to the filter panel
    And user picks "Min / max" from the indicator menu of the "AGE" filter card
    And user closes the context menu
    And user enters "88.2" into the min field of the "AGE" filter card
    And user enters "88.8" into the max field of the "AGE" filter card
    Then 0 rows should pass the filter
    And the "filtering of AGE" reading of filter panel should be "true"
    And counter of filter panel should have text "1"
    When user picks "Remove All" from the filter panel menu
    Then all rows should pass the filter
    And no errors should have been logged

  @realizes:public#1343
  Scenario: Range handles brought back to the ends keep every row (public#1343)
    When user adds a card for "AGE" to the filter panel
    And user hovers over the "view" area of histogram viewer in filter panel
    And user drags the "range min handle" area of histogram viewer in filter panel to the "bin 10" area
    Then 306 rows should pass the filter
    When user hovers over the "view" area of histogram viewer in filter panel
    And user drags the "range max handle" area of histogram viewer in filter panel to the "bin 15" area
    Then 260 rows should pass the filter
    When user hovers over the "view" area of histogram viewer in filter panel
    And user drags the "range min handle" area of histogram viewer in filter panel by 400 pixels to the left
    And user hovers over the "view" area of histogram viewer in filter panel
    And user drags the "range max handle" area of histogram viewer in filter panel by 400 pixels to the right
    Then all rows should pass the filter
    And the "min of AGE" reading of filter panel should be 18
    And the "max of AGE" reading of filter panel should be 89
    When user picks "Remove All" from the filter panel menu
    Then the filter panel should have 0 filters
    And no errors should have been logged

  @realizes:public#2359
  Scenario: A second card's search icon takes the typing while the first search is open (public#2359)
    When user adds a card for "RACE" to the filter panel
    And user adds a card for "DIS_POP" to the filter panel
    And user hovers over "RACE" filter card
    And user clicks on search icon of "RACE" filter card
    And user hovers over "DIS_POP" filter card
    And user clicks on search icon of "DIS_POP" filter card
    And user types "RA" at the caret
    Then 434 rows should pass the filter
    And the filter should pass exactly the rows where "DIS_POP" is "RA"
    And the "filtering of DIS_POP" reading of filter panel should be "true"
    And the "filtering of RACE" reading of filter panel should be "false"
    When user picks "Remove All" from the filter panel menu
    Then all rows should pass the filter
    And no errors should have been logged

  @realizes:public#1986
  Scenario: A header dragged onto the panel stays where it was in the grid (public#1986)
    Then the "column order" reading of grid should be "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"
    When user drags the "header RACE" area of grid onto the "view" area of filter panel
    Then "RACE" filter card should be visible
    And the "cards" reading of filter panel should be "RACE"
    And the "column order" reading of grid should be "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"
    When user picks "Remove All" from the filter panel menu
    Then the filter panel should have 0 filters
    And no errors should have been logged

  @realizes:public#3027
  Scenario: Text selected in a card's search is copied as text (public#3027)
    When user adds a card for "DIS_POP" to the filter panel
    And user hovers over "DIS_POP" filter card
    And user clicks on search icon of "DIS_POP" filter card
    And user types "Ps" into the search of the "DIS_POP" filter card
    Then 242 rows should pass the filter
    When user presses Control+A in search of "DIS_POP" filter card
    And user presses Control+C
    Then the clipboard should have text "Ps"
    When user picks "Remove All" from the filter panel menu
    Then all rows should pass the filter
    And no errors should have been logged

  @realizes:core#479
  Scenario: A column holding one value is listed as a row after the cards (core#479)
    When user adds a calculated column "probe_class" with formula '"low"'
    And user adds a card for "RACE" to the filter panel
    And user drags the "header probe_class" area of grid onto the "view" area of filter panel
    Then the filter panel should have a filter on "probe_class" column
    And filter panel should have a "card probe_class" area
    And "probe_class" filter card should be hidden
    And the "card probe_class" area of filter panel should lie below the "card RACE" area
    And all rows should pass the filter
    When user picks "Remove All" from the filter panel menu
    Then the filter panel should have 0 filters
    And no errors should have been logged

  @realizes:public#2501
  Scenario: Adding viewers does not shrink the filter panel (public#2501)
    When user adds a card for "RACE" to the filter panel
    And user remembers the size of filter panel
    And user adds histogram viewer
    And user adds bar chart viewer
    And user adds line chart viewer
    Then the open tableview should have 1 line chart viewer
    And filter panel should keep its remembered width
    And no errors should have been logged

  @realizes:public#2943 @realizes:public#1242
  Scenario: The header popup still filters after the column's card left the panel (public#2943, public#1242)
    When user adds a card for "AGE" to the filter panel
    And user hovers over "AGE" filter card
    And user clicks on close of "AGE" filter card
    Then "AGE" filter card should be absent
    And all rows should pass the filter
    When user hovers over the "header AGE" area of grid
    And user clicks on "Column options" icon in grid
    Then column popup should be visible
    When user hovers over the "view" area of histogram viewer in column popup
    And user drags the "range min handle" area of histogram viewer in column popup to the "bin 10" area
    Then 326 rows should pass the filter
    When user presses Escape
    And user resets the filter
    Then all rows should pass the filter
    And no errors should have been logged
