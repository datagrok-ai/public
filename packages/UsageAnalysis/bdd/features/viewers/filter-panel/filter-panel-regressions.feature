@journey @viewers @realizes:viewers.filters
Feature: Filter panel regressions
  Closed filter panel issues, each a claim that fails if its bug comes back: a missing-values choice
  that came back with a layout can be changed back, and Reset filters puts it back to "Keep missing
  values" (public#791, public#790); a card search reset before a layout was saved comes back as no
  filtering (public#3059); the order the cards were dragged into survives the layout (public#1464);
  a view cloned after Reset filters comes up unfiltered (public#319); adding viewers does not shrink
  the panel (public#2501); the column's header popup still filters after the column's card left the
  panel (public#2943, public#1242). One journey on demog-1000 (HEIGHT blank in 128 rows; DIS_POP
  Psoriasis or PsA in 242, RA in 434; RACE Black in 27; on AGE the column popup's histogram leaves 326
  rows from bin 10); every scenario but the last two empties the panel again, and the viewers the
  public#2501 scenario adds stay.
  The round's other panel issues are in the panel's journeys: public#1986 (ladder), public#3027,
  core#479, public#1889 (categorical), public#1343 (range-handles), public#377 (persistence),
  public#1837/#1460 (cloned-view), public#344/core#1141/core#1167 (scatter-plot/scatter-plot).
  Not claimed: public#2364 (the min and max fields clamp a typed end to the column's span; the range
  itself is an ApiTests `FilterGroup` test), public#2359's focus half (MISSING.md), public#1103 (the
  per-card switches through a layout are filter-panel-cloned-view's, the master toggle the ladder's).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user opens an empty filter panel
    Then the filter panel should have 0 filters
    And all rows should pass the filter

  @realizes:public#790 @realizes:public#791
  Scenario: A missing-values choice that came back with a layout can be undone, and Reset filters puts it back to Keep (public#790, public#791)
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
    When user picks "Missing values | Filter out missing values" from the indicator menu of the "HEIGHT" filter card
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

  @realizes:public#3059
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
    When user picks "Remove All" from the filter panel menu
    Then the filter panel should have 0 filters
    And no errors should have been logged

  @realizes:public#1464
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
    When user picks "Remove All" from the filter panel menu
    Then the filter panel should have 0 filters
    And no errors should have been logged

  @realizes:public#319
  Scenario: A view cloned after Reset filters comes up unfiltered (public#319)
    When user adds a card for "RACE" to the filter panel
    And user clicks on the "category Black of RACE" area of filter panel
    Then 27 rows should pass the filter
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    When user picks "View > Layout > Clone View" from the top menu
    Then the "demog-1000 copy" view should be current
    And filter panel should be visible
    And "RACE" filter card should be visible
    And all rows should pass the filter
    And the "filtering of RACE" reading of filter panel should be "false"
    When user closes the current view
    And user switches to the "demog-1000" view
    And user picks "Remove All" from the filter panel menu
    Then the filter panel should have 0 filters
    And no errors should have been logged

  @realizes:public#2501
  Scenario: Adding viewers does not shrink the filter panel (public#2501)
    When user adds a card for "RACE" to the filter panel
    And user remembers the size of filter panel
    And user remembers the size of grid
    And user adds histogram viewer
    And user adds bar chart viewer
    And user adds line chart viewer
    Then the open tableview should have 1 line chart viewer
    And grid should be narrower than remembered
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
