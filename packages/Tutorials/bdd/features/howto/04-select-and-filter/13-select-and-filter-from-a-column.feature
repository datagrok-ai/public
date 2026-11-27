@guide @tutorial @viewport:1280x720
Feature: How-to · Select and filter
  For the maintainer: selecting and filtering from a column of the grid, in one clip (it took in the former clip
  "Select rows with the same value"). With a cell current, Shift+Enter selects every row holding the same value in
  that column; in the table's Search box (Ctrl+F) Shift+Enter selects the rows a condition matches where Enter would
  filter them; the Column options icon a column header shows under the pointer opens the column's popup, whose
  Filter pane holds a filter for that column alone (a histogram with a range for AGE): dragging its handle filters
  the table at once, Add filter moves the filter to the filter panel, and Escape closes the popup. The 171 rows
  selected by the search stay selected through the filter (every one of them is older than the new range's start).
  demog (5850 rows: RACE Other 354, AGE over 70 in 171) is opened off camera; nothing is made on the server.

  Scenario: Select and filter from a column
    Select every row that shares a cell's value or matches a condition, and filter a column right from its header.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    Then no rows should be selected
    And all rows should pass the filter
    # caption: Click a cell, then press Shift+Enter
    When user clicks on the "cell 4 of RACE" area of grid
    And user presses Shift+Enter
    Then 354 rows should be selected
    And only rows where "RACE" is "Other" should be selected
    When user hovers over "Selected: 354" text in status bar
    And the video says "Every row with RACE = Other is selected"
    # caption: Esc clears the selection
    When user presses Escape in grid
    Then no rows should be selected
    # caption: In the Search box (Ctrl+F), Shift+Enter selects the matches
    When user presses Control+f in grid overlay
    Then table search should be focused
    When user types "AGE > 70" into table search
    And user presses Shift+Enter in table search
    Then 171 rows should be selected
    And only rows where "AGE" is between 71 and 200 should be selected
    And all rows should pass the filter
    When the video says "Enter would filter; Shift+Enter selects"
    And user hovers over the "header AGE" area of grid
    # caption: The icon on a column header opens its options
    And user clicks on "Column options" icon in grid
    Then column popup should be visible
    And title of column popup should have text "AGE"
    When user hovers over the "bin 1" area of histogram viewer in column popup
    # caption: Drag an end of the range to filter by this column
    And user drags the "range min handle" area of histogram viewer in column popup to the "bin 10" area
    Then 1889 rows should pass the filter
    And the "range min" reading of histogram viewer in column popup should be higher than before
    And 171 rows should be selected
    # caption: Add filter moves it to the filter panel
    When user clicks on "Add filter" action in column popup
    Then "AGE" filter card should be visible
    And 1889 rows should pass the filter
    # caption: Esc closes the popup; the filter stays
    When user presses Escape
    Then column popup should be absent
    And 1889 rows should pass the filter
    And no errors should have been logged
