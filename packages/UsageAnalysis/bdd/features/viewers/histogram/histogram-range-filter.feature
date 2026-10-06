@journey @viewers @realizes:viewers.histogram
Feature: Histogram range filter
  The histogram filters the table through the range under its plot: the two range inputs write the
  slider's values, the slider clamps a bound outside the column's extent, an inverted range
  collapses the filter and says so, a handle dragged with the pointer moves the same range, and a
  narrowed range zooms the horizontal axis, rescales the vertical one under Normalize To Filter and
  survives a stacked split; a saved layout carries the look but not the range. One journey on
  demog-1000 with a histogram of AGE; every scenario puts the full range back.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a histogram viewer with:
      | Value             | AGE  |
      | Show Range Inputs | true |
      | Filtering Enabled | true |
    And user resizes histogram viewer to 500 by 400
    Then all rows should pass the filter
    And histogram viewer should show 1000 rows
    And the "range min" reading of histogram viewer should be 18
    And the "range max" reading of histogram viewer should be 89
    And the "bins shown" reading of histogram viewer should be 20
    And histogram viewer should have a "range min input" area
    And histogram viewer should have a "range max input" area
    And histogram viewer should have a "bin 8" area

  Scenario: Typed bounds filter the table and zoom the axis to them
    When user enters "30" into the "range min input" area of histogram viewer
    Then 861 rows should pass the filter
    And the "range min" reading of histogram viewer should be 30
    And the "axis min" reading of histogram viewer should be higher than before
    And the "bins shown" reading of histogram viewer should be lower than before
    When user enters "60" into the "range max input" area of histogram viewer
    Then 708 rows should pass the filter
    And the filter should pass exactly the rows where "AGE" is between 30 and 60
    And histogram viewer should show 708 rows
    And the "range max" reading of histogram viewer should be 60
    And the "bins shown" reading of histogram viewer should be lower than before
    And histogram viewer should have a "range bar" area
    When user enters "18" into the "range min input" area of histogram viewer
    And user enters "89" into the "range max input" area of histogram viewer
    Then all rows should pass the filter
    And the "bins shown" reading of histogram viewer should be 20
    And the "range min" reading of histogram viewer should be 18
    And no errors should have been logged

  Scenario: A stacked split keeps the range filter
    When user enters "30" into the "range min input" area of histogram viewer
    And user enters "60" into the "range max input" area of histogram viewer
    Then 708 rows should pass the filter
    When user sets properties of histogram viewer:
      | Split       | SEX  |
      | Split Stack | true |
    Then 708 rows should pass the filter
    And histogram viewer should show 708 rows
    And histogram viewer should have a "bin 8 | F" area
    And histogram viewer should have a "bin 8 | M" area
    And the "bin 8 | F" and "bin 8 | M" areas of histogram viewer should be painted in different colors
    When user sets properties of histogram viewer:
      | Split       |       |
      | Split Stack | false |
    And user enters "18" into the "range min input" area of histogram viewer
    And user enters "89" into the "range max input" area of histogram viewer
    Then all rows should pass the filter
    And histogram viewer should have a "bin 8" area
    And no errors should have been logged

  Scenario: Zoom To Range decides whether the axis follows the range
    When user enters "60" into the "range min input" area of histogram viewer
    Then 171 rows should pass the filter
    And the "axis min" reading of histogram viewer should be higher than before
    And the "bins shown" reading of histogram viewer should be lower than before
    When user sets "Zoom To Range" property of histogram viewer to "false"
    Then the "axis min" reading of histogram viewer should be 18
    And the "bins shown" reading of histogram viewer should be 20
    And histogram viewer should have repainted by at least 500 pixels
    And 171 rows should pass the filter
    When user sets "Zoom To Range" property of histogram viewer to "true"
    Then the "axis min" reading of histogram viewer should be higher than before
    And the "bins shown" reading of histogram viewer should be lower than before
    When user enters "18" into the "range min input" area of histogram viewer
    Then all rows should pass the filter
    And the "bins shown" reading of histogram viewer should be 20
    And no errors should have been logged

  Scenario: Normalize To Filter scales the bars to the bins the range leaves
    When user enters "60" into the "range min input" area of histogram viewer
    Then 171 rows should pass the filter
    And the "y axis max" reading of histogram viewer should be 68
    When user sets "Normalize To Filter" property of histogram viewer to "false"
    Then the "y axis max" reading of histogram viewer should be 99
    And histogram viewer should have repainted by at least 500 pixels
    When user sets "Normalize To Filter" property of histogram viewer to "true"
    Then the "y axis max" reading of histogram viewer should be 68
    When user enters "18" into the "range min input" area of histogram viewer
    Then all rows should pass the filter
    And the "y axis max" reading of histogram viewer should be 99
    And no errors should have been logged

  Scenario: A maximum below the minimum collapses the filter and reports it
    When user enters "40" into the "range min input" area of histogram viewer
    Then 665 rows should pass the filter
    And the "error" reading of histogram viewer should be ""
    When user enters "20" into the "range max input" area of histogram viewer
    Then 0 rows should pass the filter
    And the "error" reading of histogram viewer should be "max should be greater than min"
    And histogram viewer should not have a "bin 8" area
    When user enters "60" into the "range max input" area of histogram viewer
    Then 512 rows should pass the filter
    And the "error" reading of histogram viewer should be ""
    And histogram viewer should have a "bin 8" area
    When user enters "18" into the "range min input" area of histogram viewer
    And user enters "89" into the "range max input" area of histogram viewer
    Then all rows should pass the filter
    And no errors should have been logged

  Scenario: Bounds outside the column extent are clamped, the typed text is kept
    When user enters "40" into the "range min input" area of histogram viewer
    Then fewer than 1000 rows should pass the filter
    And the "range min" reading of histogram viewer should be 40
    When user enters "-999" into the "range min input" area of histogram viewer
    Then all rows should pass the filter
    And the "range min" reading of histogram viewer should be 18
    And the range min input of histogram viewer should read "-999"
    When user enters "60" into the "range max input" area of histogram viewer
    Then 847 rows should pass the filter
    And the "range min" reading of histogram viewer should be 18
    When user enters "999" into the "range max input" area of histogram viewer
    Then all rows should pass the filter
    And the "range max" reading of histogram viewer should be 89
    And the range max input of histogram viewer should read "999"
    When user enters "18" into the "range min input" area of histogram viewer
    And user enters "89" into the "range max input" area of histogram viewer
    Then all rows should pass the filter
    And the range min input of histogram viewer should read "18"
    And no errors should have been logged

  Scenario: A handle dragged with the pointer moves the range, a double click resets it
    When user hovers over the "view" area of histogram viewer
    Then histogram viewer should have a "range min handle" area
    And histogram viewer should have a "range max handle" area
    When user drags the "range max handle" area of histogram viewer to the "bin 12" area
    Then fewer than 1000 rows should pass the filter
    And the "range max" reading of histogram viewer should be lower than before
    And histogram viewer should show fewer rows than before
    And histogram viewer should have a "range bar" area
    When user double-clicks on the "range slider" area of histogram viewer
    Then all rows should pass the filter
    And the "range max" reading of histogram viewer should be 89
    And the "bins shown" reading of histogram viewer should be 20
    And no errors should have been logged

  Scenario: Filtering Enabled off releases the table and keeps the zoom
    When user enters "40" into the "range min input" area of histogram viewer
    Then 665 rows should pass the filter
    And the "bins shown" reading of histogram viewer should be 14
    When user sets "Filtering Enabled" property of histogram viewer to "false"
    Then all rows should pass the filter
    And the "range min" reading of histogram viewer should be 40
    And the "bins shown" reading of histogram viewer should be 14
    When user sets "Filtering Enabled" property of histogram viewer to "true"
    Then 665 rows should pass the filter
    When user enters "18" into the "range min input" area of histogram viewer
    Then all rows should pass the filter
    And no errors should have been logged

  Scenario: A saved layout brings the histogram back at the full range
    When user enters "40" into the "range min input" area of histogram viewer
    And user enters "60" into the "range max input" area of histogram viewer
    Then 512 rows should pass the filter
    When user saves the layout of the current table view to the server
    And user clicks on close icon of histogram viewer
    Then histogram viewer should be absent
    And all rows should pass the filter
    When user loads the saved layout
    Then histogram viewer should be visible
    And "Value" property of histogram viewer should be "AGE"
    And histogram viewer should have a "range min input" area
    And the "range min" reading of histogram viewer should be 18
    And the "range max" reading of histogram viewer should be 89
    And the "bins shown" reading of histogram viewer should be 20
    And all rows should pass the filter
    And no errors should have been logged
