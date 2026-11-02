@journey @viewers @realizes:viewers.histogram
Feature: Histogram range filter
  The histogram filters the table through the range under its plot: the two range inputs write the
  slider's values, the slider clamps a bound outside the column's extent, an inverted range
  collapses the filter and says so, a handle dragged with the pointer moves the same range, and a
  narrowed range zooms the horizontal axis, rescales the vertical one under Normalize To Filter and
  survives a stacked split; a saved layout carries the look but not the range. Text typed into a
  range input is not taken and the filter stays as it was. Switching Value while a range filters
  releases the range. A range normalized to the filter logs no error when every row is filtered out,
  on AGE and on a date column (degenerate ranges made the axis throw on every frame until core#998).
  On the date column STARTED a maximum below the minimum empties the filter and is reported, and a
  bound past the dates is taken (483 rows start on or after 1/1/1991). A Bins value typed in the
  context panel outside the property's declared bounds (201) is not taken; the property keeps its
  value. The panel's number input marks it invalid but writes it into the property (core#1200). One
  journey on demog-1000 with a histogram of AGE (18 to 89); every scenario puts the full range back.

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
    Then 665 rows should pass the filter
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
    When user enters "40" into the "range min input" area of histogram viewer
    Then 665 rows should pass the filter
    When user enters "abc" into the "range min input" area of histogram viewer
    Then 665 rows should pass the filter
    And the "range min" reading of histogram viewer should be 40
    And no errors should have been logged
    When user enters "18" into the "range min input" area of histogram viewer
    Then all rows should pass the filter

  @realizes:public#2301
  Scenario: A handle dragged with the pointer moves the range, a double click resets it, and the min handle dragged past the max leaves a valid range (public#2301)
    When user hovers over the "view" area of histogram viewer
    Then histogram viewer should have a "range min handle" area
    And histogram viewer should have a "range max handle" area
    When user drags the "range max handle" area of histogram viewer to the "bin 12" area
    Then 809 rows should pass the filter
    And the "range max" reading of histogram viewer should be lower than before
    And histogram viewer should show fewer rows than before
    And histogram viewer should have a "range bar" area
    When user double-clicks on the "range slider" area of histogram viewer
    Then all rows should pass the filter
    And the "range max" reading of histogram viewer should be 89
    And the "bins shown" reading of histogram viewer should be 20
    When user moves the pointer away from histogram viewer
    And user hovers over the "view" area of histogram viewer
    Then histogram viewer should have a "range min handle" area
    When user drags the "range min handle" area of histogram viewer by 450 pixels to the right
    Then the "range min" reading of histogram viewer should be 89
    And the "range max" reading of histogram viewer should be 89
    And 1 row should pass the filter
    And the "error" reading of histogram viewer should be ""
    When user double-clicks on the "range slider" area of histogram viewer
    Then all rows should pass the filter
    And the "range min" reading of histogram viewer should be 18
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

  @realizes:core#998
  Scenario: A range normalized to the filter logs no error when every row is filtered out (core#998)
    When user sets "Normalize To Filter" property of histogram viewer to "true"
    And user enters "60" into the "range min input" area of histogram viewer
    Then the filter should pass exactly the rows where "AGE" is between 60 and 89
    When user adds a categorical filter on "SEX" keeping ""
    Then 0 rows should pass the filter
    And histogram viewer should show 0 rows
    And no errors should have been logged
    When user adds a categorical filter on "SEX" keeping "F, M"
    And user enters "18" into the "range min input" area of histogram viewer
    And user sets "Normalize To Filter" property of histogram viewer to "false"
    Then all rows should pass the filter
    And no errors should have been logged

  Scenario: Switching Value while a range filters releases the range
    When user enters "30" into the "range min input" area of histogram viewer
    And user enters "50" into the "range max input" area of histogram viewer
    Then the filter should pass exactly the rows where "AGE" is between 30 and 50
    When user sets "Value" property of histogram viewer to "WEIGHT"
    Then all rows should pass the filter
    And no errors should have been logged
    When user sets "Value" property of histogram viewer to "AGE"
    Then all rows should pass the filter
    And the "range min" reading of histogram viewer should be 18

  Scenario: On a date column a maximum below the minimum empties the filter and is reported, and a bound past the dates is taken
    When user sets "Value" property of histogram viewer to "STARTED"
    And user enters "1/1/1991" into the "range min input" area of histogram viewer
    Then 483 rows should pass the filter
    And the "error" reading of histogram viewer should be ""
    When user enters "1/1/1990" into the "range max input" area of histogram viewer
    Then 0 rows should pass the filter
    And the "error" reading of histogram viewer should be "max should be greater than min"
    When user enters "1/1/1992" into the "range max input" area of histogram viewer
    Then 483 rows should pass the filter
    And the "error" reading of histogram viewer should be ""
    And the range max input of histogram viewer should read "1/1/1992"
    And no errors should have been logged
    When user sets "Value" property of histogram viewer to "AGE"
    Then all rows should pass the filter

  @realizes:core#998
  Scenario: On a date column a range normalized to the filter logs no error when every row is filtered out (core#998)
    When user sets "Value" property of histogram viewer to "STARTED"
    And user sets "Normalize To Filter" property of histogram viewer to "true"
    And user enters "1/1/1991" into the "range min input" area of histogram viewer
    Then 483 rows should pass the filter
    When user adds a categorical filter on "SEX" keeping ""
    Then 0 rows should pass the filter
    And histogram viewer should show 0 rows
    And no errors should have been logged
    When user adds a categorical filter on "SEX" keeping "F, M"
    Then 483 rows should pass the filter
    When user sets properties of histogram viewer:
      | Normalize To Filter | false |
      | Value               | AGE   |
    Then all rows should pass the filter
    And no errors should have been logged

  @known-failure:278 @realizes:core#1200
  Scenario: A Bins value typed in the context panel outside its declared bounds is not taken (core#1200)
    When user clicks on settings icon of histogram viewer
    Then context panel should be visible
    Given "X Axis" category in context panel is expanded
    When user enters "201" into "Bins" property in context panel
    Then "Bins" property of histogram viewer should be "20"
