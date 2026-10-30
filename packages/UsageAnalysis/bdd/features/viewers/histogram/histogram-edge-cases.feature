@viewers @realizes:viewers.histogram
Feature: Histogram edge cases
  Guards for the histogram at the edges of its settings and its data, on demog-1000 with a histogram
  of AGE (18 to 89). Bins out of its bounds is a row of `viewer-rules-edge-cases.feature` (core#1200). A range
  normalized to the filter with every row filtered out logs no error: degenerate ranges made the
  histogram's axis throw on every frame until the histogram audit core#998 fixed them.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a histogram viewer with:
      | Value             | AGE  |
      | Show Range Inputs | true |
      | Filtering Enabled | true |
    And user resizes histogram viewer to 500 by 400
    Then histogram viewer should show 1000 rows
    And the "error" reading of histogram viewer should be ""

  @edge-cases @realizes:core#998
  Scenario: A range normalized to the filter logs no error when every row is filtered out
    When user sets "Normalize To Filter" property of histogram viewer to "true"
    And user enters "60" into the "range min input" area of histogram viewer
    Then the filter should pass exactly the rows where "AGE" is between 60 and 89
    When user adds a categorical filter on "SEX" keeping ""
    Then 0 rows should pass the filter
    And histogram viewer should show 0 rows
    And no errors should have been logged

  @edge-cases
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

  @edge-cases
  Scenario: On a date column a range normalized to the filter logs no error when every row is filtered out
    When user sets "Value" property of histogram viewer to "STARTED"
    And user sets "Normalize To Filter" property of histogram viewer to "true"
    And user enters "1/1/1991" into the "range min input" area of histogram viewer
    Then 483 rows should pass the filter
    When user adds a categorical filter on "SEX" keeping ""
    Then 0 rows should pass the filter
    And histogram viewer should show 0 rows
    And no errors should have been logged

  @edge-cases
  Scenario: Text typed into the range min input is not taken and the filter stays as it was
    When user enters "abc" into the "range min input" area of histogram viewer
    Then all rows should pass the filter
    And histogram viewer should show 1000 rows
    And no errors should have been logged

  @edge-cases
  Scenario: Switching Value while a range filters releases the range
    When user enters "30" into the "range min input" area of histogram viewer
    And user enters "50" into the "range max input" area of histogram viewer
    Then the filter should pass exactly the rows where "AGE" is between 30 and 50
    When user sets "Value" property of histogram viewer to "WEIGHT"
    Then all rows should pass the filter
    And no errors should have been logged

  @edge-cases
  Scenario: A layout whose Value column the table no longer has leaves the histogram empty with a message
    When user sets "Value" property of histogram viewer to "HEIGHT"
    And user saves the layout of the current table view
    And user removes "HEIGHT" column
    And user loads the saved layout
    Then histogram viewer should contain the text "does not exist"
    And histogram viewer should not have a "bin 1" area
    And no errors should have been logged
