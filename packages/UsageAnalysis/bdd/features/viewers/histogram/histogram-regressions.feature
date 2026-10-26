@journey @viewers @realizes:viewers.histogram
Feature: Histogram regressions from closed issues
  Closed histogram issues, each claimed where the bug showed. From the histogram audit (core#998): a
  bin click under the viewer's own Filter selects only the rows the bar counts, a split with one bin
  can be hovered, and a logarithmic X axis over a column with no positive value and a Y Min above the
  Y Max are reported in the viewer's error instead of throwing on every frame. A min handle dragged
  past the max handle leaves a valid range (public#2301), and the split lines take the split column's
  category colours (public#130). One journey on demog-1000 with a histogram of AGE, 500 by 400 so the
  range slider and the selectors are drawn; every scenario puts back what it changed.
  Not translated: the audit's distribution lines after a filter made elsewhere (no reading of a
  line's statistics, see MISSING-regressions.md), its empty-table and JS-set Bins crashes (API tests) and its
  performance half.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a histogram viewer with:
      | Value | AGE |
    And user resizes histogram viewer to 500 by 400
    Then histogram viewer should show 1000 rows
    And the "bins shown" reading of histogram viewer should be 20
    And the "error" reading of histogram viewer should be ""

  # core#998 (a)
  Scenario: Under the viewer's own Filter a bin click selects only the rows the bar counts (core#998)
    When user sets "Filter" property of histogram viewer to "${SEX} == \"F\""
    Then histogram viewer should show 553 rows
    When user clicks on the "bin 8" area of histogram viewer
    Then 53 rows should be selected
    And no rows where "SEX" is "M" should be selected
    When user clears the row selection
    And user sets "Filter" property of histogram viewer to ""
    Then histogram viewer should show 1000 rows
    And no errors should have been logged

  # core#998 (b)
  Scenario: A split with a single bin can be hovered without errors (core#998)
    When user sets "Split" property of histogram viewer to "SEX"
    Then histogram viewer should have a "line F" area
    When user sets "Bins" property of histogram viewer to "1"
    Then the "bins shown" reading of histogram viewer should be 1
    And histogram viewer should not have a "line F" area
    When user hovers over the "view" area of histogram viewer
    And user moves the pointer away from histogram viewer
    Then the "error" reading of histogram viewer should be ""
    And no errors should have been logged
    When user sets properties of histogram viewer:
      | Bins  | 20 |
      | Split |    |
    Then the "bins shown" reading of histogram viewer should be 20

  # core#998 (c)
  Scenario: A logarithmic X axis over a column with no positive value is reported, not thrown (core#998)
    When user adds a calculated column "NEG_AGE" with formula "0 - ${AGE}"
    And user sets "Value" property of histogram viewer to "NEG_AGE"
    Then histogram viewer should have a "bin 1" area
    When user sets "X Axis Type" property of histogram viewer to "logarithmic"
    Then the "error" reading of histogram viewer should be "no positive values to show on the logarithmic X axis"
    And histogram viewer should not have a "bin 1" area
    When user hovers over the "view" area of histogram viewer
    Then no errors should have been logged
    When user sets "X Axis Type" property of histogram viewer to "linear"
    Then the "error" reading of histogram viewer should be ""
    And histogram viewer should have a "bin 1" area
    When user sets "Value" property of histogram viewer to "AGE"
    And user removes "NEG_AGE" column
    Then no errors should have been logged

  # core#998 (d)
  Scenario: A Y Min above the Y Max is reported, not thrown (core#998)
    When user sets properties of histogram viewer:
      | Y Min | 80 |
      | Y Max | 20 |
    Then the "error" reading of histogram viewer should be "Y Max should be at least 80"
    And histogram viewer should not have a "bin 8" area
    When user hovers over the "view" area of histogram viewer
    Then no errors should have been logged
    When user sets properties of histogram viewer:
      | Y Min |  |
      | Y Max |  |
    Then the "error" reading of histogram viewer should be ""
    And histogram viewer should have a "bin 8" area
    And no errors should have been logged

  # public#2301
  Scenario: The min handle dragged past the max handle leaves a valid range (public#2301)
    When user hovers over the "view" area of histogram viewer
    Then histogram viewer should have a "range min handle" area
    When user drags the "range min handle" area of histogram viewer by 450 pixels to the right
    Then the "range min" reading of histogram viewer should be 89
    And the "range max" reading of histogram viewer should be 89
    And 1 row should pass the filter
    And the "error" reading of histogram viewer should be ""
    And no errors should have been logged
    When user double-clicks on the "range slider" area of histogram viewer
    Then all rows should pass the filter
    And the "range min" reading of histogram viewer should be 18

  # public#130
  Scenario: Split lines take the split column's category colours (public#130)
    When user colors "SEX" column categorically:
      | F | #00FF00 |
      | M | #FF00FF |
    And user sets "Split" property of histogram viewer to "SEX"
    Then histogram viewer should have a "line F" area
    And the "line F" area of histogram viewer should contain the color "#00FF00"
    And the "line M" area of histogram viewer should contain the color "#FF00FF"
    When user sets "Split" property of histogram viewer to ""
    And user removes the coloring of "SEX" column
    Then no errors should have been logged
