@journey @viewers @realizes:viewers.density-plot
Feature: Density plot chrome, binding and the viewer filter
  What the plot shows around the bins: the colour scale, the two axes, the three on-viewer
  selectors, the title and the description, and what it draws when it is pointed at another table
  or given a filter of its own.
  Every claim about a piece of chrome is that the plot stops REPORTING it, not that some pixels
  moved: `showColorScale` off sets `_colorScaleBounds` to null, and the axes and selectors are the
  computed auto-layout results, so `should not have a "color scale" area` is a real claim and
  `x selector shown` says what the layout decided rather than what the look asked for.
  The fixture matters and is asserted in the Background: demog-1000 has 1000 rows of which
  **128 have a blank HEIGHT**, and a row with a blank in X or Y is never binned — so a plot on
  AGE × HEIGHT reports `rows shown` 872. That reading is the rows the binning pass counted, not
  `combinedFilter.trueCount`; the scatter plot means the same thing by it.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a density plot viewer with:
      | xColumnName | AGE    |
      | yColumnName | HEIGHT |
    Then 1000 rows should pass the filter
    And the "rows shown" reading of density plot viewer should be 872
    And the "x column" reading of density plot viewer should be "AGE"
    And the "y column" reading of density plot viewer should be "HEIGHT"
    And density plot viewer should be painted

  Scenario: A blank in either column keeps the row out of the bins
    Then the "rows shown" reading of density plot viewer should be 872
    When user sets "yColumnName" property of density plot viewer to "WEIGHT"
    Then the "rows shown" reading of density plot viewer should be 1000
    And 1000 rows should pass the filter
    When user sets "yColumnName" property of density plot viewer to "HEIGHT"
    Then the "rows shown" reading of density plot viewer should be 872
    And no errors should have been logged

  Scenario: Show Color Scale decides whether the scale is drawn at all
    Then density plot viewer should have a "color scale" area
    And the "color scale" area of density plot viewer should be painted
    When user sets "showColorScale" property of density plot viewer to "false"
    Then density plot viewer should not have a "color scale" area
    And density plot viewer should have less ink than before
    When user sets "showColorScale" property of density plot viewer to "true"
    Then density plot viewer should have a "color scale" area
    And no errors should have been logged

  Scenario: The colour scale bounds are reported whether or not the scale is drawn
    Then the "color scale min" reading of density plot viewer should be 0
    And the "color scale max" and "max bin count" readings of density plot viewer should be the same
    When user sets "showColorScale" property of density plot viewer to "false"
    Then the "color scale max" and "max bin count" readings of density plot viewer should be the same
    When user sets "showColorScale" property of density plot viewer to "true"
    Then no errors should have been logged

  Scenario: Show X Axis and Show Y Axis take the axis boxes away
    Then density plot viewer should have an "x axis" area
    And density plot viewer should have a "y axis" area
    When user sets "showXAxis" property of density plot viewer to "false"
    Then density plot viewer should not have an "x axis" area
    And density plot viewer should have a "y axis" area
    When user sets "showYAxis" property of density plot viewer to "false"
    Then density plot viewer should not have a "y axis" area
    When user sets properties of density plot viewer:
      | showXAxis | true |
      | showYAxis | true |
    Then density plot viewer should have an "x axis" area
    And density plot viewer should have a "y axis" area
    And no errors should have been logged

  Scenario: Show X Selector and Show Y Selector are read as the layout resolved them
    Then the "x selector shown" reading of density plot viewer should be "true"
    And the "y selector shown" reading of density plot viewer should be "true"
    When user sets "showXSelector" property of density plot viewer to "false"
    Then the "x selector shown" reading of density plot viewer should be "false"
    And the "y selector shown" reading of density plot viewer should be "true"
    When user sets "showYSelector" property of density plot viewer to "false"
    Then the "y selector shown" reading of density plot viewer should be "false"
    When user sets properties of density plot viewer:
      | showXSelector | true |
      | showYSelector | true |
    Then the "x selector shown" reading of density plot viewer should be "true"
    And no errors should have been logged

  Scenario: The description sits above the bins and Description Position moves it below them
    When user sets "description" property of density plot viewer to "Age against height"
    Then the description of density plot viewer should be above its content
    When user sets "descriptionPosition" property of density plot viewer to "Bottom"
    Then the description of density plot viewer should be below its content
    When user sets properties of density plot viewer:
      | descriptionPosition | Top |
      | description         |     |
    Then no errors should have been logged

  Scenario: The viewer's own filter narrows what it bins and leaves the table alone
    When user sets "filter" property of density plot viewer to "${AGE} > 30"
    Then the "rows shown" reading of density plot viewer should be lower than before
    And 1000 rows should pass the filter
    And density plot viewer should have repainted
    When user sets "filter" property of density plot viewer to ""
    Then the "rows shown" reading of density plot viewer should be 872
    And no errors should have been logged

  Scenario: A filter on the table moves what the plot bins
    When user adds a categorical filter on "SEX" keeping "M"
    Then 447 rows should pass the filter
    And the "rows shown" reading of density plot viewer should be lower than before
    And density plot viewer should have repainted
    When user hovers over "SEX" filter card
    And user clicks on close of "SEX" filter card
    Then 1000 rows should pass the filter
    And the "rows shown" reading of density plot viewer should be 872
    And no errors should have been logged

  Scenario: Bound to another table, the plot bins that table's rows
    Given user opens spgi dataset
    And user switches to the "demog-1000" table view
    When user sets "table" property of density plot viewer to "spgi-100"
    Then density plot viewer should be bound to table "spgi-100"
    And the "rows shown" reading of density plot viewer should be at least 1
    And the "rows shown" reading of density plot viewer should differ from before
    And density plot viewer should be painted
    And no errors should have been logged
