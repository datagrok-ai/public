@journey @viewers @realizes:viewers.scatter-plot
Feature: Scatter plot — the "(no value)" legend item as a click target
  `scatter-plot-legend.feature` shows that a categorical Color with empty values lists them as an
  item of their own; this one clicks that item. The click narrows the plot to the rows whose colour
  value is empty and leaves the table's filter alone, and a filter-panel card on top narrows those
  rows further without bringing any of the hidden categories back; a second click on the item gives
  the plot every row the filter passes again. The empty values come from a calculated column that
  blanks DIS_POP where RACE is Other. One journey on demog-1000, X = WEIGHT, Y = HEIGHT (872 of the
  1000 rows have a HEIGHT).
  Not translated: GROK-17282, an unsupported line formula (`Y - X = 1`) in the Formula Lines dialog,
  and the dimmed look of an item the click took out — see MISSING.md.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | X | WEIGHT |
      | Y | HEIGHT |
    Then scatter plot viewer should show 872 rows

  Scenario: A click on "(no value)" draws only the rows with no colour value, and a panel card narrows them further
    When user adds a calculated column "DIS_POP but Other" with formula "if(${RACE}=='Other', null, ${DIS_POP})"
    And user sets "Color" property of scatter plot viewer to "DIS_POP but Other"
    Then the legend of scatter plot viewer should list 6 items
    And "(no value)" legend item in legend of scatter plot viewer should be visible
    When user clicks on "(no value)" item in the legend of scatter plot viewer
    Then scatter plot viewer should show 59 rows
    And all rows should pass the filter
    And "(no value)" legend item in legend of scatter plot viewer should be selected
    When user adds a categorical filter on "SEX" keeping "M"
    Then 447 rows should pass the filter
    And scatter plot viewer should show 12 rows
    When user clicks on "(no value)" item in the legend of scatter plot viewer
    Then scatter plot viewer should show 373 rows
    And 447 rows should pass the filter
    When user hovers over "SEX" filter card
    And user clicks on close of "SEX" filter card
    Then all rows should pass the filter
    And scatter plot viewer should show 872 rows
    When user sets "Color" property of scatter plot viewer to ""
    And user removes "DIS_POP but Other" column
    Then legend of scatter plot viewer should be hidden
    And no errors should have been logged
