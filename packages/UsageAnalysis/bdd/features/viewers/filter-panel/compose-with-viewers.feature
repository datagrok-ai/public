@journey @viewers @realizes:viewers.filters
Feature: The panel's criterion composes with the viewers
  A card's criterion and a viewer's own filtering intersect and neither loses the other: a scatter
  plot zoom narrows the rows the card left and Reset View gives them back, a bar click under
  On Click Filter keeps the bar's category inside the card's, a histogram range narrows further and
  closing the viewer releases only its share, and the grid shows exactly the rows the card keeps.
  One journey on demog-1000 with a RACE card keeping Caucasian — 896 rows, 416 of them M, 633 of
  them aged 30 to 60.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user opens an empty filter panel
    And user adds a categorical filter on "RACE" keeping "Caucasian"
    Then 896 rows should pass the filter
    And the "selected categories of RACE" reading of filter panel should be "Caucasian"
    And counter of filter panel should have text "1"

  Scenario: A scatter plot zoom narrows the rows the card left and Reset View gives them back
    When user adds a scatter plot viewer with:
      | X | AGE    |
      | Y | HEIGHT |
    Then "Zoom and Filter" property of scatter plot viewer should be "filter by zoom"
    When user scrolls the mouse wheel up over the "view" area of scatter plot viewer
    Then fewer than 896 rows should pass the filter
    And no rows where "RACE" is "Black" should pass the filter
    And the "rows shown" reading of filter panel should be at least 1
    And the "selected categories of RACE" reading of filter panel should be "Caucasian"
    And counter of filter panel should have text "1"
    When user picks "Reset View" from the context menu of scatter plot viewer
    Then 896 rows should pass the filter
    When user clicks on close icon of scatter plot viewer
    Then 896 rows should pass the filter
    And no errors should have been logged

  Scenario: A bar click keeps the bar's category inside the card's
    When user adds a bar chart viewer with:
      | Split | SEX |
    And user sets "On Click" property of bar chart viewer to "Filter"
    And user clicks on the "bar M" area of bar chart viewer
    Then 416 rows should pass the filter
    And the "rows shown" reading of grid should be 416
    And the "selected categories of RACE" reading of filter panel should be "Caucasian"
    And counter of filter panel should have text "1"
    When user clicks on close icon of bar chart viewer
    Then 896 rows should pass the filter
    And no errors should have been logged

  Scenario: A histogram range narrows the card's rows further
    When user adds a histogram viewer with:
      | Value             | AGE  |
      | Show Range Inputs | true |
      | Filtering Enabled | true |
    And user resizes histogram viewer to 500 by 400
    And user enters "30" into the "range min input" area of histogram viewer
    Then 771 rows should pass the filter
    When user enters "60" into the "range max input" area of histogram viewer
    Then 633 rows should pass the filter
    And the "selected categories of RACE" reading of filter panel should be "Caucasian"
    And counter of filter panel should have text "1"
    When user clicks on close icon of histogram viewer
    Then 896 rows should pass the filter
    And no errors should have been logged

  Scenario: The grid shows exactly the rows the card keeps
    Then the "rows shown" reading of grid should be 896
    When user clicks on the "category Asian of RACE" area of filter panel
    Then 15 rows should pass the filter
    And the "rows shown" reading of grid should be 15
    When user clicks on the "category Caucasian of RACE" area of filter panel
    Then 896 rows should pass the filter
    And the "rows shown" reading of grid should be 896
    And no errors should have been logged
