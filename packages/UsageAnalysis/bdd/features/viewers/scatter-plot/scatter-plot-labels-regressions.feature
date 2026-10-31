@viewers @realizes:viewers.scatter-plot
Feature: Scatter plot regressions — dragging a label
  Under Mouse Drag = Select a drag that starts on a marker label moves the label and selects nothing:
  with Show Labels For = Selected (public#3120) the dragged label leaves the auto-placed ones, so
  `labels shown` drops from 1 to 0 and the selection stays; with Show Labels For = All under a panel
  filter (public#3323) the plot places the other labels anew (274 before the drag, 301 after,
  measured; the rise is intended) and no row is selected. Each scenario draws its own plot, so one
  drag cannot move the other's label. The round's other scatter plot issues are in the plot's
  journeys and in viewer-chrome.feature.

  Background:
    Given user is logged in

  # public#3120
  @realizes:public#3120
  Scenario: A label is dragged under Mouse Drag = Select (public#3120)
    Given user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | X               | WEIGHT   |
      | Y               | HEIGHT   |
      | Label Columns   | AGE      |
      | Show Labels For | Selected |
      | mouseDrag       | Select   |
    When user selects the first 1 rows
    Then scatter plot viewer should have a "label of row 1" area
    And the "labels shown" reading of scatter plot viewer should be 1
    When user drags the "label of row 1" area of scatter plot viewer by 60 pixels to the right
    Then 1 row should be selected
    And the "labels shown" reading of scatter plot viewer should be 0
    And no errors should have been logged

  # public#3323
  @realizes:public#3323
  Scenario: A label is dragged when Show Labels For is All and a panel filter is set (public#3323)
    Given user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | X               | WEIGHT |
      | Y               | HEIGHT |
      | Label Columns   | AGE    |
      | Show Labels For | All    |
      | mouseDrag       | Select |
    When user adds a categorical filter on "SEX" keeping "F"
    Then 553 rows should pass the filter
    And scatter plot viewer should have a "label of row 1" area
    And the "labels shown" reading of scatter plot viewer should be 274
    When user drags the "label of row 1" area of scatter plot viewer by 60 pixels to the right
    Then the "labels shown" reading of scatter plot viewer should be 301
    And no rows should be selected
    And no errors should have been logged
