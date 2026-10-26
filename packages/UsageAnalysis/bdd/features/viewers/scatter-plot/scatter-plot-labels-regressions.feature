@viewers @realizes:viewers.scatter-plot
Feature: Scatter plot regressions — labels and missing columns
  Closed scatter plot issues about marker labels and about a plot that lost a column: a label is
  dragged, not the plot under it, when Mouse Drag = Select (public#3120) and when Show Labels For is
  All under a panel filter (public#3323) — a dragged label leaves the auto-placed ones, so `labels
  shown` drops by one, and the selection stays as it was; a molecule column as the label column draws
  labels (public#3041, spgi-100); and a layout applied after its logarithmic X column was removed
  loads, and the view's layout can be saved and loaded again without an error (public#2867). Each
  scenario opens its own table. The plot draws "x: Column "WEIGHT" does not exist." then, but its
  status carries no `error`, so the message itself is in MISSING-regressions.md (public#3091).

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
    When user remembers the "labels shown" reading of scatter plot viewer
    And user drags the "label of row 1" area of scatter plot viewer by 60 pixels to the right
    Then the "labels shown" reading of scatter plot viewer should not be as remembered
    And no rows should be selected
    And no errors should have been logged

  # public#3041, public#3040
  @realizes:public#3041 @realizes:public#3040
  Scenario: A molecule column as the label column draws labels (public#3041)
    Given user opens spgi-100 dataset
    And user adds a scatter plot viewer with:
      | X | Stereo Category |
      | Y | Chemist 521     |
    Then the "labels shown" reading of scatter plot viewer should be 0
    When user sets "Label Columns" property of scatter plot viewer to "Core"
    Then the "labels shown" reading of scatter plot viewer should be at least 1
    And scatter plot viewer should report no error
    And no errors should have been logged

  # public#2867
  @realizes:public#2867
  Scenario: A plot whose logarithmic X column is removed lets the layout be saved (public#2867)
    Given user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | X           | WEIGHT      |
      | Y           | HEIGHT      |
      | X Axis Type | logarithmic |
      | Y Axis Type | logarithmic |
    When user saves the layout of the current table view
    And user removes "WEIGHT" column
    And user loads the saved layout
    Then scatter plot viewer should be visible
    And "X" property of scatter plot viewer should be "WEIGHT"
    When user saves the layout of the current table view
    And user loads the saved layout
    Then scatter plot viewer should be visible
    And "X" property of scatter plot viewer should be "WEIGHT"
    And "X Axis Type" property of scatter plot viewer should be "logarithmic"
    And no errors should have been logged
    And no error or warning balloon should have been shown
