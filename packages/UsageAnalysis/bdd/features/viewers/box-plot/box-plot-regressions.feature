@viewers @realizes:viewers.box-plot
Feature: Box plot regressions from closed issues
  Closed box plot issues, each claimed where the bug showed: the zoom survives Show Empty Categories
  (core#307); an explicit Marker Color Column survives a Category 1 change (public#3469); a Shift-drag over a category picks that
  category when the categories to its left are filtered out (public#2764); the selection is
  highlighted under a linear and under a conditional colouring (public#3066); and Show Outside Values
  off takes the outlier markers away (public#1379). Each scenario adds its own box plot.
  core#960's hover over a box plot too small for its axes is claimed in `box-plot-edge-cases.feature`.
  Not translated: core#960's violin over a real NaN — the data steps write NaN as a missing value,
  see MISSING.md; its zero-row and leak halves (an API test, a memory audit).

  Background:
    Given user is logged in

  # core#307
  @realizes:core#307
  Scenario: Toggling Show Empty Categories keeps the zoom (core#307)
    Given user opens spgi dataset
    When user adds a box plot viewer with:
      | Value      | Average Mass |
      | Category 1 | Series       |
    Then "Show Empty Categories" property of box plot viewer should be "true"
    When user zooms into the value axis of box plot viewer
    Then box plot viewer should show a narrower value range than before
    When user remembers the value range of box plot viewer
    And user sets "Show Empty Categories" property of box plot viewer to "false"
    Then box plot viewer should show the remembered value range
    When user sets "Show Empty Categories" property of box plot viewer to "true"
    Then box plot viewer should show the remembered value range
    When user double-clicks on the "view" area of box plot viewer
    Then box plot viewer should show a wider value range than before
    And no errors should have been logged

  # public#3469
  @realizes:public#3469
  Scenario: An explicit Marker Color Column survives a Category 1 change (public#3469)
    Given user opens demog-1000 dataset
    When user adds a box plot viewer with:
      | Value      | AGE |
      | Category 1 | SEX |
    And user sets "Marker Color Column" property of box plot viewer to "HEIGHT"
    Then box plot viewer should have a "color scale" area
    When user sets "Category 1" property of box plot viewer to "RACE"
    Then box plot viewer should have a "category Caucasian" area
    And "Marker Color Column" property of box plot viewer should be "HEIGHT"
    And box plot viewer should have a "color scale" area
    And no errors should have been logged

  # public#2764
  @realizes:public#2764
  Scenario: A Shift-drag picks its category when the categories to its left are filtered out (public#2764)
    Given user opens demog-1000 dataset
    When user adds a box plot viewer with:
      | Value      | AGE  |
      | Category 1 | RACE |
    And user adds a categorical filter on "RACE" keeping "Caucasian, Other"
    Then 958 rows should pass the filter
    And box plot viewer should not have a "category Asian" area
    When user drags a selection box over the "Other values" area of box plot viewer
    Then some rows where "RACE" is "Other" should be selected
    And no rows where "RACE" is "Caucasian" should be selected
    And 51 rows should be selected
    When user clears the row selection
    And user drags a selection box over the "Caucasian values" area of box plot viewer
    Then some rows where "RACE" is "Caucasian" should be selected
    And no rows where "RACE" is "Other" should be selected
    And 768 rows should be selected
    And no errors should have been logged

  # public#3066
  @realizes:public#3066
  Scenario: The selection is highlighted under a linear and under a conditional colouring (public#3066)
    Given user opens demog-1000 dataset
    When user colors "HEIGHT" column linearly from "#00FF00" to "#0000FF"
    And user adds a box plot viewer with:
      | Value               | AGE    |
      | Category 1          | RACE   |
      | Marker Color Column | HEIGHT |
    And user clicks on the "category Caucasian" area of box plot viewer
    Then only rows where "RACE" is "Caucasian" should be selected
    And box plot viewer should show more selection highlight than before
    When user clears the row selection
    Then box plot viewer should show less selection highlight than before
    When user colors "HEIGHT" column conditionally:
      | > 170 | #FF0000 |
    And user clicks on the "category Caucasian" area of box plot viewer
    Then only rows where "RACE" is "Caucasian" should be selected
    And box plot viewer should show more selection highlight than before
    When user clears the row selection
    And user removes the coloring of "HEIGHT" column
    Then no errors should have been logged

  # public#1379
  @realizes:public#1379
  Scenario: Show Outside Values off takes the outlier markers away (public#1379)
    Given user opens demog-1000 dataset
    When user adds a box plot viewer with:
      | Value      | WEIGHT |
      | Category 1 | RACE   |
    And user sets "Show Inside Values" property of box plot viewer to "false"
    And user sets "Show Outside Values" property of box plot viewer to "false"
    Then box plot viewer should have less ink than before
    And box plot viewer should not have a "marker" area
    When user sets "Show Outside Values" property of box plot viewer to "true"
    Then box plot viewer should have more ink than before
    And box plot viewer should have a "marker" area
    And no errors should have been logged
