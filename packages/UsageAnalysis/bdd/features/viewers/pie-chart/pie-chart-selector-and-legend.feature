@journey @viewers @realizes:viewers.pie-chart
Feature: Pie chart on-chart selector and legend
  The category column can be re-picked on the chart itself, and the legend that goes with it: the
  selector re-splits the disc and the legend follows it, the legend lists exactly the column's
  categories, Legend Position moves it to each of the four sides and Legend Visibility takes it
  away, and a colour set on a legend item through its own right-click dialog reaches the wedge and
  the column's own categorical coloring. One journey on demog-1000 — RACE has 4 categories
  (Caucasian 896, Other 62, Black 27, Asian 15), SEX has 2 (F 553, M 447), and the default palette
  gives Asian #1F77B4, Black #FFBB78, Caucasian #2CA02C, Other #D62728.
  Not claimed here, because the product does not do it: the dialog does NOT set the column's
  `.color-coding-type` to Categorical (only the category map is written, so the platform's own
  coding-type tag stays unset), and Color coding > Off leaves the map behind — the pie keeps
  painting the custom colours. So the wedge is taken back the way it was set, through the dialog.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a pie chart viewer with:
      | Category          | RACE   |
      | Legend Visibility | Always |
    Then the "slices" reading of pie chart viewer should be 4
    And legend of pie chart viewer should be visible

  Scenario: The legend lists exactly the category column's categories
    Then the legend of pie chart viewer should list 4 items
    And "Asian" legend item in legend of pie chart viewer should be visible
    And "Black" legend item in legend of pie chart viewer should be visible
    And "Caucasian" legend item in legend of pie chart viewer should be visible
    And "Other" legend item in legend of pie chart viewer should be visible
    And the "Asian" and "Caucasian" items in the legend of pie chart viewer should be colored differently
    And no errors should have been logged

  Scenario: Picking SEX in the on-chart selector re-splits the disc and the legend
    When user hovers over pie chart viewer
    And user picks "SEX" in the category selector of pie chart viewer
    Then "Category" property of pie chart viewer should be "SEX"
    And the "slices" reading of pie chart viewer should be 2
    And pie chart viewer should have a "slice F" area
    And pie chart viewer should have a "slice M" area
    And the "angle value of F" reading of pie chart viewer should be 553
    And the "angle value of M" reading of pie chart viewer should be 447
    And the legend of pie chart viewer should list 2 items
    And pie chart viewer should have repainted by at least 1000 pixels
    And no errors should have been logged

  Scenario: Picking RACE back restores it
    When user hovers over pie chart viewer
    And user picks "RACE" in the category selector of pie chart viewer
    Then "Category" property of pie chart viewer should be "RACE"
    And the "slices" reading of pie chart viewer should be 4
    And the "angle value of Caucasian" reading of pie chart viewer should be 896
    And pie chart viewer should not have a "slice F" area
    And the legend of pie chart viewer should list 4 items
    And pie chart viewer should have repainted by at least 1000 pixels
    And no errors should have been logged

  Scenario: Legend Position lays the legend out on each of the four sides
    When user sets "Legend Position" property of pie chart viewer to "Left"
    Then the legend of pie chart viewer should be on the left
    And the "pie radius" reading of pie chart viewer should be lower than before
    When user sets "Legend Position" property of pie chart viewer to "Right"
    Then the legend of pie chart viewer should be on the right
    When user sets "Legend Position" property of pie chart viewer to "Top"
    Then the legend of pie chart viewer should be on the top
    When user sets "Legend Position" property of pie chart viewer to "Bottom"
    Then the legend of pie chart viewer should be on the bottom
    And legend of pie chart viewer should be visible
    When user sets "Legend Position" property of pie chart viewer to "Auto"
    Then legend of pie chart viewer should be visible
    And no errors should have been logged

  Scenario: Legend Visibility Never takes it away and Always brings it back
    Then legend of pie chart viewer should be visible
    When user sets "Legend Visibility" property of pie chart viewer to "Never"
    Then legend of pie chart viewer should be hidden
    When user sets "Legend Visibility" property of pie chart viewer to "Always"
    Then legend of pie chart viewer should be visible
    And the legend of pie chart viewer should list 4 items
    And no errors should have been logged

  Scenario: A colour picked on a legend item reaches the wedge and the column
    Given the "pie" area of pie chart viewer should not contain the color "#9467BD"
    And the "Asian" item in the legend of pie chart viewer should be colored "#1F77B4"
    When user right-clicks on "Asian" legend item in legend of pie chart viewer
    Then "Asian" dialog should be visible
    When user picks the color "#9467BD" in the color picker dialog
    And user clicks on OK button in "Asian" dialog
    Then the categorical color of "Asian" in "RACE" column should be "#9467BD"
    And the "Asian" item in the legend of pie chart viewer should be colored "#9467BD"
    And the "pie" area of pie chart viewer should contain the color "#9467BD"
    And the "pie" area of pie chart viewer should not contain the color "#1F77B4"
    And no errors should have been logged

  Scenario: Picking the default colour back takes the wedge back
    When user right-clicks on "Asian" legend item in legend of pie chart viewer
    And user picks the color "#1F77B4" in the color picker dialog
    And user clicks on OK button in "Asian" dialog
    Then the categorical color of "Asian" in "RACE" column should be "#1F77B4"
    And the "Asian" item in the legend of pie chart viewer should be colored "#1F77B4"
    And the "pie" area of pie chart viewer should contain the color "#1F77B4"
    And the "pie" area of pie chart viewer should not contain the color "#9467BD"
    When user removes the coloring of "RACE" column
    Then "RACE" column should have no color coding
    And no errors should have been logged
