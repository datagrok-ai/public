@journey @viewers @realizes:viewers.legend
Feature: Legend regressions on the scatter plot
  Closed legend issues turned into claims on the scatter plot's legend: a viewer Filter that empties
  the plot and then lets rows back brings the legend back with only the categories left, and a saved
  layout keeps it that way (public#1209); a viewer Filter set after Markers drops the filtered-out
  categories from the marker section (public#2505); with Color and Markers on the same column a
  color can still be picked from the legend (public#2965); and a Color column with a thousand
  categories leaves the Markers section on the screen (public#3289). One journey on demog-1000,
  X = WEIGHT, Y = HEIGHT (RACE: Caucasian, Other, Black, Asian; the viewer Filter keeps Asian and
  Black); every scenario puts back what it changed, and the color picked in the legend is taken off
  the RACE column again.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | X | WEIGHT |
      | Y | HEIGHT |
    Then scatter plot viewer should show 872 rows

  # public#1209
  @realizes:public#1209
  Scenario: A viewer Filter that empties the plot and lets rows back brings the legend back (public#1209)
    When user sets properties of scatter plot viewer:
      | Color             | RACE |
      | Legend Visibility | Auto |
    Then the legend of scatter plot viewer should list 4 items
    When user sets "Filter" property of scatter plot viewer to '${RACE} == "none"'
    Then scatter plot viewer should show 0 rows
    When user sets "Filter" property of scatter plot viewer to '${RACE} in ["Asian", "Black"]'
    Then scatter plot viewer should show 38 rows
    And the legend of scatter plot viewer should list 2 items
    And "Asian" legend item in legend of scatter plot viewer should be visible
    And "Caucasian" legend item in legend of scatter plot viewer should be absent
    When user saves the layout of the current table view
    And user loads the saved layout
    Then scatter plot viewer should show 38 rows
    And the legend of scatter plot viewer should list 2 items
    And "Black" legend item in legend of scatter plot viewer should be visible
    And "Caucasian" legend item in legend of scatter plot viewer should be absent
    When user sets properties of scatter plot viewer:
      | Filter |  |
      | Color  |  |
    Then scatter plot viewer should show 872 rows
    And no errors should have been logged

  # public#2505
  @realizes:public#2505
  Scenario: A viewer Filter set after Markers drops the filtered-out categories from the marker legend (public#2505)
    When user sets "Markers" property of scatter plot viewer to "RACE"
    Then the legend of scatter plot viewer should list 4 items
    When user sets "Filter" property of scatter plot viewer to '${RACE} in ["Asian", "Black"]'
    Then scatter plot viewer should show 38 rows
    And the legend of scatter plot viewer should list 2 items
    And "Caucasian" legend item in legend of scatter plot viewer should be absent
    And "Other" legend item in legend of scatter plot viewer should be absent
    When user sets properties of scatter plot viewer:
      | Filter  |  |
      | Markers |  |
    Then scatter plot viewer should show 872 rows
    And no errors should have been logged

  # public#2965
  @realizes:public#2965
  Scenario: With Color and Markers on the same column a color is picked from the legend (public#2965)
    When user sets properties of scatter plot viewer:
      | Color             | RACE   |
      | Markers           | RACE   |
      | Legend Visibility | Always |
    Then the legend of scatter plot viewer should list 4 items
    When user hovers over "Asian" legend item in legend of scatter plot viewer
    And user clicks on color picker icon
    Then "Asian" dialog should be visible
    When user picks the color "#9467BD" in the color picker dialog
    And user clicks on OK button in "Asian" dialog
    Then "Asian" dialog should be absent
    And the categorical color of "Asian" in "RACE" column should be "#9467BD"
    And the "Asian" item in the legend of scatter plot viewer should be colored "#9467BD"
    When user removes the coloring of "RACE" column
    And user sets properties of scatter plot viewer:
      | Color             |      |
      | Markers           |      |
      | Legend Visibility | Auto |
    Then "RACE" column should have no color coding
    And no errors should have been logged

  # public#3289
  @realizes:public#3289
  Scenario: A Color column of a thousand categories leaves the Markers section on the screen (public#3289)
    When user sets properties of scatter plot viewer:
      | Color             | USUBJID |
      | Markers           | SEX     |
      | Legend Visibility | Always  |
    Then the legend of scatter plot viewer should be docked
    And "F" legend item in legend of scatter plot viewer should be visible
    And "M" legend item in legend of scatter plot viewer should be visible
    When user sets properties of scatter plot viewer:
      | Color             |      |
      | Markers           |      |
      | Legend Visibility | Auto |
    Then scatter plot viewer should show 872 rows
    And no errors should have been logged
