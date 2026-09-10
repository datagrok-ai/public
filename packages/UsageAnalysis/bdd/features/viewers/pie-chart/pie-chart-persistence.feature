@journey @viewers @realizes:viewers.pie-chart
Feature: Pie chart persistence
  A configured pie chart survives a layout saved on the server and a project saved and reopened:
  the properties come back, the disc is drawn from them again, the saved layout restores the viewer
  set it was saved with (a viewer added afterwards is gone), and the custom category colour the
  legend dialog wrote into the column travels with both. One journey on demog-1000 with a pie of
  RACE by sum(AGE), started at 45 degrees, shifted 5 pixels, showing values, titled, with Asian
  recoloured — sum(AGE) makes Caucasian 89.59 % of the disc against 89.6 % under count, so the
  restored chart is the configured one and not a default. The coloring is put back by naming
  Asian's default colour rather than by Color coding > Off: turning the coding off leaves the
  category map on the column and the pie goes on painting from it.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a pie chart viewer with:
      | Category                | RACE     |
      | Segment Angle Aggr Type | sum      |
      | Show Value              | true     |
      | Start Angle             | 45       |
      | Shift                   | 5        |
      | Legend Visibility       | Always   |
      | Show Title              | true     |
      | Title                   | Race mix |
    And user colors "RACE" column categorically:
      | Asian | #9467BD |
    Then the "slices" reading of pie chart viewer should be 4
    And the "start angle of Asian" reading of pie chart viewer should be 45
    And the "share of Caucasian" reading of pie chart viewer should be between 89.55 and 89.6

  Scenario: The configured chart is the one on screen
    Then title of pie chart viewer should have text "Race mix"
    And the slices of pie chart viewer should sit 5 pixels off the centre
    And "RACE" column should be color-coded categorically
    And the "Asian" item in the legend of pie chart viewer should be colored "#9467BD"
    And the "pie" area of pie chart viewer should contain the color "#9467BD"
    And no errors should have been logged

  Scenario: A layout saved on the server restores the configuration and the viewer set
    When user saves the layout of the current table view to the server
    And user clicks on close icon of pie chart viewer
    Then pie chart viewer should be absent
    When user adds a scatter plot viewer
    Then scatter plot viewer should be visible
    When user loads the saved layout
    Then pie chart viewer should be visible
    And scatter plot viewer should be absent
    And properties of pie chart viewer should be:
      | Category                | RACE     |
      | Segment Angle Aggr Type | sum      |
      | Show Value              | true     |
      | Start Angle             | 45       |
      | Shift                   | 5        |
      | Title                   | Race mix |
    And the "start angle of Asian" reading of pie chart viewer should be 45
    And the "share of Caucasian" reading of pie chart viewer should be between 89.55 and 89.6
    And the slices of pie chart viewer should sit 5 pixels off the centre
    And no errors should have been logged

  Scenario: The custom category colour comes back with the layout
    Then "RACE" column should be color-coded categorically
    And the categorical color of "Asian" in "RACE" column should be "#9467BD"
    And the "pie" area of pie chart viewer should contain the color "#9467BD"
    And no errors should have been logged

  Scenario: A project saved, closed and reopened brings all of it back
    When user saves the current view as project "bdd pie chart round trip"
    And user closes all views
    And user opens the "bdd pie chart round trip" project
    Then pie chart viewer should be visible
    And properties of pie chart viewer should be:
      | Category                | RACE     |
      | Segment Angle Aggr Type | sum      |
      | Show Value              | true     |
      | Start Angle             | 45       |
      | Shift                   | 5        |
      | Title                   | Race mix |
    And the "slices" reading of pie chart viewer should be 4
    And the "start angle of Asian" reading of pie chart viewer should be 45
    And the "share of Caucasian" reading of pie chart viewer should be between 89.55 and 89.6
    And "RACE" column should be color-coded categorically
    And the "pie" area of pie chart viewer should contain the color "#9467BD"
    And no errors should have been logged

  Scenario: The coloring goes back where it was found
    When user colors "RACE" column categorically:
      | Asian | #1F77B4 |
    Then the "pie" area of pie chart viewer should contain the color "#1F77B4"
    And the "pie" area of pie chart viewer should not contain the color "#9467BD"
    When user removes the coloring of "RACE" column
    Then "RACE" column should have no color coding
    And no errors should have been logged
