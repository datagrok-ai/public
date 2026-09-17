@journey @viewers @realizes:viewers.scatter-plot
Feature: Scatter plot legend
  What the scatter plot puts in its legend and what the legend does back to the plot: a categorical
  Color builds it, Markers adds a second section to the same legend, a numeric Color leaves only
  the marker section, clearing Markers takes it away again, a table filter shrinks the marker
  section and giving the filter back restores it, a click on an entry hides that category on the
  canvas without touching the table filter, and Legend Visibility and Legend Position decide
  whether and where it shows. The item counts are what the legend lists (every section, rendered
  or not — a corner legend renders only the rows that fit). The item-level rendering — glyphs,
  colors, the cross — belongs to the legend's own feature. One journey on demog-1000, X = WEIGHT,
  Y = HEIGHT (872 of the 1000 rows have a HEIGHT); every scenario puts back what it changed.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | X | WEIGHT |
      | Y | HEIGHT |
    Then scatter plot viewer should show 872 rows
    And legend of scatter plot viewer should be hidden

  Scenario: Color and Markers build one legend of two sections
    When user sets "Color" property of scatter plot viewer to "RACE"
    Then legend of scatter plot viewer should be visible
    And the legend of scatter plot viewer should list 4 items
    And legend of scatter plot viewer should contain text "Asian"
    And legend of scatter plot viewer should contain text "Caucasian"
    When user sets "Markers" property of scatter plot viewer to "SEX"
    Then the legend of scatter plot viewer should list 6 items
    When user sets "Color" property of scatter plot viewer to "AGE"
    Then the legend of scatter plot viewer should list 2 items
    When user sets "Color" property of scatter plot viewer to "RACE"
    Then the legend of scatter plot viewer should list 6 items
    And legend of scatter plot viewer should contain text "Asian"
    And no errors should have been logged

  Scenario: Clearing Markers leaves the color entries alone
    When user sets properties of scatter plot viewer:
      | Color   | SEX |
      | Markers | SEX |
    Then the legend of scatter plot viewer should list 2 items
    When user sets "Markers" property of scatter plot viewer to ""
    Then "Markers" property of scatter plot viewer should be ""
    And the legend of scatter plot viewer should list 2 items
    And legend of scatter plot viewer should contain text "F"
    When user sets properties of scatter plot viewer:
      | Color   | RACE |
      | Markers | SEX  |
    Then the legend of scatter plot viewer should list 6 items
    And no errors should have been logged

  Scenario: A table filter drops the filtered-out categories and Markers on the color column adds no section
    When user sets "Markers" property of scatter plot viewer to ""
    Then the legend of scatter plot viewer should list 4 items
    When user opens the filter panel
    And user adds a categorical filter on "RACE" keeping "Asian, Caucasian"
    Then 911 rows should pass the filter
    And the legend of scatter plot viewer should list 2 items
    When user sets "Markers" property of scatter plot viewer to "RACE"
    Then the legend of scatter plot viewer should list 2 items
    When user sets "Markers" property of scatter plot viewer to "SEX"
    Then the legend of scatter plot viewer should list 4 items
    When user resets the filter
    Then all rows should pass the filter
    And the legend of scatter plot viewer should list 6 items
    When user sets "Markers" property of scatter plot viewer to "RACE"
    Then the legend of scatter plot viewer should list 4 items
    When user sets "Markers" property of scatter plot viewer to "SEX"
    Then the legend of scatter plot viewer should list 6 items
    And no errors should have been logged

  Scenario: A click on a legend entry hides that category on the canvas, not in the table
    Then all rows should pass the filter
    And scatter plot viewer should show 872 rows
    When user clicks on "Asian" item in the legend of scatter plot viewer
    Then scatter plot viewer should show fewer rows than before
    And scatter plot viewer should have less ink than before
    And all rows should pass the filter
    When user clicks on "Asian" item in the legend of scatter plot viewer
    Then scatter plot viewer should show 872 rows
    And all rows should pass the filter
    And no errors should have been logged

  Scenario: Legend Visibility and Legend Position decide whether and where it shows
    When user sets "Legend Visibility" property of scatter plot viewer to "Never"
    Then legend of scatter plot viewer should be hidden
    When user sets "Legend Visibility" property of scatter plot viewer to "Always"
    Then legend of scatter plot viewer should be visible
    And the legend of scatter plot viewer should list 6 items
    When user sets "Legend Position" property of scatter plot viewer to "Left"
    Then the legend of scatter plot viewer should be on the left
    When user sets "Legend Position" property of scatter plot viewer to "Right"
    Then the legend of scatter plot viewer should be on the right
    When user sets properties of scatter plot viewer:
      | Legend Position   | Auto |
      | Legend Visibility | Auto |
      | Color             |      |
      | Markers           |      |
    Then legend of scatter plot viewer should be hidden
    And scatter plot viewer should show 872 rows
    And no errors should have been logged
