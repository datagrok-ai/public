@journey @viewers @realizes:viewers.pc-plot
Feature: PC plot colouring, legend and colour scale
  Colouring the polylines by a column. A categorical column gets a legend listing its categories; a
  numerical one gets no legend but the colour scale the plot draws on its overlay (the `color scale`
  hit area), which the Color Scheme menu inverts and edits and which Color Min / Color Max clamp. A
  column the grid colours conditionally hands the plot its bins instead, and a DateTime column with
  a Color Map is split into the categories that map names.
  demog-1000: RACE has exactly Caucasian, Asian, Black and Other; HEIGHT has 128 blanks, so a
  conditional coding on it shows its two bins plus "no value"; STARTED spans 1989..1991, so the
  colour map yields 3 years, 4 quarters and 12 months. Every scenario clears the colour column it
  set.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a pc plot viewer with:
      | Column Names | AGE, HEIGHT, WEIGHT |
    Then pc plot viewer should show 1000 rows
    And "Color" property of pc plot viewer should be ""
    And legend of pc plot viewer should be hidden
    And pc plot viewer should not have a "color scale" area

  Scenario: A categorical colour column lists its categories in the legend
    When user sets "Color" property of pc plot viewer to "RACE"
    Then "Color" property of pc plot viewer should be "RACE"
    And legend of pc plot viewer should be visible
    And the legend of pc plot viewer should list 4 items
    And "Caucasian" legend item in legend of pc plot viewer should be visible
    And "Asian" legend item in legend of pc plot viewer should be visible
    And "Black" legend item in legend of pc plot viewer should be visible
    And "Other" legend item in legend of pc plot viewer should be visible
    And the "Caucasian" and "Asian" items in the legend of pc plot viewer should be colored differently
    And pc plot viewer should have repainted by at least 2000 pixels
    And pc plot viewer should not have a "color scale" area
    When user sets "Color" property of pc plot viewer to ""
    Then legend of pc plot viewer should be hidden
    And pc plot viewer should have repainted by at least 2000 pixels
    And no errors should have been logged

  Scenario: Legend Position moves the legend and Legend Visibility takes it away
    When user sets "Color" property of pc plot viewer to "RACE"
    Then the legend of pc plot viewer should list 4 items
    When user sets "Legend Position" property of pc plot viewer to "Left"
    Then the legend of pc plot viewer should be on the left
    And the legend of pc plot viewer should list the same items as before
    When user sets "Legend Position" property of pc plot viewer to "Right"
    Then the legend of pc plot viewer should be on the right
    And the legend of pc plot viewer should list the same items as before
    When user sets "Legend Position" property of pc plot viewer to "Top"
    Then the legend of pc plot viewer should be on the top
    When user sets "Legend Position" property of pc plot viewer to "Bottom"
    Then "Legend Position" property of pc plot viewer should be "Bottom"
    And the legend of pc plot viewer should list 4 items
    When user sets "Legend Visibility" property of pc plot viewer to "Never"
    Then legend of pc plot viewer should be hidden
    When user sets "Legend Visibility" property of pc plot viewer to "Auto"
    Then legend of pc plot viewer should be visible
    And the legend of pc plot viewer should list 4 items
    When user sets properties of pc plot viewer:
      | Legend Position | Auto |
      | Color           |      |
    Then legend of pc plot viewer should be hidden
    And no errors should have been logged

  Scenario: A numerical colour column draws the colour scale instead of a legend
    When user sets "Color" property of pc plot viewer to "AGE"
    Then legend of pc plot viewer should be hidden
    And pc plot viewer should have a "color scale" area
    And the "color scale" area of pc plot viewer should be painted
    And the "color scale" area of pc plot viewer should contain the color "#FF0000"
    And pc plot viewer should have repainted by at least 2000 pixels
    When user sets "Color" property of pc plot viewer to ""
    Then pc plot viewer should not have a "color scale" area
    And no errors should have been logged

  Scenario: Color Scheme > Invert Color Scheme repaints the scale and the lines
    When user sets "Color" property of pc plot viewer to "AGE"
    Then "Invert Color Scheme" property of pc plot viewer should be "false"
    When user picks "Color Scheme > Invert Color Scheme" from the context menu of pc plot viewer
    Then "Invert Color Scheme" property of pc plot viewer should be "true"
    And the "color scale" area of pc plot viewer should have repainted
    And pc plot viewer should have repainted by at least 2000 pixels
    When user picks "Color Scheme > Invert Color Scheme" from the context menu of pc plot viewer
    Then "Invert Color Scheme" property of pc plot viewer should be "false"
    And the "color scale" area of pc plot viewer should have repainted
    And pc plot viewer should have repainted by at least 2000 pixels
    And no errors should have been logged

  Scenario: Color Scheme > Edit... opens the column's colour-coding dialog
    When user picks "Color Scheme > Edit..." from the context menu of pc plot viewer
    Then "Color-coding: AGE" dialog should be visible
    When user presses Escape
    Then "Color-coding: AGE" dialog should be absent
    And no errors should have been logged

  Scenario: Color Min and Color Max clamp the scale and recolour the lines
    When user sets properties of pc plot viewer:
      | Color Min | 30 |
      | Color Max | 60 |
    Then "Color Min" property of pc plot viewer should be "30"
    And the "color scale" area of pc plot viewer should have repainted
    And pc plot viewer should have repainted by at least 2000 pixels
    When user sets "Color Axis Type" property of pc plot viewer to "logarithmic"
    Then the "color scale" area of pc plot viewer should have repainted
    When user sets properties of pc plot viewer:
      | Color Axis Type | linear |
      | Color Min       |        |
      | Color Max       |        |
    Then pc plot viewer should have repainted by at least 2000 pixels
    When user sets "Color" property of pc plot viewer to ""
    Then pc plot viewer should not have a "color scale" area
    And no errors should have been logged

  Scenario: A conditional colour coding on the column hands the plot its bins
    When user sets "Color" property of pc plot viewer to "HEIGHT"
    Then pc plot viewer should have a "color scale" area
    When user colors "HEIGHT" column conditionally:
      | 20-150  | #00FF00 |
      | 150-250 | #FFA500 |
    Then "HEIGHT" column should be color-coded conditionally
    And legend of pc plot viewer should be visible
    And the legend of pc plot viewer should list 3 items
    And "20-150" legend item in legend of pc plot viewer should be visible
    And "150-250" legend item in legend of pc plot viewer should be visible
    And pc plot viewer should not have a "color scale" area
    When user removes the coloring of "HEIGHT" column
    Then "HEIGHT" column should have no color coding
    And legend of pc plot viewer should be hidden
    And pc plot viewer should have a "color scale" area
    When user sets "Color" property of pc plot viewer to ""
    Then no errors should have been logged

  Scenario: A DateTime colour column is split by its Color Map
    When user sets properties of pc plot viewer:
      | Color     | STARTED |
      | Color Map | year    |
    Then legend of pc plot viewer should be visible
    And the legend of pc plot viewer should list 3 items
    And "1990" legend item in legend of pc plot viewer should be visible
    And pc plot viewer should show 1000 rows
    When user sets "Color Map" property of pc plot viewer to "quarter"
    Then the legend of pc plot viewer should list 4 items
    And "Q1" legend item in legend of pc plot viewer should be visible
    And pc plot viewer should have repainted by at least 1000 pixels
    When user sets "Color Map" property of pc plot viewer to "month"
    Then the legend of pc plot viewer should list 12 items
    And "January" legend item in legend of pc plot viewer should be visible
    When user sets properties of pc plot viewer:
      | Color Map |  |
      | Color     |  |
    Then legend of pc plot viewer should be hidden
    And pc plot viewer should show 1000 rows
    And no errors should have been logged
