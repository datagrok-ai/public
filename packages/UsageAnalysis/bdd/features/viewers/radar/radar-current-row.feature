@viewers @realizes:charts.viewer.radar
Feature: Radar highlights the current row's line with Show Current Row
  The Radar draws the current row's line again on top in its current-row color, green (#00FF00),
  while Selection > Show Current Row is on (the default); off, that line is drawn in the ordinary
  light-blue line color and no green is left on the chart. The claim is the color inside the
  viewer's `view` area; with Show Mouse Over Row off (so a grid hover draws nothing) a click on
  another grid row repaints the chart, so the line moves with the current row. Which row the green
  line belongs to has no reading yet (MISSING.md).
  Translated from the TestTrack case Charts/radar (step 10), whose old spec only read the property
  back.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens demog-1000 dataset
    When user clicks on "Add viewer" icon
    Then "Add Viewer" dialog should be visible
    When user clicks on first "Radar" card in "Add Viewer" dialog
    Then radar viewer should be visible
    And "Values" property of radar viewer should be "AGE, HEIGHT, WEIGHT"

  Scenario: Show Current Row draws the current row's line in green, and off it is gone
    When user sets "Show Mouse Over Row" property of radar viewer to "false"
    And user clicks on the "cell 5 of AGE" area of grid
    Then row 5 should be current
    And the "view" area of radar viewer should contain the color "#00FF00"
    When user takes a snapshot of radar viewer
    And user clicks on the "cell 9 of AGE" area of grid
    Then row 9 should be current
    And radar viewer should have repainted
    When user clicks on settings icon of radar viewer
    Given "Selection" category in context panel is expanded
    When user unchecks "Show Current Row" property in context panel
    Then "Show Current Row" property of radar viewer should be "false"
    And the "view" area of radar viewer should not contain the color "#00FF00"
    When user checks "Show Current Row" property in context panel
    Then the "view" area of radar viewer should contain the color "#00FF00"
    And no errors should have been logged
