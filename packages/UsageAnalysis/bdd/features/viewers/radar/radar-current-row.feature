@journey @viewers @realizes:charts.viewer.radar
Feature: Radar: Show Current Row, and Pick Up / Apply between two Radars
  The Radar draws the current row's line again on top in its current-row color, green (#00FF00),
  while Selection > Show Current Row is on (the default); off, that line is drawn in the ordinary
  light-blue line color and no green is left on the chart. The claim is the color inside the
  viewer's `view` area; with Show Mouse Over Row off (so a grid hover draws nothing) a click on
  another grid row repaints the chart, so the line moves with the current row. Which row the green
  line belongs to has no reading yet (MISSING.md).
  Pick Up on one Radar and Apply on a second copies a setting, as it does for the core viewers,
  while a later change of the first leaves the second as it was (public#1472).
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

  @realizes:public#1472
  Scenario: Pick Up / Apply copies a setting between two Radar viewers, and the two stay apart (public#1472)
    Given user adds a radar viewer
    Then the open tableview should have 2 radar viewers
    When user sets "normalization" property of first radar viewer to "Global"
    Then "normalization" property of second radar viewer should be "Column"
    When user picks "Pick Up / Apply > Pick Up" from the viewer menu of first radar viewer
    And user picks "Pick Up / Apply > Apply" from the viewer menu of second radar viewer
    Then "normalization" property of second radar viewer should be "Global"
    When user takes a snapshot of first radar viewer
    And user sets "normalization" property of first radar viewer to "Column"
    Then first radar viewer should have repainted
    And "normalization" property of second radar viewer should be "Global"
    And no errors should have been logged
