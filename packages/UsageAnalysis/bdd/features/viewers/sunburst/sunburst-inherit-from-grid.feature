@viewers @realizes:charts.viewer.sunburst
Feature: Sunburst takes its segment colors from the grid's color coding with Inherit From Grid
  The Sunburst colors a segment by its column's categorical color coding while Color > Inherit From
  Grid is on (the default); off, it uses its own palette, where the first category, F, is #1F77B4.
  SEX on demog-1000 is colored F red and M blue through the column's color coding, as the grid's
  Color Coding > Edit dialog would set it; the hierarchy is SEX alone, picked in the Select columns
  dialog of the Hierarchy property. A segment's color is read inside its `segment <path>` area.
  Translated from the TestTrack case Charts/sunburst (step 3.3, Inherit From Grid), whose old spec
  only read the property back.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens demog-1000 dataset
    When user clicks on "Add viewer" icon
    Then "Add Viewer" dialog should be visible
    When user clicks on first "Sunburst" card in "Add Viewer" dialog
    Then sunburst viewer should be visible
    When user picks "Properties..." from the context menu of the "view" area of sunburst viewer
    And user clicks on "..." button in "Hierarchy" property in context panel
    Then "Select columns..." dialog should be visible
    When user clicks on "None" link in "Select columns..." dialog
    Then "0 checked" text in "Select columns..." dialog should be visible
    When user toggles the "SEX" column in the column list of "Select columns..." dialog
    Then "1 checked" text in "Select columns..." dialog should be visible
    When user clicks on OK button in "Select columns..." dialog
    Then "Select columns..." dialog should be absent
    And the "hierarchy columns" reading of sunburst viewer should be "SEX"
    And the "segments" reading of sunburst viewer should be 2

  Scenario: Segments follow the grid's categorical colors, and their own palette with the option off
    Then the "segment F" area of sunburst viewer should contain the color "#1F77B4"
    And the "segment F" area of sunburst viewer should not contain the color "#FF0000"
    When user colors "SEX" column categorically:
      | F | #FF0000 |
      | M | #0000FF |
    Then "SEX" column should be color-coded categorically
    And the "segment F" area of sunburst viewer should contain the color "#FF0000"
    And the "segment M" area of sunburst viewer should contain the color "#0000FF"
    And the "segment F" area of sunburst viewer should not contain the color "#1F77B4"
    Given "Color" category in context panel is expanded
    When user unchecks "Inherit From Grid" property in context panel
    Then the "segment F" area of sunburst viewer should contain the color "#1F77B4"
    And the "segment F" area of sunburst viewer should not contain the color "#FF0000"
    And the "segment M" area of sunburst viewer should not contain the color "#0000FF"
    When user checks "Inherit From Grid" property in context panel
    Then the "segment F" area of sunburst viewer should contain the color "#FF0000"
    And the "segment M" area of sunburst viewer should contain the color "#0000FF"
    And no errors should have been logged
