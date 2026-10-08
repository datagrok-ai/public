@viewers @realizes:charts.viewer.sankey @realizes:charts.viewer.chord
Feature: Sankey and Chord columns, and redrawing on every filter change
  Sankey draws flows from Source values to Target values weighted by Value; Chord draws links
  between From and To categories. Both draw SVG and report what they drew: the Sankey its `node
  <name>` and `link <source> -> <target>` areas and its `nodes`, `node names` and `links` readings,
  the Chord its `category <name>` areas and its `categories` and `chords` readings. The viewers are
  added from the ribbon's Add viewer gallery on demog (5850 rows; SEX F 3243) and set up in the
  Context Panel; the filters are set on the Filter Panel's cards. Which columns Source and Target
  offer is read from their column pickers. A filter that no row passes leaves the Sankey empty and
  logs nothing (GROK-21110), and the Chord takes From set to the column To holds (GROK-21111).
  Translated from the TestTrack case Charts/charts-flow-viewers.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens demog dataset

  Scenario: Sankey starts on SEX, RACE and AGE and takes other columns (GROK-18048)
    When user clicks on "Add viewer" icon
    Then "Add Viewer" dialog should be visible
    When user clicks on first "Sankey" card in "Add Viewer" dialog
    Then sankey viewer should be visible
    And sankey viewer should be bound to table "demog"
    And the "node names" reading of sankey viewer should be "F, M, Caucasian, Other, Asian, Black"
    And the "links" reading of sankey viewer should be 5850
    When user clicks on grid
    And user clicks on settings icon of sankey viewer
    Then "Source" property in context panel should be visible
    And properties of sankey viewer should be:
      | Source | SEX  |
      | Target | RACE |
      | Value  | AGE  |
    And "Source" property in context panel should offer the columns "USUBJID, SEX, RACE, DIS_POP, DEMOG, SEVERITY"
    And "Target" property in context panel should offer the columns "USUBJID, SEX, RACE, DIS_POP, DEMOG, SEVERITY"
    And no errors should have been logged
    When user selects "RACE" in "Source" property in context panel
    Then "Source" property of sankey viewer should be "RACE"
    When user selects "DIS_POP" in "Target" property in context panel
    Then "Target" property of sankey viewer should be "DIS_POP"
    When user selects "WEIGHT" in "Value" property in context panel
    Then properties of sankey viewer should be:
      | Source | RACE    |
      | Target | DIS_POP |
      | Value  | WEIGHT  |
    And the "node names" reading of sankey viewer should contain "Psoriasis"
    And the "node names" reading of sankey viewer should contain "Caucasian"
    And the "node names" reading of sankey viewer should not contain "M"
    And no errors should have been logged

  Scenario: Sankey follows the filter (GROK-18035)
    When user clicks on "Add viewer" icon
    And user clicks on first "Sankey" card in "Add Viewer" dialog
    Then the "node names" reading of sankey viewer should contain "M"
    When user clicks on filter icon in toolbar
    Then filter panel should be visible
    When user clicks on the "category F of SEX" area of filter panel
    Then 3243 rows should pass the filter
    And the filter should pass exactly the rows where "SEX" is "F"
    And the "links" reading of sankey viewer should be 3243
    And the "node names" reading of sankey viewer should contain "F"
    And the "node names" reading of sankey viewer should not contain "M"
    And no errors should have been logged
    When user hovers over the "link F -> Caucasian" area of sankey viewer
    Then tooltip should contain text "2823 rows"
    When user hovers over the "link F -> Asian" area of sankey viewer
    Then tooltip should contain text "37 rows"
    When user moves the pointer away from sankey viewer
    Then no errors should have been logged
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then 5850 rows should pass the filter
    And the "links" reading of sankey viewer should be 5850
    And the "node names" reading of sankey viewer should contain "M"
    And no errors should have been logged

  Scenario: Sankey with a filter no row passes draws nothing and logs no error (GROK-21110)
    When user clicks on "Add viewer" icon
    And user clicks on first "Sankey" card in "Add Viewer" dialog
    Then the "links" reading of sankey viewer should be 5850
    When user clicks on filter icon in toolbar
    Then filter panel should be visible
    When user clicks on the "category true of CONTROL" area of filter panel
    And user clicks on the "category Asian of RACE" area of filter panel
    Then 0 rows should pass the filter
    And the "links" reading of sankey viewer should be 0
    And the "nodes" reading of sankey viewer should be 0
    And no errors should have been logged

  Scenario: Chord takes From set to the column To holds, and redraws on a filter change without a click (GROK-21111, GROK-17772)
    When user clicks on "Add viewer" icon
    Then "Add Viewer" dialog should be visible
    When user clicks on first "Chord" card in "Add Viewer" dialog
    Then chord viewer should be visible
    And the "categories" reading of chord viewer should be 6
    When user clicks on grid
    And user clicks on settings icon of chord viewer
    Then "From" property in context panel should be visible
    And properties of chord viewer should be:
      | From | SEX  |
      | To   | RACE |
    And no errors should have been logged
    When user selects "RACE" in "From" property in context panel
    Then "From" property of chord viewer should be "RACE"
    And the "categories" reading of chord viewer should be 4
    And no errors should have been logged
    When user selects "DIS_POP" in "To" property in context panel
    Then properties of chord viewer should be:
      | From | RACE    |
      | To   | DIS_POP |
    And chord viewer should have a "category Black" area
    And chord viewer should have a "category PsA" area
    And chord viewer should not have a "category F" area
    And no errors should have been logged
    When user clicks on filter icon in toolbar
    Then filter panel should be visible
    When user clicks on the "category Asian of RACE" area of filter panel
    Then the filter should pass exactly the rows where "RACE" is "Asian"
    And chord viewer should have a "category Asian" area
    And chord viewer should not have a "category Black" area
    And no errors should have been logged
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then 5850 rows should pass the filter
    And chord viewer should have a "category Black" area
    And no errors should have been logged
