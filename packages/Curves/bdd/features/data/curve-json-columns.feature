@journey
Feature: Curve JSON typed into a formula, a series without a name, and a broken curve cell
  Translated from TestTrack Curves/curves-data-robustness.md, scenarios 1 and 2, and
  curves-statistics.md, scenario 3. A formula that returns JSON adds its column and leaves Add New
  Column usable (public#3792, GROK-20232); curve JSON needs only series and points to be detected as
  a curve column (GROK-18073); a series without a name still gives its statistic column (GROK-17492);
  a curve cell holding truncated JSON leaves the rest of the table drawn (GROK-15258, public#2748).

  Not claimed: the red cross a broken cell draws — the renderer reports no readings.

  Background:
    Given user is logged in
    And simple mode is off
    And the package autostarts have completed
    And the "Curves" package is installed

  Scenario: A JSON formula adds its column and Add New Column opens again
    Given user opens demog-1000 dataset
    When user picks "Edit > Add New Column..." from the top menu
    Then "Add New Column" dialog should be visible
    When user pastes "{}" into code editor in "Add New Column" dialog
    And user types "jsoncol" into column name input
    And user clicks on OK button in "Add New Column" dialog
    Then the "Add New Column" dialog should close
    And the table should have a column "jsoncol"
    And no errors should have been logged
    When user picks "Edit > Add New Column..." from the top menu
    Then "Add New Column" dialog should be visible

  Scenario: Curve JSON without a fit function becomes a curve column
    When user pastes "'{\"series\":[{\"points\":[{\"x\":0.1,\"y\":1},{\"x\":1,\"y\":5},{\"x\":10,\"y\":9}]}]}'" into code editor in "Add New Column" dialog
    And user types "curvecol" into column name input
    And user clicks on OK button in "Add New Column" dialog
    Then the "Add New Column" dialog should close
    And the table should have a column "curvecol"
    And "curvecol" column should have semantic type "fit"
    And no errors should have been logged

  Scenario: The statistic of a series without a name becomes a column
    Given the context panel is open
    When user scrolls the grid to the "curvecol" column
    And user clicks on the "cell 1 of curvecol" area of grid
    Then the current column should be "curvecol"
    Given "Fit" pane in context panel is expanded
    When user clicks on first "Calculate rSquared for the whole column" icon in context panel
    Then a new column matching "^curvecol .*rSquared$" should have been added
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A truncated curve cell leaves the other rows drawn
    Given user closes all views
    And user opens a table "broken-curves" with:
      | id | curve                                                                                  |
      | 1  | {"series":[{"points":[{"x":0.1,"y":1},{"x":1,"y":5},{"x":10,"y":9}]}]}                 |
      | 2  | {"series":[{"points":[{"x":0.1                                                         |
      | 3  | {"series":[{"points":[{"x":0.1,"y":2},{"x":1,"y":6},{"x":10,"y":8}]}]}                 |
    Then "curve" column should have semantic type "fit"
    And the table should have 3 rows
    When user clicks on the "cell 2 of curve" area of grid
    Then the current column should be "curve"
    And no errors should have been logged
    And no error or warning balloon should have been shown
