@journey
Feature: Curve JSON typed into a formula, a series without a name, and a broken curve cell
  Translated from TestTrack Curves/curves-data-robustness.md, scenarios 1 and 2, and
  curves-statistics.md, scenario 3. A formula that returns JSON adds its column and leaves Add New
  Column usable (public#3792, GROK-20232); a series without a name still gives its statistic column
  through the "+" of the Fit pane (GROK-17492); a click on a curve cell holding truncated JSON logs no
  error (GROK-15258, public#2748).

  Curve JSON with no fit function is detected as a curve column (GROK-18073) in the package's
  detector-tests ("curve JSON needs only series with points"); what a broken cell draws (the red cross)
  is the renderer's, tested by converter-tests (parseCellValue.malformedJson).

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
    When user presses Escape
    Then the "Add New Column" dialog should close

  Scenario: The statistic of a series without a name becomes a column
    Given user closes all views
    And user opens a table "unnamed-series" with:
      | id | curve                                                                                  |
      | 1  | {"series":[{"points":[{"x":0.1,"y":1},{"x":1,"y":5},{"x":10,"y":9}]}]}                 |
      | 2  | {"series":[{"points":[{"x":0.1,"y":2},{"x":1,"y":6},{"x":10,"y":8}]}]}                 |
    And the context panel is open
    When user clicks on the "cell 1 of curve" area of grid
    Then the current column should be "curve"
    Given "Fit" pane in context panel is expanded
    When user clicks on first "Calculate rSquared for the whole column" icon in context panel
    Then the table should have the columns "id, curve, curve 1 rSquared"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A click on a truncated curve cell logs no error
    Given user closes all views
    And user opens a table "broken-curves" with:
      | id | curve                                                                                  |
      | 1  | {"series":[{"points":[{"x":0.1,"y":1},{"x":1,"y":5},{"x":10,"y":9}]}]}                 |
      | 2  | {"series":[{"points":[{"x":0.1                                                         |
      | 3  | {"series":[{"points":[{"x":0.1,"y":2},{"x":1,"y":6},{"x":10,"y":8}]}]}                 |
    Then "curve" column should have semantic type "fit"
    When user clicks on the "cell 2 of curve" area of grid
    Then the current column should be "curve"
    And no errors should have been logged
    And no error or warning balloon should have been shown
