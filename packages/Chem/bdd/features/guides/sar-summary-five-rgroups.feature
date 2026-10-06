@guide
Feature: SAR Matrix Summary over five R-group columns
  Six components on 4,172 compounds, built from known additive coefficients so the fit can be checked
  against a truth: R2 should move pIC50 most, then the Scaffold, then R1, R5, R4, and R3 barely at all.
  One of the four scaffolds carries no R3, so that column also has a blank level, which is a level like
  any other rather than a missing value.

  Scenario: six components on four thousand compounds
    Given user is logged in
    And simple mode is off
    And user opens five-rgroups dataset
    Given user watches the task bar
    When user calls "Chem:SarMatrixAnalysis" function with:
      | table               | table                   |
      | molecules           | column:Compound         |
      | activity            | column:pIC50            |
      | scaling             | none                    |
      | activityDirection   | Higher is better        |
      | fragmentCutoff      | 0.4                     |
      | fragmentationLevels | 3                       |
      | predictVirtual      | true                    |
      | useMcsAnchors       | false                   |
      | seriesColumn        |                         |
      | coreColumn          | column:Scaffold         |
      | fragmentColumns     | columns:R1,R2,R3,R4,R5  |
      | columnAxis          | R2                      |
    Then the task bar should have finished "Building SAR matrices"
    When user clicks on "Summary" tab
    And user clicks on "Effects" summary segment
    Then tab panel should contain text "Changing R2 moves pIC50 most: R2 > Scaffold > R1 > R5"
    And tab panel should contain text "R2 — offsets from the additive fit"
    When user clicks on "R5" effects tab
    Then tab panel should contain text "R5 — offsets from the additive fit"
    When user clicks on "R3" effects tab
    Then tab panel should contain text "nothing at this position"
    When user clicks on "Method" summary segment
    Then tab panel should contain text "Fit quality"
    And tab panel should contain text "Where it holds best"
    And no error or warning balloon should have been shown
