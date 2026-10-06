@guide
Feature: SAR Matrix Summary over four R-group columns
  A decomposition is not limited to a core and two R-groups. This set carries a scaffold column and
  four R-group columns — five components in all — and was built from known additive coefficients, so
  what the Summary reports can be checked against a truth rather than against itself: R2 should move
  pIC50 most, then the Scaffold, then R1, then R4, and R3 should barely move it.

  Five components stacked would make the Effects segment five cards tall before the per-series evidence
  below them, so the segment carries its own tabs — one per component, then the per-series card — and
  each tab shows its component's spread, so which component matters most is readable without opening
  any of them.

  Scenario: five components, one tab each
    Given user is logged in
    And simple mode is off
    And user opens four-rgroups dataset
    Given user watches the task bar
    When user calls "Chem:SarMatrixAnalysis" function with:
      | table               | table                 |
      | molecules           | column:Compound       |
      | activity            | column:pIC50          |
      | scaling             | none                  |
      | activityDirection   | Higher is better      |
      | fragmentCutoff      | 0.4                   |
      | fragmentationLevels | 3                     |
      | predictVirtual      | true                  |
      | useMcsAnchors       | false                 |
      | seriesColumn        |                       |
      | coreColumn          | column:Scaffold       |
      | fragmentColumns     | columns:R1,R2,R3,R4   |
      | columnAxis          | R2                    |
    Then the task bar should have finished "Building SAR matrices"
    When user clicks on "Summary" tab
    And user clicks on "Effects" summary segment
    Then tab panel should contain text "Changing R2 moves pIC50 most: R2 > Scaffold > R1 > R4"
    And tab panel should contain text "R2 — offsets from the additive fit"
    When user clicks on "Scaffold" effects tab
    Then tab panel should contain text "Scaffold — offsets from the additive fit"
    When user clicks on "R3" effects tab
    Then tab panel should contain text "R3 — offsets from the additive fit"
    When user clicks on "Measured in series" effects tab
    Then tab panel should contain text "R2 — within-series ranking"
    And no error or warning balloon should have been shown
