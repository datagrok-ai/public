@guide
Feature: SAR Matrix Summary without fragment columns
  The component cards rest on a decomposition the table already holds: only there does one value of a
  column mean the same thing in every row. Under the default fragmentation the molecules are cut
  automatically, a substituent label is local to the series it was cut from, and one pooled offset
  would average unrelated quantities — so no component card may appear, and the rest of the Summary
  must work exactly as it did before they existed.

  Demo: 9,999 ChEMBL compounds with hERG pIC50 on 3,109 of them, already a log scale.

  Scenario: no component cards under fragmentation
    Given user is logged in
    And simple mode is off
    And user opens sar-matrix-demo dataset
    Given user watches the task bar
    When user calls "Chem:SarMatrixAnalysis" function with:
      | table               | table              |
      | molecules           | column:smiles      |
      | activity            | column:hERG_pIC50  |
      | scaling             | none               |
      | activityDirection   | Higher is better   |
      | fragmentCutoff      | 0.4                |
      | fragmentationLevels | 3                  |
      | predictVirtual      | true               |
      | useMcsAnchors       | false              |
      | seriesColumn        |                    |
    Then the task bar should have finished "Building SAR matrices"
    When user clicks on "Summary" tab
    Then tab panel should contain text "Start here"
    And tab panel should contain text "113 L1 · 14 L2"
    And tab panel should contain text "open the tab to detect"
    When user clicks on "SAR Transfer" tab
    Then tab panel should contain text "transfer"
    When user clicks on "Summary" tab
    Then tab panel should not contain text "open the tab to detect"
    When user clicks on "Effects" summary segment
    Then tab panel should contain text "R-groups — within-series ranking"
    And tab panel should contain text "Cores are not comparable across series here"
    And tab panel should not contain text "offsets from the additive fit"
    And tab panel should not contain text "moves pIC50 most"
    When user clicks on "Method" summary segment
    Then tab panel should contain text "Fit quality"
    And tab panel should contain text "Two different R²"
    And tab panel should not contain text "Additive fit over"
    And no error or warning balloon should have been shown
