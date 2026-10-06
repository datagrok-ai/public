@guide @help:datagrok/solutions/domains/chem
Feature: SAR Matrix
  A guide: the answer to "how do I run a SAR Matrix analysis?". Chem | Analyze | SAR Matrix groups
  a compound set into analog series and lays each one out as cores against the substituents they
  share, coloring every cell by potency and filling the combinations nobody has made with a
  Free-Wilson prediction. The launch dialog plots the activity column beside the inputs, so the
  scaling can be read off the distribution before the analysis runs. Demo: the set the SAR Matrix
  demo opens, read on hERG_pIC50 — already a pIC50, so the scaling stays none and the direction is
  declared higher-is-better, which is what tells the analysis to read the column as a log.

  The analysis opens on its Summary: what it concluded, before the grids it concluded it from. The
  walkthrough goes through that tab segment by segment — the scale everything rests on, the coverage
  of the run, which component to change and what was measured, the per-component rankings, the
  analogs worth making, and the gate they had to clear — and only then into the matrices, the
  transfers between them and the list of what to make.

  Scenario: SAR matrix demo
    Given user is logged in
    And simple mode is off
    And user opens sar-matrix-demo dataset
    When user picks "Chem > Analyze > SAR Matrix..." from the top menu
    Then "SAR Matrix" dialog should be visible
    And Activity input in "SAR Matrix" dialog should contain text "CYP3A4"
    When user selects "hERG_pIC50" in Activity input in "SAR Matrix" dialog
    And user selects "none" in Scaling input in "SAR Matrix" dialog
    And user selects "Higher is better" in Direction input in "SAR Matrix" dialog
    Given user watches the task bar
    When user clicks on OK button in "SAR Matrix" dialog
    Then the task bar should have finished "Building SAR matrices"
    And "Summary" tab should be visible
    Then tab panel should contain text "untransformed"
    And tab panel should contain text "Cores and substituents found by cutting the molecules"
    And tab panel should contain text "Start here"
    When user clicks on "L1" tier chip
    Then tab panel should contain text "Start here"
    When user clicks on "All" tier chip
    And user clicks on "Effects" summary segment
    Then tab panel should contain text "within-series ranking"
    And tab panel should contain text "Cores are not comparable across series here"
    And tab panel should contain text "What buys potency inside these series"
    When user clicks on "Worth making" summary segment
    Then tab panel should contain text "Already in the table, never assayed"
    And tab panel should contain text "Worth making"
    When user clicks on "Method" summary segment
    Then tab panel should contain text "Fit quality"
    And tab panel should contain text "Trust gate"
    When user clicks on "Overview" summary segment
    And user clicks on "SAR Matrix" tab
    Then "SAR Matrix" tab should be visible
    When user clicks on "SAR Transfer" tab
    Then tab panel should not contain text "Detecting SAR transfers"
    When user clicks on "Make list" tab
    And user clicks on "SAR Matrix" tab
    Then no error or warning balloon should have been shown
