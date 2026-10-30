@realizes:chem.cp.sar-matrix
Feature: SAR Matrix offers the compounds the table holds but never assayed
  Sixteen para-disubstituted benzenes, four substituents on each side, with pIC50 on all but one:
  PARA-10, 4-chlorophenol, is in the table with a blank activity. Cut from the molecules they make one
  series, a 4 × 4 matrix of 15 measured compounds whose only hole is that compound, so nothing is left
  to make. Its row and its column each hold three measured compounds, so the model predicts it, and the
  Summary lists it above Worth making as already in the table, never assayed. Queueing it puts it on
  the Make list as untested, not as an analog to synthesize.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens sar-matrix-para dataset

  Scenario: A compound never assayed is offered for testing and queued on the Make list
    When user picks "Chem > Analyze > SAR Matrix..." from the top menu
    Then "SAR Matrix" dialog should be visible
    When user selects "smiles" in Molecules input in "SAR Matrix" dialog
    And user selects "pIC50" in Activity input in "SAR Matrix" dialog
    And user selects "none" in Scaling input in "SAR Matrix" dialog
    And user selects "Higher is better" in Direction input in "SAR Matrix" dialog
    And user clicks on OK button in "SAR Matrix" dialog
    Then SAR Matrix Viewer viewer should have built its matrices
    And the "matrices" reading of SAR Matrix Viewer viewer should be 1
    And the "rows" reading of SAR Matrix Viewer viewer should be 4
    And the "columns" reading of SAR Matrix Viewer viewer should be 4
    And the "measured" reading of SAR Matrix Viewer viewer should be 15
    And the "predicted" reading of SAR Matrix Viewer viewer should be 0
    And Never assayed card should contain text "Already in the table, never assayed · 1"
    And first never assayed row should contain text "n=3"
    And Queue for testing button should have text "Queue this compound for testing"
    When user clicks on Queue for testing button
    Then an info balloon containing "Added 1 to the Make list (1 total)" should have been shown
    And the "make list" reading of SAR Matrix Viewer viewer should be 1
    When user clicks on Make list tab
    Then Make list title should have text "1 compound"
    And the "text of cell 1 of Structure" reading of Make list grid viewer should be "Oc1ccc(Cl)cc1"
    And the "text of cell 1 of Status" reading of Make list grid viewer should be "Untested"
    And no errors should have been logged
