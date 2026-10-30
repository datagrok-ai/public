@journey @realizes:chem.cp.sar-matrix
Feature: SAR Matrix cut from the molecules of a small table
  The first 600 molecules of the SAR Matrix demo set, with hERG pIC50: small enough that the dialog's
  default path, which cuts the molecules into cores and substituents, builds in seconds. They make one
  series, a 2 × 5 matrix of 6 measured compounds and 4 predictions, none of them beating the series' best,
  and the Summary says so. Every scenario starts by putting its own tab on screen and puts back the
  settings it changed. Predictions turned off leave the Worth making card saying they are off. Rows
  deleted after a build rebuild the matrices over the rows that remain, so a measured cell still makes
  its own compound current.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the context panel is open
    And user opens sar-matrix-demo dataset keeping the first 600 rows

  Scenario: The dialog cuts the molecules and opens on the Summary
    When user picks "Chem > Analyze > SAR Matrix..." from the top menu
    Then "SAR Matrix" dialog should be visible
    When user selects "smiles" in Molecules input in "SAR Matrix" dialog
    And user selects "hERG_pIC50" in Activity input in "SAR Matrix" dialog
    And user selects "none" in Scaling input in "SAR Matrix" dialog
    And user selects "Higher is better" in Direction input in "SAR Matrix" dialog
    And user clicks on OK button in "SAR Matrix" dialog
    Then SAR Matrix Viewer viewer should have built its matrices
    And the "source" reading of SAR Matrix Viewer viewer should be "fragments"
    And the "matrices" reading of SAR Matrix Viewer viewer should be 1
    And the "measured" reading of SAR Matrix Viewer viewer should be 6
    And the "predicted" reading of SAR Matrix Viewer viewer should be 4
    And the "tab" reading of SAR Matrix Viewer viewer should be "Summary"
    And the "series" reading of SAR Matrix Viewer viewer should be the text of name of first key series row
    And Worth making card should contain text "No prediction resting on 2 or more measured compounds on both axes beats its series' best measured compound."
    And no errors should have been logged

  Scenario: Predictions turned off leave Worth making saying so
    When user clicks on Summary tab in SAR Matrix Viewer viewer
    And user sets "Predict virtual analogs" property of SAR Matrix Viewer viewer to "false"
    Then SAR Matrix Viewer viewer should have built its matrices
    And the "predicted" reading of SAR Matrix Viewer viewer should be 0
    And Worth making card should contain text "Predicted analogs are off — turn on"
    When user sets "Predict virtual analogs" property of SAR Matrix Viewer viewer to "true"
    Then SAR Matrix Viewer viewer should have built its matrices
    And the "predicted" reading of SAR Matrix Viewer viewer should be 4
    And Worth making card should not contain text "Predicted analogs are off"
    And no errors should have been logged

  Scenario: Rows deleted after a build rebuild the matrices, and a cell still opens its own compound
    When user clicks on Summary tab in SAR Matrix Viewer viewer
    And user selects the first 10 rows
    And user deletes the selected rows
    Then the table should have 590 rows
    And SAR Matrix Viewer viewer should have built its matrices
    And the "matrices" reading of SAR Matrix Viewer viewer should be 1
    And the "measured" reading of SAR Matrix Viewer viewer should be 6
    When user clicks on SAR Matrix tab in SAR Matrix Viewer viewer
    And user clicks on the "cell 1 of c0" area of grid viewer
    Then "CMPD_CHEMBLID" of the current row should be "CHEMBL1934525"
    And no errors should have been logged
