@guide @help:datagrok/solutions/domains/chem
Feature: SAR Matrix on a compound set
  A guide: the answer to "how do I run a SAR Matrix analysis?". The analysis cuts the molecules into cores
  and substituents, groups them into series and opens on the Summary. Filmed for the help page's How to use
  steps. Demo: the SAR Matrix demo set, 10,000 molecules with hERG pIC50.
  Filmed in a 1000x625 window (BDD_GUIDE_VIEWPORT=1000x625): at 800x500 neither the Summary nor the
  matrix beside the Context Panel fits.

  Scenario: Run a SAR Matrix analysis
    Given user is logged in
    And simple mode is off
    And the package autostarts have completed
    And user opens sar-matrix-demo dataset
    When user picks "Chem > Analyze > SAR Matrix..." from the top menu
    Then "SAR Matrix" dialog should be visible
    When user selects "smiles" in Molecules input in "SAR Matrix" dialog
    And user selects "hERG_pIC50" in Activity input in "SAR Matrix" dialog
    # caption: hERG pIC50 is already a log value, so no scaling
    And user selects "none" in Scaling input in "SAR Matrix" dialog
    Then Scaling input in "SAR Matrix" dialog should have value "none"
    When user selects "Higher is better" in Direction input in "SAR Matrix" dialog
    And user clicks on OK button in "SAR Matrix" dialog
    Then SAR Matrix Viewer viewer should have built its matrices
    # caption: The analysis opens on the Summary
    And Key series card should be visible
    And Worth making card should be visible
    And no errors should have been logged
