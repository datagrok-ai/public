@guide @help:datagrok/solutions/domains/chem
Feature: SAR Matrix on a PROTAC set
  A guide: the answer to "my compounds already come in parts, can the SAR Matrix use them?". A PROTAC is
  a warhead, a linker and an E3 ligand, and a patent set usually carries the three as columns of their
  own. Given them, the analysis never cuts a molecule, and one additive fit on the Summary ranks all
  three parts at once. Filmed for the help page's walkthrough. Demo: 2,792 degraders from
  PROTAC-PatentDB, with the linker as the core and the warhead across the columns; the endpoint is a
  predicted solubility in log units, so the analysis leaves it untransformed.
  Filmed in a 1000x625 window (BDD_GUIDE_VIEWPORT=1000x625): at 800x500 neither the Summary nor the
  matrix beside the Context Panel fits.

  Scenario: Warhead, linker and E3 ligand in one run
    Given user is logged in
    And simple mode is off
    And the package autostarts have completed
    And user opens protac-degraders dataset
    When user picks "Chem > Analyze > SAR Matrix..." from the top menu
    Then "SAR Matrix" dialog should be visible
    When user selects "Compound" in Molecules input in "SAR Matrix" dialog
    And user selects "Solubility logS (pred)" in Activity input in "SAR Matrix" dialog
    # caption: logS is negative, so it cannot be log-scaled
    Then Scaling input in "SAR Matrix" dialog should be disabled
    And Scaling input in "SAR Matrix" dialog should have value "none"
    When user selects "Higher is better" in Direction input in "SAR Matrix" dialog
    And user checks "Use existing R-groups" input in "SAR Matrix" dialog
    Then Core input in "SAR Matrix" dialog should be visible
    When user selects "Linker" in Core input in "SAR Matrix" dialog
    And user clicks on editor of R-groups input in "SAR Matrix" dialog
    Then "Select columns..." dialog should be visible
    When user clicks on None label in "Select columns..." dialog
    And user clicks on the "cell 1 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on the "cell 3 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    And user selects "Warhead" in "Matrix columns" input in "SAR Matrix" dialog
    And user clicks on OK button in "SAR Matrix" dialog
    Then SAR Matrix Viewer viewer should have built its matrices
    # caption: Scroll down to What to change
    When user scrolls E3 Ligand component row into view
    # caption: What to change ranks the warhead, the linker and the E3 ligand in one fit
    Then E3 Ligand component row should be visible
    And What to change title should have text "What to change · 3 components"
    And no errors should have been logged
