@guide @help:datagrok/solutions/domains/chem
Feature: SAR Transfer tab
  A guide: the answer to "which of my cores could borrow what another core has learned?". The SAR Transfer
  tab pairs cores whose potency trends run in parallel across the substituents both have tried. A source
  that transfers to several targets is one card, and the Transfer to list above the grid switches between
  its targets while the source stays. Filmed for the help page's SAR Transfer tab. Demo: the SAR Matrix
  demo project, where one core of Series 69 transfers to two others.
  Filmed in a 1000x625 window (BDD_GUIDE_VIEWPORT=1000x625): at 800x500 neither the Summary nor the
  matrix beside the Context Panel fits.

  Scenario: Switch a transfer to its other target
    Given user is logged in
    And simple mode is off
    And the package autostarts have completed
    And user opens the "SarMatrixDemo" project
    Then SAR Matrix Viewer viewer should have built its matrices
    # caption: Open the SAR Transfer tab: each card is a source core
    When user clicks on SAR Transfer tab
    Then "Transfer to" input should have value "Series 33 · Core 1 · R1"
    # caption: Switch the source to its other target
    When user selects "Series 98 · Core 1 · R1" in "Transfer to" input
    Then transfer title should have text "Series 69 · Core 3 → Series 98 · Core 1"
    And no errors should have been logged
