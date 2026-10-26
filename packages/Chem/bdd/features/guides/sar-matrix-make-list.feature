@guide @help:datagrok/solutions/domains/chem
Feature: SAR Matrix Make list
  A guide: the answer to "how do I collect the analogs I want to make?". The cart on an analog of the
  Summary's Worth making list adds it to the Make list tab, with its potency, the support behind it and the
  series, core and substituent it came from. Filmed for the help page's Make list tab. Demo: the SAR Matrix
  demo project.
  Filmed in a 1000x625 window (BDD_GUIDE_VIEWPORT=1000x625): at 800x500 neither the Summary nor the
  matrix beside the Context Panel fits.

  Scenario: Collect an analog into the Make list
    Given user is logged in
    And simple mode is off
    And the package autostarts have completed
    And user opens the "SarMatrixDemo" project
    Then SAR Matrix Viewer viewer should have built its matrices
    # caption: Add an analog worth making with its cart
    When user clicks on cart of first analog worth making row
    # caption: Open the Make list
    And user clicks on Make list tab
    # caption: The analog, with its potency and where it came from
    Then Make list title should have text "1 compound"
    And the "text of cell 1 of Status" reading of Make list grid viewer should be "Virtual"
    And no errors should have been logged
