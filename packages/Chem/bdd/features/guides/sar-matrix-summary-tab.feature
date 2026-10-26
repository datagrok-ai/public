@guide @help:datagrok/solutions/domains/chem
Feature: SAR Matrix Summary tab
  A guide: the answer to "what does the SAR Matrix Summary tell me, and what does each part open?". A click
  on a key series opens its best core, on an answer the evidence behind it, and on an analog worth making
  what its prediction adds up from. Filmed for the help page's Summary tab. Demo: the SAR Matrix demo
  project.
  Filmed in a 1000x625 window (BDD_GUIDE_VIEWPORT=1000x625): at 800x500 neither the Summary nor the
  matrix beside the Context Panel fits.

  Scenario: Read the Summary
    Given user is logged in
    And simple mode is off
    And the package autostarts have completed
    And user opens the "SarMatrixDemo" project
    And SAR Matrix Viewer viewer should have built its matrices
    # caption: Open a series to see its best core
    When user clicks on first key series row
    Then key series detail should be visible
    # caption: Scroll down to the answers
    When user scrolls What improves potency answer tile into view
    # caption: Open an answer to see the evidence behind it
    And user clicks on What improves potency answer tile
    Then Summary evidence should contain text "What buys potency"
    And first evidence row should contain text "weakest pair"
    # caption: Open an analog to see what its prediction adds up from
    When user clicks on first analog worth making row
    Then analog worth making detail should be visible
    And no errors should have been logged
