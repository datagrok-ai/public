@guide @help:datagrok/solutions/domains/chem
Feature: SAR Matrix tab
  A guide: the answer to "how do I look at one series as a matrix, and at the compound behind a cell?". The
  SAR Matrix tab lists the series on the left and draws the selected one as a matrix: cores down the rows,
  substituents across the columns, a solid cell a compound of the table, a dashed one a predicted analog. A
  click on a card opens that series; a click on a cell opens its compound in the Context Panel. Filmed for
  the help page's SAR Matrix tab. Demo: the SAR Matrix demo project.
  Filmed in a 1000x625 window (BDD_GUIDE_VIEWPORT=1000x625): at 800x500 neither the Summary nor the
  matrix beside the Context Panel fits.

  Scenario: Open a series and one of its cells
    Given user is logged in
    And simple mode is off
    And the package autostarts have completed
    And the context panel is open
    And user opens the "SarMatrixDemo" project
    Then SAR Matrix Viewer viewer should have built its matrices
    When user clicks on SAR Matrix tab in SAR Matrix Viewer viewer
    Then Series 69 series card should be visible
    # caption: Click a card to open that series
    When user clicks on Series 1 series card
    Then the "series" reading of SAR Matrix Viewer viewer should be "Series 1"
    # caption: Click a cell to open its compound in the Context Panel
    When user clicks on the "cell 1 of c0" area of grid viewer
    Then context panel should contain text "SAR analysis"
    And no errors should have been logged
