@guide @help:datagrok/solutions/domains/chem
Feature: SAR Matrix on a PROTAC set
  A guide: the answer to "my compounds already come in parts — can the SAR Matrix use them?". A
  PROTAC is a warhead, a linker and an E3 ligand, and a patent set usually carries the three as
  columns of their own. Given them, the analysis never cuts a molecule: the matrices are a grouping
  of the table, one component runs across the columns and the rest fold into the row.

  What that buys is on the Summary. One additive fit ranks all three parts at once and says which of
  them moves the endpoint furthest; the measured half pools matched pairs for each part separately,
  because a pair alike in everything but the linker is a grouping on the warhead and the ligand. No
  rebuild is needed to look at another component.

  Demo: 2,792 degraders from PROTAC-PatentDB, with the linker as the core and the warhead across the
  columns. The endpoint is a predicted solubility in log units — already a log, already
  higher-is-better, so the analysis must leave it untransformed.

  Scenario: warhead, linker and E3 ligand in one run
    Given user is logged in
    And simple mode is off
    And user opens protac-degraders dataset
    When user picks "Chem > Analyze > SAR Matrix..." from the top menu
    Then "SAR Matrix" dialog should be visible
    When user selects "Compound" in Molecules input in "SAR Matrix" dialog
    And user selects "Solubility logS (pred)" in Activity input in "SAR Matrix" dialog
    Then Scaling input in "SAR Matrix" dialog should be disabled
    And Scaling input in "SAR Matrix" dialog should have value "none"
    When user selects "Higher is better" in Direction input in "SAR Matrix" dialog
    And user checks "Use existing R-groups" input in "SAR Matrix" dialog
    Then Core input in "SAR Matrix" dialog should be visible
    When user selects "Linker" in Core input in "SAR Matrix" dialog
    And user clicks on editor of R-groups input in "SAR Matrix" dialog
    Then "Select columns..." dialog should be visible
    And the "text of cell 1 of __name" reading of grid viewer in "Select columns..." dialog should be "Warhead"
    And the "text of cell 3 of __name" reading of grid viewer in "Select columns..." dialog should be "E3 Ligand"
    When user clicks on None label in "Select columns..." dialog
    And user clicks on the "cell 1 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on the "cell 3 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    And user selects "Warhead" in "Matrix columns" input in "SAR Matrix" dialog
    Then "SAR Matrix" dialog should contain text "Columns: Warhead"
    When user clicks on OK button in "SAR Matrix" dialog
    Then SAR Matrix Viewer viewer should be visible
    And SAR Matrix Viewer viewer should have built its matrices
    And the "source" reading of SAR Matrix Viewer viewer should be "R-group columns"
    And "Summary" tab should be visible
    And tab panel should contain text "core: Linker"
    And tab panel should contain text "across the matrix columns: Warhead"
    And tab panel should contain text "folded into the row: E3 Ligand"
    And tab panel should contain text "What to change"
    And tab panel should contain text "Best measured swap"
    When user clicks on "Effects" summary segment
    Then tab panel should contain text "Warhead — offsets from the additive fit"
    When user clicks on "E3 Ligand" effects tab
    Then tab panel should contain text "E3 Ligand — offsets from the additive fit"
    And tab panel should contain text "The SAR Matrix columns enumerate Warhead"
    When user clicks on "Measured in series" effects tab
    Then tab panel should contain text "Swapping Warhead"
    And tab panel should contain text "Swapping Linker"
    And tab panel should contain text "Swapping E3 Ligand"
    When user clicks on "Worth making" summary segment
    Then tab panel should contain text "Nothing clears the trust gate"
    When user clicks on "Overview" summary segment
    And user clicks on "SAR Matrix" tab
    Then "SAR Matrix" tab should be visible
    And no error or warning balloon should have been shown
