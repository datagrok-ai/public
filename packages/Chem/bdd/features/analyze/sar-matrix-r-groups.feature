@journey @realizes:chem.cp.sar-matrix
Feature: SAR Matrix from the R-group columns SPGI already has
  SPGI carries its own R-group decomposition: a Core column and R1, R2, R3, R100 and R101, with
  attachment points written [*:n] (R100 and R101 both hold one bridge). Its first 1,500 rows hold every
  one of these columns, the bridge on 25 compounds, and build the same 6 series and 5 components as the
  whole set, in a fifth of the time. With Use existing R-groups on, the SAR Matrix dialog says what the
  set-up still lacks and keeps OK disabled until it has a core and an R-group, offers only those columns
  as the core and the R-groups, drops the core when it is also picked as an R-group, puts R2 on the
  matrix columns, and says what the rows and columns will be. OK builds one series per core from those
  columns, with Average Mass as the value. The Summary ranks the bridge once: R100 and R101 hold it whole
  in every compound, so they are one component.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens SPGI-full dataset keeping the first 1500 rows

  Scenario: The dialog says what an R-group set-up lacks and states the layout
    When user picks "Chem > Analyze > SAR Matrix..." from the top menu
    Then "SAR Matrix" dialog should be visible
    And Molecules input in "SAR Matrix" dialog should contain text "Structure"
    When user checks "Use existing R-groups" input in "SAR Matrix" dialog
    Then "SAR Matrix" dialog should contain text "Pick the core column."
    And OK button in "SAR Matrix" dialog should be disabled
    When user selects "Core" in Core input in "SAR Matrix" dialog
    Then "SAR Matrix" dialog should contain text "Pick the R-groups."
    And OK button in "SAR Matrix" dialog should be disabled
    When user clicks on editor of R-groups input in "SAR Matrix" dialog
    Then "Select columns..." dialog should be visible
    And the "text of cell 1 of __name" reading of grid viewer in "Select columns..." dialog should be "Core"
    And the "text of cell 2 of __name" reading of grid viewer in "Select columns..." dialog should be "R1"
    And the "text of cell 3 of __name" reading of grid viewer in "Select columns..." dialog should be "R2"
    And the "text of cell 6 of __name" reading of grid viewer in "Select columns..." dialog should be "R101"
    When user clicks on None label in "Select columns..." dialog
    And user clicks on the "cell 1 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on the "cell 2 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on the "cell 3 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on the "cell 4 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on the "cell 5 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on the "cell 6 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then an info balloon containing "Core is the core column, so it was dropped" should have been shown
    And "Matrix columns" input in "SAR Matrix" dialog should have value "R2"
    And "Matrix columns" input in "SAR Matrix" dialog should offer "R1, R2, R3, R100, R101"
    And "SAR Matrix" dialog should contain text "Rows: Core + R1 + R3 + R100 + R101"
    And "SAR Matrix" dialog should contain text "Columns: R2"
    And OK button in "SAR Matrix" dialog should be enabled
    When user selects "Average Mass" in Activity input in "SAR Matrix" dialog
    And user selects "none" in Scaling input in "SAR Matrix" dialog
    Then no errors should have been logged

  Scenario: OK builds one series per core from those columns
    When user clicks on OK button in "SAR Matrix" dialog
    Then SAR Matrix Viewer viewer should be visible
    And SAR Matrix Viewer viewer should have built its matrices
    And the "source" reading of SAR Matrix Viewer viewer should be "R-group columns"
    And the "matrices" reading of SAR Matrix Viewer viewer should be 6
    And What to change title should have text "What to change · 5 components"
    When user clicks on more components fold
    Then R100 + R101 component row should be visible
    And no errors should have been logged
