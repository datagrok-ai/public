@realizes:chem.cp.sar-matrix
Feature: SAR Matrix series from the user's own Series column
  spgi-100 carries its R-group decomposition and a Series column of its own (64 Triazoles, 21 Pyrrolidines,
  4 Diazabicyclooctane, 4 Aminopiperidines, 7 blank). Picked in the dialog, the Series column groups the
  compounds instead of their cores, and each value becomes a matrix named after it when its compounds
  connect: the Triazoles do, 49 of them sharing one R2. Among the 21 Pyrrolidines only two share an R2
  (20 different R2 groups), so the largest block they form holds two compounds, under the three a matrix
  needs; the other two series have four compounds and four different R2 groups each.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And user opens spgi dataset

  Scenario: A Series column names one matrix per value
    When user picks "Chem > Analyze > SAR Matrix..." from the top menu
    Then "SAR Matrix" dialog should be visible
    When user selects "Average Mass" in Activity input in "SAR Matrix" dialog
    And user selects "none" in Scaling input in "SAR Matrix" dialog
    And user checks "Use existing R-groups" input in "SAR Matrix" dialog
    And user selects "Core" in Core input in "SAR Matrix" dialog
    And user clicks on editor of R-groups input in "SAR Matrix" dialog
    Then "Select columns..." dialog should be visible
    When user clicks on None label in "Select columns..." dialog
    And user clicks on the "cell 2 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on the "cell 3 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on the "cell 4 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on the "cell 5 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on the "cell 6 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    And user selects "Series" in Series input in "SAR Matrix" dialog
    Then "SAR Matrix" dialog should contain text "One matrix per Series"
    When user clicks on OK button in "SAR Matrix" dialog
    Then SAR Matrix Viewer viewer should have built its matrices
    And the "matrices" reading of SAR Matrix Viewer viewer should be 1
    And the "series" reading of SAR Matrix Viewer viewer should be "Triazoles"
    And no errors should have been logged
