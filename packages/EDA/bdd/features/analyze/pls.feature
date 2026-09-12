@journey @eda @realizes:ml.menu.analyze.pls
Feature: Partial least squares regression
  ML | Analyze | PLS... over cars, predicting price with three components. Translated from
  files/TestTrack/EDA/pls.md and the package's playwright/pls.test.ts, which stopped at the dialog.

  The case says to select every column in Using. The picker's All offers the numeric columns, price
  among them, and a column cannot be both predicted and a predictor: RUN stays disabled until price
  is unchecked, so the feature takes price back out. The picker lists the columns in the table's
  order, which puts price sixteenth.

  Background:
    Given user is logged in
    And user opens cars dataset

  Scenario: The dialog opens on what to predict, the predictors and the components
    When user picks "ML > Analyze > PLS..." from the top menu
    Then "PLS" dialog should be visible
    And editor of Predict input in "PLS" dialog should have text "price"
    And editor of Using input in "PLS" dialog should contain text "(15)"
    And Components input in "PLS" dialog should have value "3"
    And Quadratic input in "PLS" dialog should be unchecked
    And RUN button in "PLS" dialog should be enabled

  Scenario: Every column as a predictor includes the predicted one, which RUN does not accept
    When user clicks on editor of Using input in "PLS" dialog
    And user clicks on All label in "Select columns..." dialog
    Then "Select columns..." dialog should contain text "16 checked"
    When user clicks on OK button in "Select columns..." dialog
    Then editor of Using input in "PLS" dialog should contain text "(16)"
    And RUN button in "PLS" dialog should be disabled

  Scenario: With price taken back out, RUN adds three PLS components
    When user clicks on editor of Using input in "PLS" dialog
    And user types "price" into Search input in "Select columns..." dialog
    Then the "text of cell 16 of __name" reading of grid viewer in "Select columns..." dialog should be "price"
    When user clicks on the "cell 16 of x" area of grid viewer in "Select columns..." dialog
    Then the "text of cell 16 of x" reading of grid viewer in "Select columns..." dialog should be "false"
    And "Select columns..." dialog should contain text "15 checked"
    When user clicks on OK button in "Select columns..." dialog
    Then editor of Using input in "PLS" dialog should contain text "(15)"
    And RUN button in "PLS" dialog should be enabled
    When user clicks on RUN button in "PLS" dialog
    Then the top menu command should have completed
    And "PLS" dialog should be hidden
    And 3 new columns should have been added
    And the table should have a column "PLS1"
    And the table should have a column "PLS2"
    And the table should have a column "PLS3"
    And "PLS1" column should have no missing values
    And no error or warning balloon should have been shown
    And no errors should have been logged
