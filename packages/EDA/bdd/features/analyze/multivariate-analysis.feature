@journey @eda @realizes:ml.menu.analyze.multivariate-analysis
Feature: Multivariate analysis
  ML | Analyze | Multivariate Analysis... over cars with the dialog's own choices: price predicted
  from the other fifteen numeric columns, three components. Translated from
  files/TestTrack/EDA/multivariate-analysis.md and the package's playwright/multivariate-analysis.test.ts,
  which counted viewer types and left the interactivity of the case to a human.

  The analysis docks its charts over two tables: the scores and the prediction are columns of cars,
  so the grid, Observed vs. Predicted and Scores share one selection; the loadings and the
  coefficients are rows of "cars(Features Analysis)", one per predictor, so Loadings and Regression
  Coefficients share another. Regression Coefficients shares a tab stack with Variable Importance.
  Each title is claimed before a chart is read, so a change in the docking fails loudly.
  In predictor order, horsepower is row eleven of the features table.
  Scores draws its 95% confidence ellipse in orange, so selection paint is checked at the markers.

  Background:
    Given user is logged in
    And user opens cars dataset

  Scenario: Running the analysis adds the scores and the prediction, and docks its charts
    When user picks "ML > Analyze > Multivariate Analysis..." from the top menu
    Then "Multivariate Analysis (PLS)" dialog should be visible
    And editor of Predict input in "Multivariate Analysis (PLS)" dialog should have text "price"
    And editor of Using input in "Multivariate Analysis (PLS)" dialog should contain text "(15)"
    And Components input in "Multivariate Analysis (PLS)" dialog should have value "3"
    When user clicks on RUN button in "Multivariate Analysis (PLS)" dialog
    Then the top menu command should have completed
    And "Multivariate Analysis (PLS)" dialog should be hidden
    And 7 new columns should have been added
    And the table should have a column "price (predicted)"
    And the table should have a column "x.score.t3"
    And "price (predicted)" column should have no missing values
    And "price (predicted)" column should have at least 2 distinct values
    And table "cars" should have 30 rows
    And table "cars(Features Analysis)" should have 15 rows
    And table "cars(Explained Variance)" should have 3 rows
    And the open tableview should have 1 grid viewer
    And the open tableview should have 3 scatter plot viewers
    And the open tableview should have 3 bar chart viewers
    And "Regression Coefficients" tab should be visible
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: A selection in the grid reaches Observed vs. Predicted and Scores
    Then title of second scatter plot viewer should have text "Observed vs. Predicted"
    And title of third scatter plot viewer should have text "Scores"
    And grid viewer should be bound to table "cars"
    And second scatter plot viewer should be bound to table "cars"
    And third scatter plot viewer should be bound to table "cars"
    When user clears the row selection
    Then second scatter plot viewer should show no selection highlight
    And the "rows selected" reading of third scatter plot viewer should be 0
    And the "marker of row 1" area of third scatter plot viewer should not contain the color "#FF8C00"
    And the "marker of row 2" area of third scatter plot viewer should not contain the color "#FF8C00"
    And the "marker of row 3" area of third scatter plot viewer should not contain the color "#FF8C00"
    When user clicks on the "row header 1" area of grid viewer
    And user clicks on the "row header 3" area of grid viewer holding Shift
    Then only rows where "model" is one of "alfaromeo, audi, bmw" should be selected
    And the "rows selected" reading of second scatter plot viewer should be 3
    And the "rows selected" reading of third scatter plot viewer should be 3
    And the "marker of row 1" area of second scatter plot viewer should contain the color "#FF8C00"
    And the "marker of row 2" area of second scatter plot viewer should contain the color "#FF8C00"
    And the "marker of row 3" area of second scatter plot viewer should contain the color "#FF8C00"
    And the "marker of row 1" area of third scatter plot viewer should contain the color "#FF8C00"
    And the "marker of row 2" area of third scatter plot viewer should contain the color "#FF8C00"
    And the "marker of row 3" area of third scatter plot viewer should contain the color "#FF8C00"
    When user clears the row selection
    Then the "rows selected" reading of second scatter plot viewer should be 0
    And the "rows selected" reading of third scatter plot viewer should be 0
    And second scatter plot viewer should show no selection highlight
    And the "marker of row 1" area of third scatter plot viewer should not contain the color "#FF8C00"
    And the "marker of row 2" area of third scatter plot viewer should not contain the color "#FF8C00"
    And the "marker of row 3" area of third scatter plot viewer should not contain the color "#FF8C00"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: A bar of Regression Coefficients selects its predictor in Loadings
    Then title of first scatter plot viewer should have text "Loadings"
    And first scatter plot viewer should be bound to table "cars(Features Analysis)"
    And the "rows selected" reading of first scatter plot viewer should be 0
    When user clicks on "Regression Coefficients" tab
    Then title of first bar chart viewer should have text "Regression Coefficients"
    And first bar chart viewer should be bound to table "cars(Features Analysis)"
    And the "bars" reading of first bar chart viewer should be 15
    When user clicks on the "bar horsepower" area of first bar chart viewer
    Then the "rows selected" reading of first scatter plot viewer should be 1
    And the "marker of row 11" area of first scatter plot viewer should contain the color "#FF8C00"
    And no errors should have been logged
