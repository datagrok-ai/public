@guide @tutorial
Feature: How-to · ML
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). ML > Analyze > Multivariate Analysis... is EDA's PLS
  regression, run in the browser (Crux, no script); the Tutorials ML track walks the same dialog. On cars: Predict
  price, Using every other numeric column (All, then price unticked through the search box: it is below the list's
  fold), three components, Names = model; RUN adds "price (predicted)" and the component columns and docks six
  charts: Observed vs. Predicted, Regression Coefficients and Variable Importance (tabs of one panel), Loadings,
  Scores and Explained Variance. The clip points at their titles; the Loadings and Scores lines are pinned at the
  bottom (placed by themselves they sat at the top, over the Observed vs. Predicted title they don't describe). cars
  is opened off camera; nothing is made on the server.

  Scenario: Find what drives a value with multivariate analysis
    Multivariate analysis (PLS) predicts one column from many others and shows which of them matter most.

    Given user is logged in
    And simple mode is off
    And user opens cars dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    When the video says "A car's price, and fifteen specs that might explain it" at the bottom
    # caption: Multivariate Analysis is under ML > Analyze
    And user picks "ML > Analyze > Multivariate Analysis..." from the top menu
    Then "Multivariate Analysis (PLS)" dialog should be visible
    # caption: Predict: the value to explain
    When user selects "price" in Predict input in "Multivariate Analysis (PLS)" dialog
    # caption: Using: the columns that might explain it
    And user clicks on editor of Using input in "Multivariate Analysis (PLS)" dialog
    Then "Select columns..." dialog should be visible
    When user clicks on "All" link in "Select columns..." dialog
    # caption: Leave out price itself: it is the value to explain
    And user types "price" into Search input in "Select columns..." dialog
    And user toggles the "price" column in the column list of "Select columns..." dialog
    Then "Select columns..." dialog should contain text "15 checked"
    When user clicks on OK button in "Select columns..." dialog
    # caption: Components: how many hidden factors to fit
    And user enters "3" into Components input in "Multivariate Analysis (PLS)" dialog
    # caption: Names: label the points with the car model
    And user selects "model" in Names input in "Multivariate Analysis (PLS)" dialog
    And user clicks on RUN button in "Multivariate Analysis (PLS)" dialog
    Then the table should have a column "price (predicted)"
    And the open tableview should have 3 scatter plot viewers
    And the open tableview should have 3 bar chart viewers
    # caption: Observed vs. Predicted: close to the line means a good fit
    When user hovers over the text "Observed vs. Predicted" in "Observed vs. Predicted" dock panel
    # caption: Variable Importance: the specs that matter most for price
    And user hovers over the text "Variable Importance" in "Variable Importance" dock panel
    When the video says "Loadings: specs close together move together" at the bottom
    And user hovers over the text "Loadings" in "Loadings" dock panel
    When the video says "Scores: each car placed by the fitted factors" at the bottom
    And user hovers over the text "Scores" in "Scores" dock panel
    # caption: Explained Variance: share of price explained by 1, 2 and 3 components
    And user hovers over the text "Explained Variance" in "Explained Variance" dock panel
    Then no error or warning balloon should have been shown
    And no errors should have been logged
