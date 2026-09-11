@journey @eda @realizes:ml.menu.models.train-model @realizes:eda.model.linear-regression @realizes:eda.model.pls-regression @realizes:eda.model.xg-boost
Feature: Training a model to predict the price of a car
  ML | Models | Train Model... over cars: price predicted from every column but price and model, by
  the Linear Regression, PLS Regression and XGBoost engines of the package, one after another in
  the same view. Translated from files/TestTrack/EDA/MLMethods/linear-regression.md,
  pls-regression.md and xgboost2.md, and the package's playwright/MLMethods specs, which trained
  through grok.functions.call and asserted only that something came back.

  The view trains whenever its inputs change; there is no Train button. What it trained is on the
  model card: the engine as its heading, the parameters it was trained with and the scores it
  reached, so a claim on the card is a claim on the model, not on the input that was typed. The
  Features pane is a PC plot over the features, price and the prediction, so its axes count what
  the model was given. The view remembers the hyperparameters of the last session: a scenario sets
  every value it reads.

  The Features picker lists the table's columns in order, model first and price seventeenth. The
  charts PLS adds to the view are its own bar charts, one bar per feature for the coefficients and
  one per component for the explained variance, and its loadings and scores scatter plots, after
  the two every engine shows. The result charts are interactive: rows boxed in Predicted vs Actual are
  selected and highlighted there.

  Background:
    Given user is logged in
    And user opens cars dataset

  Scenario: The view takes price as the target and fifteen columns as the features
    When user picks "ML > Models > Train Model..." from the top menu
    Then the "Predictive model" view should be current
    When user selects "price" in Predict input
    Then editor of Predict input should have text "price"
    When user clicks on editor of Features input
    And user clicks on All label in "Select columns..." dialog
    Then "Select columns..." dialog should contain text "17 checked"
    When user types "price" into Search input in "Select columns..." dialog
    Then the "text of cell 17 of __name" reading of grid viewer in "Select columns..." dialog should be "price"
    When user clicks on the "cell 17 of x" area of grid viewer in "Select columns..." dialog
    And user types "model" into Search input in "Select columns..." dialog
    Then the "text of cell 1 of __name" reading of grid viewer in "Select columns..." dialog should be "model"
    When user clicks on the "cell 1 of x" area of grid viewer in "Select columns..." dialog
    Then "Select columns..." dialog should contain text "15 checked"
    When user clicks on OK button in "Select columns..." dialog
    Then editor of Features input should contain text "(15)"
    And "Model Engine" input should be visible

  Scenario: Linear Regression trains on the fifteen features
    When user selects "Eda: Linear Regression" in "Model Engine" input
    Then "Eda: Linear Regression" heading should be visible
    And "R squared" table row should be visible
    And "Predicted price vs Actual" label should be visible
    And the "axes" reading of pc plot viewer should be 17
    And the "rows selected" reading of first scatter plot viewer should be 0
    When user drags a selection box over the "view" area of first scatter plot viewer
    Then the "rows selected" reading of first scatter plot viewer should be at least 1
    And first scatter plot viewer should show a selection highlight
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: PLS Regression trains with its components and shows its own charts
    When user selects "Eda: PLS Regression" in "Model Engine" input
    And user enters "3" into Components input
    Then "Eda: PLS Regression" heading should be visible
    And "components" table row should contain text "3"
    And "R squared" table row should be visible
    And the "axes" reading of pc plot viewer should be 17
    And the "bars" reading of first bar chart viewer should be 15
    And the "bars" reading of third bar chart viewer should be 3
    And the "rows shown" reading of third scatter plot viewer should be 15
    And the "rows shown" reading of fourth scatter plot viewer should be 30
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: XGBoost retrains as its clickers and sliders move
    When user selects "Eda: XGBoost" in "Model Engine" input
    And user enters "20" into Iterations input
    And user enters "6" into "Max Depth" input
    Then "Eda: XGBoost" heading should be visible
    And "iterations" table row should contain text "20"
    And "maxDepth" table row should contain text "6"
    When user hovers over Iterations input
    And user clicks on plus icon in Iterations input
    Then "iterations" table row should contain text "21"
    When user hovers over "Max Depth" input
    And user clicks on minus icon in "Max Depth" input
    Then "maxDepth" table row should contain text "5"
    When user enters "0.3" into Rate input
    And user enters "1" into Lambda input
    And user enters "0" into Alpha input
    Then "eta" table row should contain text "0.30"
    When user drags the slider of Rate input to 0.5
    Then Rate input should have a value between 0.48 and 0.52
    And "eta" table row should not contain text "eta0.30"
    When user drags the slider of Lambda input to 50
    Then Lambda input should have a value between 48 and 52
    And "lambda" table row should not contain text "lambda1"
    When user drags the slider of Alpha input to 40
    Then Alpha input should have a value between 38 and 42
    And "alpha" table row should not contain text "alpha0"
    And "R squared" table row should be visible
    And no error or warning balloon should have been shown
    And no errors should have been logged
