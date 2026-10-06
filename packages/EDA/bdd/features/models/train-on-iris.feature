@journey @eda @realizes:ml.menu.models.train-model @realizes:eda.model.softmax @realizes:eda.model.xg-boost
Feature: Training a model to classify iris species
  ML | Models | Train Model... over iris: Species predicted from the four measurements, by the
  Softmax and XGBoost engines of the package. Translated from
  files/TestTrack/EDA/MLMethods/softmax.md and xgboost1.md, and the package's playwright/MLMethods
  specs, which trained through grok.functions.call and asserted only that something came back.

  The preview reports ready only after the latest training, predictions, charts and history are
  complete. A card alone is insufficient: its parameters can change before the result arrives.
  The view remembers hyperparameters, so each scenario sets every value it reads.

  Background:
    Given user is logged in
    And user opens iris dataset

  Scenario: The view takes Species as the target and the four measurements as the features
    When user picks "ML > Models > Train Model..." from the top menu
    Then the "Predictive model" view should be current
    When user selects "Species" in Predict input
    Then editor of Predict input should have text "Species"
    When user clicks on editor of Features input
    And user clicks on All label in "Select columns..." dialog
    Then "Select columns..." dialog should contain text "6 checked"
    And the "text of cell 1 of __name" reading of grid viewer in "Select columns..." dialog should be "col 1"
    And the "text of cell 6 of __name" reading of grid viewer in "Select columns..." dialog should be "Species"
    When user clicks on the "cell 1 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on the "cell 6 of x" area of grid viewer in "Select columns..." dialog
    Then "Select columns..." dialog should contain text "4 checked"
    When user clicks on OK button in "Select columns..." dialog
    Then editor of Features input should contain text "(4)"
    And "One-hot encoding" input should be hidden
    And "Model Engine" input should be visible

  Scenario: Softmax classifies the species and retrains as its hyperparameters move
    When user selects "Eda: Softmax" in "Model Engine" input
    And user enters "100" into Iterations input
    And user enters "2" into Rate input
    And user enters "0.1" into Penalty input
    Then model preview should be ready
    And "Eda: Softmax" heading should be visible
    And "Accuracy" table row should be visible
    And "iterations" table row should contain text "100"
    And "rate" table row should contain text "rate2"
    And "penalty" table row should contain text "penalty0.10"
    And "Predicted Species vs Actual" label should be visible
    And "Confusions" label should be visible
    When user drags the slider of Rate input to 10
    Then model preview should be ready
    And Rate input should have a value between 9.5 and 10.5
    And "rate" table row should not contain text "rate2"
    When user drags the slider of Penalty input to 0.5
    Then model preview should be ready
    And Penalty input should have a value between 0.48 and 0.52
    And "penalty" table row should not contain text "penalty0.10"
    And "Accuracy" table row should be visible
    And the "rows shown" reading of first scatter plot viewer should be 150
    When user drags a selection box over the "view" area of first scatter plot viewer
    Then the "rows selected" reading of first scatter plot viewer should be at least 1
    And first scatter plot viewer should show a selection highlight
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: XGBoost classifies the species and retrains as its clickers and sliders move
    When user selects "Eda: XGBoost" in "Model Engine" input
    And user enters "20" into Iterations input
    And user enters "6" into "Max Depth" input
    And user enters "0.3" into Rate input
    And user enters "1" into Lambda input
    And user enters "0" into Alpha input
    Then model preview should be ready
    And "Eda: XGBoost" heading should be visible
    And "Accuracy" table row should be visible
    And "iterations" table row should contain text "20"
    And "eta" table row should contain text "eta0.30"
    And "lambda" table row should contain text "lambda1"
    And "alpha" table row should contain text "alpha0"
    When user hovers over Iterations input
    And user clicks on plus icon in Iterations input
    Then model preview should be ready
    And "iterations" table row should contain text "21"
    When user hovers over "Max Depth" input
    And user clicks on minus icon in "Max Depth" input
    Then model preview should be ready
    And "maxDepth" table row should contain text "5"
    When user drags the slider of Rate input to 0.5
    Then model preview should be ready
    And Rate input should have a value between 0.48 and 0.52
    And "eta" table row should not contain text "eta0.30"
    When user drags the slider of Lambda input to 50
    Then model preview should be ready
    And Lambda input should have a value between 48 and 52
    And "lambda" table row should not contain text "lambda1"
    When user drags the slider of Alpha input to 40
    Then model preview should be ready
    And Alpha input should have a value between 38 and 42
    And "alpha" table row should not contain text "alpha0"
    And "Accuracy" table row should be visible
    And the "rows shown" reading of first scatter plot viewer should be 150
    When user drags a selection box over the "view" area of first scatter plot viewer
    Then the "rows selected" reading of first scatter plot viewer should be at least 1
    And first scatter plot viewer should show a selection highlight
    And no error or warning balloon should have been shown
    And no errors should have been logged
