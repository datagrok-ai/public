@journey @serial
Feature: Training a regressor of wine quality in the Train view
  ML | Forge | Train... over winequality opened from the demo files: a classifier of the default target first, then
  the quality regressor. Translated from the manual scenarios 2.4, 5.2 (step 4), 5.3 and 5.12 (steps 1-3 and 5) of
  core/docs/features/platform/predictive-modeling/forge/manual-checks.md. A training is claimed by the run its results
  show (data-forge-run-id): a new run's result, then that run as the table's newest on the server with the method
  and the settings the step chose; SVM on 6497 rows gets a budget of its own.

  Not claimed here: the panel widths a row starts with and the titles' one style (claimable by their boxes, left out
  / visual); Wrong predictions of the type classifier (it shows only when the model gets a training row wrong, which
  depends on the trees); the Insights and panel title tooltips, the regression line and its equation, the histogram's
  two colours (claimable or pixel claims, left out); 5.12 step 1's classifier metric rows Accuracy ... Negative
  Predicted Value, Class balance without column selectors and a bar's count on hover; 5.12 step 2's Insights lines,
  the failure colour of a weak metric, the Statistics header tooltips, dragging the dividers between panels and rows,
  the area's vertical scroll and the symmetric Residuals range; 2.4 step 3's `type` listed in Features but unchecked
  (claimable, left out); 5.12 step 3's second half, a hovered table row lighting up its
  point (the table view and the Train view are not on screen together, so no person can do it); 5.12 step 4, sizes
  and viewers kept across a retrain to XGBoost (it needs a size the user set and drops the PLS Regression panels;
  step 5's view switch is claimed for every panel); 5.2 step 3, the Table change from iris, which needs the iris
  journey's view; 2.4 step 6, the catalog listing both models (the iris feature claims the catalog). The superseded
  training of 5.3 is claimed three ways: two "Training SVM model" entries in the task bar, exactly one more run, and
  the shown run as the newest, completed with the second value.

  Background:
    Given user is logged in
    And no Forge models named "forge-test-quality-{time}" are on the server
    And no Forge training runs of the table "forge-test-wine-{time}" are on the server
    And user opens the demo file "winequality.csv" as "forge-test-wine-{time}"
    And user remembers the shown training
    When user picks "ML > Forge > Train..." from the top menu
    Then the "Predictive model" view should be current
    And a new training's result should be shown

  Scenario: A classifier of type first, then a regressor of quality; the conflict of Target and Features is shown on Features (2.4, 5.12)
    Then the shown training should be the newest training run of the table "forge-test-wine-{time}", with:
      | status      | completed      |
      | engine_name | XGBoost        |
      | task        | classification |
      | target_name | type           |
      | row_count   | 6497           |
    And editor of Table input should have text "forge-test-wine-{time}"
    And editor of Target input should have text "type"
    And editor of Features input should contain text "(12)"
    And Method input should have value "XGBoost"
    And "Predict probability" input should be visible
    And Save button should be enabled
    And "Predicted type vs Actual" results panel should be visible
    And "Performance" results panel should be visible
    And "Confusions" results panel should be visible
    And "Features" results panel should be visible
    And "Features correlation" results panel should be visible
    And "Class balance" results panel should be visible
    And "ROC Curve" results panel should be absent
    And "Residuals" results panel should be absent
    And the "bars" reading of bar chart viewer in "Class balance" results panel should be 2
    And "Insights & Tips" results panel should contain text "The classes are unbalanced."
    When user remembers the training run count of the table "forge-test-wine-{time}"
    And user selects "quality" in Target input
    Then editor of Features input should be invalid
    And editor of Features input should contain text "(12)"
    And results area should be invalid
    And "Performance" results panel should contain text "Choose the target and the features; the model trains as you change them."
    And Train button should be hidden
    And Save button should be disabled
    And the training run count of the table "forge-test-wine-{time}" should be as remembered
    When user hovers over label of Features input
    Then tooltip should contain text "The target 'quality' is also a feature. Uncheck it in Features."
    And no balloon of any colour should have been shown
    When user clicks on editor of Features input
    And user types "quality" into column search in "Select columns..." dialog
    And user toggles the "quality" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then a new training's result should be shown
    And editor of Features input should be valid
    And Method input should have value "PLS Regression"
    And the shown training should be the newest training run of the table "forge-test-wine-{time}", with:
      | status      | completed      |
      | engine_name | PLS Regression |
      | task        | regression     |
      | target_name | quality        |
    Given user counts the task bar entries
    When user selects "XGBoost" in Method input
    Then a new training's result should be shown
    And the task bar should have shown 1 entry "Fold 1 of 5"
    And the shown training should be the newest training run of the table "forge-test-wine-{time}", with:
      | status      | completed  |
      | engine_name | XGBoost    |
      | task        | regression |
    And the "text of cell 1 of Metric" reading of grid viewer in "Performance" results panel should be "MSE"
    And the "text of cell 2 of Metric" reading of grid viewer in "Performance" results panel should be "RMSE"
    And the "text of cell 3 of Metric" reading of grid viewer in "Performance" results panel should be "MAE"
    And the "text of cell 4 of Metric" reading of grid viewer in "Performance" results panel should be "R2"
    And the "text of cell 4 of Validation" reading of grid viewer in "Performance" results panel should be between 0.3 and 0.7
    And the "text of cell 4 of Train" reading of grid viewer in "Performance" results panel should be higher than the "text of cell 4 of Validation" reading of grid viewer in "Performance" results panel
    And the "text of cell 2 of Validation" reading of grid viewer in "Performance" results panel should be between 0.6 and 0.8
    And "Performance" results panel should contain text "Validation: 5-fold cross-validation on 6497 rows"
    And "Performance" results panel should contain text "Seed: "
    When user clicks on Save button
    Then Name input in "Save model" dialog should have value "Predict quality by fixed.acidity, volatile.acidit..."
    When user selects "None" in "Data storage" input in "Save model" dialog
    And user enters "forge-test-quality-{time}" into Name input in "Save model" dialog
    And user clicks on OK button in "Save model" dialog
    Then the "Save model" dialog should close
    And an info balloon containing 'Model "forge-test-quality-{time}" saved.' should have been shown
    And Save button should be disabled
    And the Forge model "forge-test-quality-{time}" on the server should have:
      | task         | regression |
      | target_name  | quality    |
      | engine_name  | XGBoost    |
      | row_count    | 6497       |
      | storage_mode | none       |
    And the shown seed should be the seed of the Forge model "forge-test-quality-{time}" on the server
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Every regression method trains on the 6497 rows (5.2)
    When user selects "Linear Regression" in Method input
    Then a new training's result should be shown
    And the shown training should be the newest training run of the table "forge-test-wine-{time}", with:
      | status      | completed         |
      | engine_name | Linear Regression |
      | row_count   | 6497              |
    And the "text of cell 4 of Metric" reading of grid viewer in "Performance" results panel should be "R2"
    When user selects "SVM" in Method input
    Then a new training's result should be shown within 180 seconds
    And the shown training should be the newest training run of the table "forge-test-wine-{time}", with:
      | status      | completed |
      | engine_name | SVM       |
      | row_count   | 6497      |
    And the "text of cell 4 of Metric" reading of grid viewer in "Performance" results panel should be "R2"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: A change during a training supersedes it, and only the last one is recorded (5.3)
    Given user counts the task bar entries
    When user remembers the training run count of the table "forge-test-wine-{time}"
    And user sets "Penalty" input to "2" and then to "3" while the training runs
    Then a new training's result should be shown within 180 seconds
    And Penalty input should have value "3"
    And the task bar should have shown 2 entries "Training SVM model"
    And the training run count of the table "forge-test-wine-{time}" should be one more than remembered
    And the shown training should be the newest training run of the table "forge-test-wine-{time}", with:
      | status               | completed |
      | engine_name          | SVM       |
      | hyperparameters.cost | 3         |
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: A regressor's results in rows of panels; a hovered point is the table's row (5.12)
    When user selects "PLS Regression" in Method input
    Then a new training's result should be shown
    And the shown training should be the newest training run of the table "forge-test-wine-{time}", with:
      | status      | completed      |
      | engine_name | PLS Regression |
    And "Predicted quality vs Actual" results panel should be visible
    And "Residuals" results panel should be visible
    And "Distributions of predicted and actual quality" results panel should be visible
    And "Statistics" results panel should be visible
    And "Features" results panel should be visible
    And "Features correlation" results panel should be visible
    And "PLS Regression" results panel should be visible
    And "Class balance" results panel should be absent
    And "Confusions" results panel should be absent
    And the "rows shown" reading of scatter plot viewer in "Predicted quality vs Actual" results panel should be 6497
    And the "rows" reading of grid viewer in "Statistics" results panel should be 11
    And the "axes" reading of pc plot viewer in "Features" results panel should be 13
    And the table "forge-test-wine-{time}" should have no mouse-over row
    When user hovers over the "marker of row 5" area of scatter plot viewer in "Predicted quality vs Actual" results panel
    Then the mouse-over row of the table "forge-test-wine-{time}" should be the "hovered row" reading of scatter plot viewer in "Predicted quality vs Actual" results panel
    When user remembers the sizes of the results panels
    And user switches to the "forge-test-wine-{time}" table view
    And user switches to the "Predictive model" view
    Then the results panels should keep their remembered sizes
    And "Residuals" results panel should be visible
    And no error or warning balloon should have been shown
    And no errors should have been logged
