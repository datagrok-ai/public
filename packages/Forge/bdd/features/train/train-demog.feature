@journey @serial
Feature: Training on demog in the Train view: dates, one-hot, ids, Predict probability, the ROC curve, a slow method
  ML | Forge | Train... over demog opened from the demo files (HEIGHT has 751 gaps, WEIGHT one, AGE one; STARTED is a
  date; USUBJID an id), then a 50000-row demo table. Translated from the manual scenarios 3.14 (step 1), 3.15 (step
  1), 5.8, 5.9, 5.10, 5.14 and 5.4 of core/docs/features/platform/predictive-modeling/forge/manual-checks.md. A
  training is claimed by the run its results show (data-forge-run-id): a new run's result, then that run as the
  table's newest on the server with the method and the preparation it recorded; nothing trained is claimed by the
  table's run count once the view has settled. The one-hot and the probability models are read back on the server by
  their names (features, options, metrics); the ids model is counted there, and its features and skipped columns are
  claimed on its training run. The three are chosen by name in the Apply dialog of an open demog: the one-hot and the
  probability models are applied there, the ids model's dialog is read for its rows and cancelled. The ROC curve is claimed by the
  area under the points it draws against the AUC in its title, so a curve left from an earlier model fails.
  A cutoff re-cuts the shown result without a training: claimed as the same run shown, the same run count, and the
  label metrics of the Train column moving against the ones noted before (Sensitivity down, Specificity up, AUC-ROC
  the same). With Predict probability the Performance panel is narrow next to ROC Curve, so its grid is scrolled to
  the Train column before its cells are read. The 50000-row table's default target is its last column, the boolean
  control; the scenario sets Target to sex first. A disabled Train shows its reason through an overlay that answers a
  moving pointer, not one that lands in a single move: its hovers slide the pointer onto it.

  Not claimed here: 3.15 step 1's Save of "Demog gaps" (claimed
  inside 5.8 without the Save: the same gestures; apply-demog's Background trains the saved model); the exact values
  of oneHotCategories (present); 5.8 step 3's "an unseen value counts as none of the known ones" (claimed as a
  prediction in every row without a gap, the Martian row 1 included); the order of the Columns rows of 5.9 (visible
  and absent only); 5.10's "no progress bar, no loader" (claimed as the same run shown and the run count unchanged),
  the slider (0.9 is typed) and the AUC-ROC tooltip beyond its first sentence; 5.14's hover of the diagonal (a
  formula line is no hit area), the curve that "does not jump or redraw", the ROC panel "in place" and the viewers
  "not drawn from scratch" (pixel claims; the sizes of the panels of every row that keeps its panels across the cutoff
  are claimed: a row is laid out again when Insights & Tips joins it at 0.9, as the manual allows), Confusions and Wrong
  predictions changing at 0.9 (claimable as repaints, left out), the Insights line (it depends on the model), the FPR
  and TPR axis captions, the colours of Class balance and that it counts the rows used only; the grid highlight of a
  hovered point (the table's mouse-over row is read through the API); 5.4's loader next to Performance, Train at the
  right edge of the inputs (visual) and "the training takes a couple of seconds" (the cross is clicked as soon as the
  entry shows).

  Where the product moved on from the manual text: 3.15 step 1 is claimed without saving "Demog gaps"; the Apply
  dialogs of 5.8, 5.9 and 5.10 choose the feature's model by name, since the user's own demog models fit too and the
  dialog presets the newest; Features lists the columns in the order they were checked (5.10 ends on HEIGHT, WEIGHT,
  AGE); the ROC Curve plot has its X and Y column selectors (added after 5.14's "no column selectors" was written).
  The hovered point of 5.14 step 3 is a row of the plot past demog's first skipped row (965), claimed by its values.
  5.4 runs on a demo table of 50000 rows, not the manual's 20000: SVM takes about half a second on 20000 rows on the
  stand, a narrow margin for the cancel, and about 1.5 s on 50000 (the package test 'a training completed after a
  cancelled one clears the invalid mark' cancels SVM on 50000 too).

  Background:
    Given user is logged in
    And no Forge models named "forge-test-demog-onehot-{time}, forge-test-demog-ids-{time}, forge-test-sex-prob-{time}" are on the server
    And no Forge training runs of the tables "forge-test-demog-{time}, forge-test-demog50k-{time}" are on the server
    And user opens the demo file "demog.csv" as "forge-test-demog-{time}"
    And user remembers the shown training
    When user picks "ML > Forge > Train..." from the top menu
    Then the "Predictive model" view should be current
    And a new training's result should be shown
    And no balloon of any colour should have been shown
    And no errors should have been logged

  Scenario: Dates are not default features; a feature or a target the method cannot use is flagged, and nothing trains (3.14)
    Then editor of Target input should have text "SEVERITY"
    And editor of Features input should have text "(3) AGE, HEIGHT, WEIGHT"
    And Method input should have value "XGBoost"
    And the groups of train form should be "Data, Preparation, Method"
    And the "Preparation" group should show exactly the inputs "Missing values"
    And "Missing values" input in "Preparation" group should have "Skip rows" chosen
    When user remembers the training run count of the table "forge-test-demog-{time}"
    And user clicks on editor of Features input
    And user types "STARTED" into column search in "Select columns..." dialog
    Then the "STARTED" column should not be checked in the column list of "Select columns..." dialog
    When user toggles the "STARTED" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then editor of Features input should be invalid
    And Train button in train form should be hidden
    And Method input should have value "XGBoost"
    When user hovers over label of Features input
    Then tooltip should contain text "XGBoost needs numerical features. Uncheck: STARTED."
    And the training run count of the table "forge-test-demog-{time}" should be as remembered
    When user clicks on editor of Features input
    And user types "STARTED" into column search in "Select columns..." dialog
    And user toggles the "STARTED" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then a new training's result should be shown
    When user clicks on editor of Features input
    And user types "SEX" into column search in "Select columns..." dialog
    And user toggles the "SEX" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then "One-hot encoding" input in "Preparation" group should be checked
    And editor of Features input should be valid
    And a new training's result should be shown
    When user remembers the training run count of the table "forge-test-demog-{time}"
    And user unchecks "One-hot encoding" input
    Then editor of Features input should be invalid
    And Method input should have value "XGBoost"
    When user hovers over label of Features input
    Then tooltip should contain text "XGBoost needs numerical features. Check One-hot encoding in Preparation, or uncheck: SEX."
    And the training run count of the table "forge-test-demog-{time}" should be as remembered
    When user clicks on editor of Features input
    And user types "SEX" into column search in "Select columns..." dialog
    And user toggles the "SEX" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then "One-hot encoding" input should be hidden
    And editor of Features input should be valid
    And a new training's result should be shown
    When user remembers the training run count of the table "forge-test-demog-{time}"
    And user selects "STARTED" in Target input
    Then editor of Target input should be invalid
    And Method input should have value "XGBoost"
    When user hovers over label of Target input
    Then tooltip should contain text "XGBoost cannot predict 'STARTED' (type datetime). Choose a numerical, text or boolean target."
    And the training run count of the table "forge-test-demog-{time}" should be as remembered
    When user selects "SEVERITY" in Target input
    Then a new training's result should be shown
    And editor of Target input should be valid
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Text features train through one-hot, and the encoding is replayed with the training categories (5.8, 3.15)
    When user remembers the training run count of the table "forge-test-demog-{time}"
    And user selects "AGE" in Target input
    Then editor of Features input should be invalid
    When user hovers over label of Features input
    Then tooltip should contain text "The target 'AGE' is also a feature. Uncheck it in Features."
    And the training run count of the table "forge-test-demog-{time}" should be as remembered
    When user clicks on editor of Features input
    And user toggles the "AGE" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then a new training's result should be shown
    And Method input should have value "Linear Regression"
    And the shown training should be the newest training run of the table "forge-test-demog-{time}", with:
      | engine_name                       | Linear Regression  |
      | target_name                       | AGE                |
      | row_count                         | 5097               |
      | options.preprocessingInfo         | ["ignore-missing"] |
      | options.missingValues.skippedRows | 753                |
    And "Performance" results panel should contain text "Rows: 5097 used, 753 skipped (missing values)"
    And the "text of cell 1 of Metric" reading of grid viewer in "Performance" results panel should be "MSE"
    And the "text of cell 4 of Metric" reading of grid viewer in "Performance" results panel should be "R2"
    When user hovers over label of "Missing values" input
    Then tooltip should contain text "Rows with missing values in: HEIGHT: 751, WEIGHT: 1."
    When user clicks on editor of Features input
    And user types "SEX" into column search in "Select columns..." dialog
    And user toggles the "SEX" column in the column list of "Select columns..." dialog
    And user types "RACE" into column search in "Select columns..." dialog
    And user toggles the "RACE" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then "One-hot encoding" input in "Preparation" group should be visible
    And "One-hot encoding" input should be checked
    And editor of Features input should be valid
    And a new training's result should be shown
    And Method input should have value "PLS Regression"
    And the "text of cell 1 of Metric" reading of grid viewer in "Performance" results panel should be "MSE"
    And the "text of cell 2 of Metric" reading of grid viewer in "Performance" results panel should be "RMSE"
    And the "text of cell 3 of Metric" reading of grid viewer in "Performance" results panel should be "MAE"
    And the "text of cell 4 of Metric" reading of grid viewer in "Performance" results panel should be "R2"
    When user hovers over label of "One-hot encoding" input
    Then tooltip should contain text "Lets the model use categories as numbers without implying an order: each category gets its own 0/1 column."
    When user remembers the training run count of the table "forge-test-demog-{time}"
    And user unchecks "One-hot encoding" input
    Then editor of Features input should be invalid
    When user hovers over label of Features input
    Then tooltip should contain text "PLS Regression needs numerical features. Check One-hot encoding in Preparation, or uncheck: SEX, RACE."
    And the training run count of the table "forge-test-demog-{time}" should be as remembered
    When user checks "One-hot encoding" input
    Then editor of Features input should be valid
    And a new training's result should be shown
    When user clicks on Save button
    And user selects "None" in "Data storage" input in "Save model" dialog
    And user enters "forge-test-demog-onehot-{time}" into Name input in "Save model" dialog
    And user clicks on OK button in "Save model" dialog
    Then the "Save model" dialog should close
    And the Forge model "forge-test-demog-onehot-{time}" on the server should have:
      | engine_name                   | PLS Regression               |
      | target_name                   | AGE                          |
      | features.columns.2.name       | SEX                          |
      | features.columns.2.type       | string                       |
      | features.columns.3.name       | RACE                         |
      | features.columns.3.type       | string                       |
      | options.preprocessingInfo     | ["ignore-missing","one-hot"] |
      | options.oneHotCategories.SEX  | present                      |
      | options.oneHotCategories.RACE | present                      |
    Given user counts the task bar entries
    And user opens the demo file "demog.csv" as "forge-test-demog2-{time}"
    When user sets "RACE" column in row 1 to "Martian"
    And user picks "ML > Forge > Apply..." from the top menu
    Then there should be 1 visible "Apply predictive model" dialog
    When user selects "forge-test-demog-onehot-{time}" in Model input in "Apply predictive model" dialog
    Then Model input in "Apply predictive model" dialog should have value "forge-test-demog-onehot-{time}"
    And summary of "Columns" group should have the text "4 of 4 matched"
    When user clicks on OK button in "Apply predictive model" dialog
    Then the "Apply predictive model" dialog should close
    And the task bar entry that showed "Predicting AGE" should have ended
    And the table should have a column "AGE (predicted)"
    And an info balloon containing 'Added the column "AGE (predicted)" to forge-test-demog2-{time}' should have been shown
    And the "AGE (predicted)" column should be empty exactly in the rows where any of "HEIGHT, WEIGHT, SEX, RACE" columns is empty
    And the value of "AGE (predicted)" column in row 1 should match "^-?\d"
    When user switches to the "Predictive model" view
    Then the same training's result should be shown
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: An id-like text column is left out, and applying the model does not ask for it (5.9)
    When user clicks on editor of Features input
    And user types "USUBJID" into column search in "Select columns..." dialog
    And user toggles the "USUBJID" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then "Skip unique categories" input in "Preparation" group should be visible
    And "Skip unique categories" input should be checked
    And a new training's result should be shown
    And Method input should have value "PLS Regression"
    And the shown training should be the newest training run of the table "forge-test-demog-{time}", with:
      | features.columns.4.name   | USUBJID                                               |
      | options.preprocessingInfo | ["ignore-missing","skip-unique-categories","one-hot"] |
      | options.skippedColumns    | ["USUBJID"]                                           |
    When user hovers over label of "Skip unique categories" input
    Then tooltip should contain text "Leaves out text features whose values are all different, such as ids."
    When user clicks on Save button
    And user selects "None" in "Data storage" input in "Save model" dialog
    And user enters "forge-test-demog-ids-{time}" into Name input in "Save model" dialog
    And user clicks on OK button in "Save model" dialog
    Then the "Save model" dialog should close
    And 1 Forge model named "forge-test-demog-ids-{time}" should be on the server
    When user switches to the "forge-test-demog2-{time}" table view
    And user picks "ML > Forge > Apply..." from the top menu
    Then there should be 1 visible "Apply predictive model" dialog
    When user selects "forge-test-demog-ids-{time}" in Model input in "Apply predictive model" dialog
    Then Model input in "Apply predictive model" dialog should have value "forge-test-demog-ids-{time}"
    And summary of "Columns" group should have the text "4 of 4 matched"
    When user expands "Columns" group
    Then "HEIGHT" input in "Apply predictive model" dialog should be visible
    And "WEIGHT" input in "Apply predictive model" dialog should be visible
    And "SEX" input in "Apply predictive model" dialog should be visible
    And "RACE" input in "Apply predictive model" dialog should be visible
    And "USUBJID" input in "Apply predictive model" dialog should be absent
    When user clicks on CANCEL button in "Apply predictive model" dialog
    Then the "Apply predictive model" dialog should close
    When user switches to the "Predictive model" view
    And user clicks on editor of Features input
    And user types "USUBJID" into column search in "Select columns..." dialog
    And user toggles the "USUBJID" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then "Skip unique categories" input should be hidden
    And a new training's result should be shown
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Predict probability trains a regressor on the two classes; the cutoff re-cuts without a training (5.10)
    When user clicks on editor of Features input
    And user types "SEX" into column search in "Select columns..." dialog
    And user toggles the "SEX" column in the column list of "Select columns..." dialog
    And user types "RACE" into column search in "Select columns..." dialog
    And user toggles the "RACE" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then "One-hot encoding" input should be hidden
    And a new training's result should be shown
    When user selects "SEX" in Target input
    Then a new training's result should be shown
    When user clicks on editor of Features input
    And user toggles the "AGE" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then a new training's result should be shown
    And editor of Features input should have text "(3) HEIGHT, WEIGHT, AGE"
    And the "Preparation" group should show exactly the inputs "Missing values, Predict probability"
    And "Predict probability" input should not be checked
    And Method input should have value "XGBoost"
    And Method input should offer "Softmax, SVM, XGBoost"
    And the "rows" reading of grid viewer in "Performance" results panel should be 6
    And the "text of cell 1 of Metric" reading of grid viewer in "Performance" results panel should be "Accuracy"
    And the "text of cell 2 of Metric" reading of grid viewer in "Performance" results panel should be "F1"
    And the "text of cell 3 of Metric" reading of grid viewer in "Performance" results panel should be "Sensitivity"
    And the "text of cell 4 of Metric" reading of grid viewer in "Performance" results panel should be "Specificity"
    And the "text of cell 5 of Metric" reading of grid viewer in "Performance" results panel should be "Precision"
    And the "text of cell 6 of Metric" reading of grid viewer in "Performance" results panel should be "Negative Predicted Value"
    And "Performance" results panel should contain text "Positive class: F"
    When user hovers over label of "Predict probability" input
    Then tooltip should contain text "Trains a regression model on the two classes and predicts the probability of F. The score is the regression output, not a calibrated probability."
    When user checks "Predict probability" input
    Then "Positive class cutoff" input in "Preparation" group should be visible
    And "Positive class cutoff" input should have value "0.5"
    And a new training's result should be shown
    And Method input should have value "Linear Regression"
    And Method input should offer "Linear Regression, PLS Regression, SVM, XGBoost"
    And the "rows" reading of grid viewer in "Performance" results panel should be 7
    And the "text of cell 7 of Metric" reading of grid viewer in "Performance" results panel should be "AUC-ROC"
    When user hovers over label of "Positive class cutoff" input
    Then tooltip should contain text "Probability from which a row is predicted as F."
    When user hovers over the "cell 7 of Metric" area of grid viewer in "Performance" results panel
    Then tooltip should contain text "Area Under the Receiver Operating Characteristic Curve"
    When user scrolls grid viewer in "Performance" results panel to the "Train" column
    Then the "text of cell 7 of Train" reading of grid viewer in "Performance" results panel should be between 0.5 and 1
    And the "text of cell 7 of Validation" reading of grid viewer in "Performance" results panel should be between 0.5 and 1
    When user notes the "text of cell 3 of Train" reading of grid viewer in "Performance" results panel
    And user notes the "text of cell 4 of Train" reading of grid viewer in "Performance" results panel
    And user remembers the "text of cell 3 of Train" reading of grid viewer in "Performance" results panel
    And user remembers the "text of cell 7 of Train" reading of grid viewer in "Performance" results panel
    And user remembers the training run count of the table "forge-test-demog-{time}"
    And user enters "0.9" into "Positive class cutoff" input
    Then the "text of cell 3 of Train" reading of grid viewer in "Performance" results panel should be lower than noted
    And the "text of cell 4 of Train" reading of grid viewer in "Performance" results panel should be higher than noted
    And the "text of cell 7 of Train" reading of grid viewer in "Performance" results panel should be as remembered
    And the same training's result should be shown
    When user enters "1.5" into "Positive class cutoff" input
    Then editor of "Positive class cutoff" input should be invalid
    And "Performance" results panel should contain text "Fix the settings."
    And Save button should be disabled
    When user enters "0.5" into "Positive class cutoff" input
    Then the "text of cell 3 of Train" reading of grid viewer in "Performance" results panel should be as remembered
    And the same training's result should be shown
    And Save button should be enabled
    When user enters "1.5" into "Positive class cutoff" input
    Then "Performance" results panel should contain text "Fix the settings."
    When user clicks on editor of Features input
    And user toggles the "AGE" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then editor of Features input should have text "(2) HEIGHT, WEIGHT"
    When user clicks on editor of Features input
    And user toggles the "AGE" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then editor of Features input should have text "(3) HEIGHT, WEIGHT, AGE"
    And results area should be invalid
    And "Performance" results panel should contain text "Fix the settings."
    And the training run count of the table "forge-test-demog-{time}" should be as remembered
    When user enters "0.5" into "Positive class cutoff" input
    Then a new training's result should be shown
    And the "text of cell 7 of Metric" reading of grid viewer in "Performance" results panel should be "AUC-ROC"
    When user clicks on Save button
    And user selects "None" in "Data storage" input in "Save model" dialog
    And user enters "forge-test-sex-prob-{time}" into Name input in "Save model" dialog
    And user clicks on OK button in "Save model" dialog
    Then the "Save model" dialog should close
    And the Forge model "forge-test-sex-prob-{time}" on the server should have:
      | task                                  | regression                |
      | engine_name                           | Linear Regression         |
      | options.postprocessingInfo            | ["binary-classification"] |
      | options.positiveClass                 | F                         |
      | options.negativeClass                 | M                         |
      | options.binaryClassificationThreshold | 0.5                       |
      | options.targetType                    | string                    |
      | metrics.validation.auc                | number                    |
    Given user counts the task bar entries
    When user switches to the "forge-test-demog-{time}" table view
    And user picks "ML > Forge > Apply..." from the top menu
    Then there should be 1 visible "Apply predictive model" dialog
    When user selects "forge-test-sex-prob-{time}" in Model input in "Apply predictive model" dialog
    Then Model input in "Apply predictive model" dialog should have value "forge-test-sex-prob-{time}"
    And summary of "Columns" group should have the text "3 of 3 matched"
    When user clicks on OK button in "Apply predictive model" dialog
    Then the "Apply predictive model" dialog should close
    And the task bar entry that showed "Predicting SEX" should have ended
    And the table should have a column "SEX (predicted)"
    And the "SEX (predicted)" column should be empty exactly in the rows where any of "AGE, HEIGHT, WEIGHT" columns is empty
    And the value of "SEX (predicted)" column in row 1 should match "^(F|M)$"
    And every value of "SEX (predicted)" column should match "^(F|M)$"
    When user switches to the "Predictive model" view
    Then the same training's result should be shown
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: The ROC curve follows the scores; a cutoff only moves its marked point (5.14)
    When user scrolls grid viewer in "Performance" results panel to the "Train" column
    Then "Predicted SEX vs Actual" results panel should be visible
    And ROC Curve panel should be visible
    And the AUC in the title of ROC Curve panel should equal the "text of cell 7 of Train" reading of grid viewer in "Performance" results panel
    And the area under the curve of scatter plot viewer in ROC Curve panel should equal the AUC in the title of ROC Curve panel within 0.02
    And the current row of scatter plot viewer in ROC Curve panel should be 51
    And "Confusions" results panel should be visible
    And "Class balance" results panel should be visible
    And "Features" results panel should be visible
    And "Features correlation" results panel should be visible
    And the "bars" reading of bar chart viewer in "Class balance" results panel should be 2
    When user hovers over the "marker of row 20" area of scatter plot viewer in ROC Curve panel
    Then tooltip should contain text "Threshold"
    When user remembers the sizes of the results panels
    And user enters "0.9" into "Positive class cutoff" input
    Then the current row of scatter plot viewer in ROC Curve panel should be 91
    And the same training's result should be shown
    And the AUC in the title of ROC Curve panel should equal the "text of cell 7 of Train" reading of grid viewer in "Performance" results panel
    And the results panels should keep their remembered sizes
    When user selects "XGBoost" in Method input
    Then a new training's result should be shown
    And Train button in train form should be hidden
    And ROC Curve panel should be visible
    When user scrolls grid viewer in "Performance" results panel to the "Train" column
    Then the AUC in the title of ROC Curve panel should equal the "text of cell 7 of Train" reading of grid viewer in "Performance" results panel
    And the area under the curve of scatter plot viewer in ROC Curve panel should equal the AUC in the title of ROC Curve panel within 0.02
    And the current row of scatter plot viewer in ROC Curve panel should be 91
    When user enters "0.5" into "Positive class cutoff" input
    Then the current row of scatter plot viewer in ROC Curve panel should be 51
    And the "rows shown" reading of scatter plot viewer in "Predicted SEX vs Actual" results panel should be 5097
    And the table "forge-test-demog-{time}" should have no mouse-over row
    When user hovers over the point of row 1000 of scatter plot viewer in "Predicted SEX vs Actual" results panel
    Then the mouse-over row of the table "forge-test-demog-{time}" should hold the "AGE, HEIGHT, WEIGHT" values of the hovered row of scatter plot viewer in "Predicted SEX vs Actual" results panel
    When user closes the current view
    Then no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: A slow method trains on Train only; a cancel is recorded (5.4)
    Given a demo table "forge-test-demog50k-{time}" of 50000 rows is open
    When user picks "ML > Forge > Train..." from the top menu
    Then the "Predictive model" view should be current
    And a new training's result should be shown
    When user selects "sex" in Target input
    Then a new training's result should be shown
    When user remembers the training run count of the table "forge-test-demog50k-{time}"
    And user selects "SVM" in Method input
    Then "Performance" results panel should contain text "SVM takes a while on 50000 rows, so it does not retrain on every change. Click Train."
    And Train button in train form should be visible
    And Train button in train form should be enabled
    And the training run count of the table "forge-test-demog50k-{time}" should be as remembered
    When user hovers over Train button in train form
    Then tooltip should contain text "Train the model on the current selection."
    Given user counts the task bar entries
    When user clicks on Train button in train form
    And user cancels the "Training SVM model" task bar entry
    Then a warning balloon containing "Training was cancelled." should have been shown
    And results area should be invalid
    And "Performance" results panel should contain text "SVM takes a while on 50000 rows, so it does not retrain on every change. Click Train."
    And Train button in train form should be enabled
    And Save button should be disabled
    And the newest training run of the table "forge-test-demog50k-{time}" should have:
      | status      | cancelled |
      | engine_name | SVM       |
      | row_count   | 50000     |
    When user remembers the training run count of the table "forge-test-demog50k-{time}"
    And user clicks on Train button in train form
    Then a new training's result should be shown within 300 seconds
    And the training run count of the table "forge-test-demog50k-{time}" should be one more than remembered
    And the shown training should be the newest training run of the table "forge-test-demog50k-{time}", with:
      | status      | completed |
      | engine_name | SVM       |
      | row_count   | 50000     |
    And Save button should be enabled
    And Train button in train form should be disabled
    When user slides the pointer onto Train button in train form
    Then tooltip should contain text "The model is trained on this selection. Change an input to train again."
    When user remembers the training run count of the table "forge-test-demog50k-{time}"
    And user enters "2" into Penalty input
    Then Train button in train form should be enabled
    And the training run count of the table "forge-test-demog50k-{time}" should be as remembered
    When user enters "5000" into Penalty input
    Then editor of Penalty input should be invalid
    And results area should be invalid
    And "Performance" results panel should contain text "Fix the settings."
    And Train button in train form should be visible
    And Train button in train form should be disabled
    When user slides the pointer onto Train button in train form
    Then tooltip should contain text "Penalty: Value must be less than 1000"
    When user enters "2" into Penalty input
    Then Train button in train form should be enabled
    When user selects "XGBoost" in Method input
    Then Train button in train form should be hidden
    And a new training's result should be shown
    And the "text of cell 1 of Metric" reading of grid viewer in "Performance" results panel should be "Accuracy"
    When user closes the current view
    Given the table "forge-test-demog50k-{time}" is closed
    Then no error or warning balloon should have been shown
    And no errors should have been logged
