@journey @serial
Feature: Training a classifier on iris in the Train view
  ML | Forge | Train... over iris opened from the demo files (a table with a creation script, so Save offers
  Reference). Translated from the manual scenarios 2.1, 6.2 (steps 1-3), 2.5, 5.1, 5.2 (steps 1-2), 5.13 (steps 1-2),
  2.2, 2.3, 2.8, 4.7 (step 1), 5.5, 5.6 and 5.7 of core/docs/features/platform/predictive-modeling/forge/manual-checks.md.
  The view trains by itself. A training is claimed by the run its results show (data-forge-run-id on the results
  area): "a new training's result should be shown" waits until the view is ready and shows a run other than the one
  remembered, and the server is asked whether that run is the table's newest, with the method and the settings the
  step chose. Server claims read the stand through grok s, scoped to this feature's table and model names.

  Not claimed here: the tab icon, the group labels in one column, the form's height, its start width of about 400 px
  and its kept width after a view switch (visual); Save disabled and the "Choose the target ..." hint before the first
  training (the view trains as it opens, so that state lasts a moment with no signal holding it; the hint is claimed
  for a rejected selection in 2.5); the loader next to Performance (it shows while aria-busy is true, which is what
  the waits read); the tooltips of Table, Target, Features (outside a problem) and the hyperparameters, the Performance
  headers and metric names, the Insights pairs and the Runs cells (claimable as Method's is, left out); the Runs cell
  "F1 (validation)" (the grid reports only the cells it draws, and at 1920 px that column is past the right edge);
  the catalog's other columns of the saved model and the catalog refreshing by itself while open (2.2 steps 1 and 4:
  only the name row is claimed); the Created date; the places of the panels in their rows (5.13: Class balance
  moving between rows 2 and 3; claimable by their boxes, left out); 5.13 step 2's hover of a wrong prediction and
  step 3 (the iris grid is not on screen with the Train view; a cleared cell is the gaps journey of M1c); 6.2 steps
  4-6 (a regression run's metric columns, a removed column, a Table change: they need a column removed and a second
  table); 6.2 step 3's other clicks: on the run already shown (trains nothing), after a cutoff-only move (re-cuts),
  on a cancelled or failed run (trains it again), and the arrow keys on the read-only grid (claimable, left out; the
  package's UI: Runs tests cover most of them); the Runs title tooltip; 5.13 step 1's "no ROC Curve" and step 2's "Predicted Species"
  header tooltip and the scatter plot kept in place, not drawn from scratch (claimable or pixel claims, left out);
  5.5 step 1's table without a known origin (a second table); 4.7 step 1's browser autofill list; the gap under the
  balloon's text (visual). The disabled OK's reason on hover (2.2 step 3) is a @known-failure scenario: a platform
  regression since 2026-10-02/03, confirmed by hand on 2026-10-15.

  Where the product moved on from the manual text: in 2.5 a checked Species turns One-hot encoding on by itself (the
  default for few categories), so the "needs numerical features" mark is claimed after One-hot is unchecked, and Train
  stays hidden (XGBoost retrains by itself); after 6.2's click on a run XGBoost is the user's choice, so 2.5 ends on
  XGBoost and 5.1's invalid setting is XGBoost's Iterations (at most 100). Features lists the columns in the order
  they were checked, so a column checked again goes last and the default model name follows that order.

  Background:
    Given user is logged in
    And no Forge models named "forge-test-iris-{time}, forge-test-iris-ref-{time}, forge-test-iris-copy-{time}, forge-test-iris-svm-{time}, forge-test-iris-link-{time}, forge-test-iris-link2-{time}" are on the server
    And no Forge training runs of the table "forge-test-iris-{time}" are on the server
    And user opens the demo file "iris.csv" as "forge-test-iris-{time}"
    And user watches the task bar
    And user remembers the shown training
    When user picks "ML > Forge > Train..." from the top menu
    Then the "Predictive model" view should be current
    And a new training's result should be shown
    And no balloon of any colour should have been shown
    And no errors should have been logged

  Scenario: The view opens on iris and trains a classifier by itself (2.1)
    Then the task bar should have finished "Training XGBoost model"
    And the shown training should be the newest training run of the table "forge-test-iris-{time}", with:
      | status      | completed      |
      | engine_name | XGBoost        |
      | task        | classification |
      | target_name | Species        |
      | row_count   | 150            |
      | source      | ui             |
      | model_id    | empty          |
    And editor of Table input should have text "forge-test-iris-{time}"
    And editor of Target input should have text "Species"
    And editor of Features input should have text "(4) Sepal.Length, Sepal.Width, Petal.Length, Petal.Width"
    When user clicks on editor of Features input
    Then the "col 1" column should not be checked in the column list of "Select columns..." dialog
    And the "Petal.Width" column should be checked in the column list of "Select columns..." dialog
    When user clicks on CANCEL button in "Select columns..." dialog
    Then the "Select columns..." dialog should close
    And the same training's result should be shown
    And Method input should have value "XGBoost"
    And Method input should offer "CatBoost, Decision Tree, Extra Trees, Hist Gradient Boosting, k-Nearest Neighbors, LDA, Naive Bayes, Neural Network, Random Forest, Softmax, SVM, XGBoost"
    And Iterations input should have value "20"
    And Rate input should have value "0.3"
    And "Max Depth" input should have value "6"
    And Lambda input should have value "1"
    And Alpha input should have value "0"
    And Train button should be hidden
    And "Data" group should be expanded
    And "Method" group should be expanded
    And "Preparation" group should be hidden
    When user collapses "Data" group
    Then "Data" group should be collapsed
    And Target input should be hidden
    When user expands "Data" group
    Then Target input should be visible
    And the "text of cell 1 of Metric" reading of grid viewer in "Performance" results panel should be "Accuracy"
    And the "text of cell 2 of Metric" reading of grid viewer in "Performance" results panel should be "F1"
    And the "rows" reading of grid viewer in "Performance" results panel should be 2
    And the "text of cell 1 of Train" reading of grid viewer in "Performance" results panel should be between 0.95 and 1
    And the "text of cell 1 of Validation" reading of grid viewer in "Performance" results panel should be between 0.85 and 1
    And "Performance" results panel should contain text "Validation: 5-fold cross-validation on 150 rows"
    And "Performance" results panel should contain text "Seed: "
    And "Insights & Tips" results panel should contain text "Features are highly correlated."
    And Save button should be enabled
    When user clicks on seed copy icon
    Then an info balloon containing "copied." should have been shown
    And the clipboard should hold the shown seed
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Every training of the view is a row of Runs; a click on an older one sets its inputs again (6.2)
    Then the "rows" reading of grid viewer in "Runs" results panel should be 1
    And the "text of cell 1 of #" reading of grid viewer in "Runs" results panel should be "1"
    And the "text of cell 1 of Method" reading of grid viewer in "Runs" results panel should be "XGBoost"
    And the "text of cell 1 of Target" reading of grid viewer in "Runs" results panel should be "Species"
    And the "text of cell 1 of Features" reading of grid viewer in "Runs" results panel should be "4"
    And the "text of cell 1 of Preparation" reading of grid viewer in "Runs" results panel should be ""
    And the "text of cell 1 of Rows" reading of grid viewer in "Runs" results panel should be "150"
    And the "text of cell 1 of Status" reading of grid viewer in "Runs" results panel should be "completed"
    And the "text of cell 1 of Seed" reading of grid viewer in "Runs" results panel should be the shown seed
    And the "text of cell 1 of Accuracy (validation)" reading of grid viewer in "Runs" results panel should equal the "text of cell 1 of Validation" reading of grid viewer in "Performance" results panel
    And the "current row" reading of grid viewer in "Runs" results panel should be 1
    When user enters "30" into Iterations input
    Then a new training's result should be shown
    And the shown training should be the newest training run of the table "forge-test-iris-{time}", with:
      | status                     | completed |
      | engine_name                | XGBoost   |
      | hyperparameters.iterations | 30        |
    And the "rows" reading of grid viewer in "Runs" results panel should be 2
    And the "text of cell 1 of #" reading of grid viewer in "Runs" results panel should be "2"
    And the "text of cell 1 of Seed" reading of grid viewer in "Runs" results panel should be the shown seed
    And the "current row" reading of grid viewer in "Runs" results panel should be 1
    When user scrolls "Runs" results panel into view
    And user clicks on the "cell 2 of #" area of grid viewer in "Runs" results panel
    Then Iterations input should have value "20"
    And a new training's result should be shown
    And the shown training should be the newest training run of the table "forge-test-iris-{time}", with:
      | status                     | completed |
      | engine_name                | XGBoost   |
      | hyperparameters.iterations | 20        |
    And the "rows" reading of grid viewer in "Runs" results panel should be 3
    And the "text of cell 1 of #" reading of grid viewer in "Runs" results panel should be "3"
    And the "text of cell 3 of #" reading of grid viewer in "Runs" results panel should be "1"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: A feature the method cannot use is flagged on Features, and nothing trains until it is fixed (2.5)
    When user remembers the training run count of the table "forge-test-iris-{time}"
    And user selects "Petal.Length" in Target input
    Then editor of Features input should be invalid
    And results area should be invalid
    And "Performance" results panel should contain text "Choose the target and the features; the model trains as you change them."
    And Train button should be hidden
    And Save button should be disabled
    And the training run count of the table "forge-test-iris-{time}" should be as remembered
    When user clicks on editor of Features input
    And user toggles the "Species" column in the column list of "Select columns..." dialog
    And user toggles the "Petal.Length" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then editor of Features input should have text "(4) Sepal.Length, Sepal.Width, Petal.Width, Species"
    And "One-hot encoding" input should be checked
    And a new training's result should be shown
    When user remembers the training run count of the table "forge-test-iris-{time}"
    And user unchecks "One-hot encoding" input
    Then editor of Features input should be invalid
    And results area should be invalid
    When user hovers over label of Features input
    Then tooltip should contain text "needs numerical features. Check One-hot encoding in Preparation, or uncheck: Species."
    And editor of Target input should be valid
    And the training run count of the table "forge-test-iris-{time}" should be as remembered
    And no balloon of any colour should have been shown
    When user clicks on editor of Features input
    And user toggles the "Species" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then a new training's result should be shown
    And editor of Features input should be valid
    And Method input should have value "XGBoost"
    And the training run count of the table "forge-test-iris-{time}" should be one more than remembered
    And the shown training should be the newest training run of the table "forge-test-iris-{time}", with:
      | status      | completed    |
      | engine_name | XGBoost      |
      | task        | regression   |
      | target_name | Petal.Length |
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: The methods follow the target; an invalid setting blocks training; the splitter resizes the form (5.1)
    When user hovers over label of Method input
    Then tooltip should contain text "The machine learning method that builds the model."
    And Method input should offer "CatBoost, Decision Tree, Extra Trees, Hist Gradient Boosting, k-Nearest Neighbors, Linear Regression, Neural Network, PLS Regression, Random Forest, SVM, XGBoost"
    And the "text of cell 1 of Metric" reading of grid viewer in "Performance" results panel should be "MSE"
    And the "text of cell 4 of Metric" reading of grid viewer in "Performance" results panel should be "R2"
    When user remembers the training run count of the table "forge-test-iris-{time}"
    And user enters "500" into Iterations input
    Then editor of Iterations input should be invalid
    And results area should be invalid
    And "Performance" results panel should contain text "Fix the settings."
    And Save button should be disabled
    And the training run count of the table "forge-test-iris-{time}" should be as remembered
    When user enters "50" into Iterations input
    Then a new training's result should be shown
    And Save button should be enabled
    And the training run count of the table "forge-test-iris-{time}" should be one more than remembered
    And the shown training should be the newest training run of the table "forge-test-iris-{time}", with:
      | status                     | completed    |
      | engine_name                | XGBoost      |
      | target_name                | Petal.Length |
      | hyperparameters.iterations | 50           |
    When user remembers the size of train form
    And user drags form splitter by 100 and 0 pixels
    Then train form should be wider than remembered
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Every method trains; the method the user chose stays while it can learn from the data (5.2)
    When user selects "Species" in Target input
    Then a new training's result should be shown
    And Method input should have value "XGBoost"
    And editor of Features input should have text "(3) Sepal.Length, Sepal.Width, Petal.Width"
    When user selects "SVM" in Method input
    Then a new training's result should be shown
    And the shown training should be the newest training run of the table "forge-test-iris-{time}", with:
      | status      | completed |
      | engine_name | SVM       |
      | target_name | Species   |
    And the "text of cell 1 of Metric" reading of grid viewer in "Performance" results panel should be "Accuracy"
    When user selects "Softmax" in Method input
    Then a new training's result should be shown
    And the shown training should be the newest training run of the table "forge-test-iris-{time}", with:
      | status      | completed |
      | engine_name | Softmax   |
    When user clicks on editor of Features input
    And user toggles the "Petal.Length" column in the column list of "Select columns..." dialog
    And user toggles the "Petal.Width" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then a new training's result should be shown
    And Method input should have value "Softmax"
    When user selects "Petal.Width" in Target input
    Then an info balloon containing "Softmax cannot be used with this selection; Linear Regression is chosen." should have been shown
    And a new training's result should be shown
    And Method input should have value "Linear Regression"
    And the shown training should be the newest training run of the table "forge-test-iris-{time}", with:
      | status      | completed         |
      | engine_name | Linear Regression |
      | target_name | Petal.Width       |
    When user selects "PLS Regression" in Method input
    Then a new training's result should be shown
    And Components input should have value "3"
    And the shown training should be the newest training run of the table "forge-test-iris-{time}", with:
      | engine_name                | PLS Regression |
      | hyperparameters.components | 3              |
    When user enters "2" into Components input
    Then a new training's result should be shown
    And the shown training should be the newest training run of the table "forge-test-iris-{time}", with:
      | engine_name                | PLS Regression |
      | hyperparameters.components | 2              |
    When user clicks on editor of Features input
    And user toggles the "Sepal.Width" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then a new training's result should be shown
    And Method input should have value "PLS Regression"
    And editor of Features input should have text "(2) Sepal.Length, Petal.Length"
    When user selects "Species" in Target input
    Then an info balloon containing "PLS Regression cannot be used with this selection; XGBoost is chosen." should have been shown
    And a new training's result should be shown
    And the shown training should be the newest training run of the table "forge-test-iris-{time}", with:
      | engine_name | XGBoost |
      | target_name | Species |
    When user clicks on editor of Features input
    And user clicks on All label in "Select columns..." dialog
    And user toggles the "col 1" column in the column list of "Select columns..." dialog
    And user toggles the "Species" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then a new training's result should be shown
    And editor of Features input should have text "(4) Sepal.Length, Petal.Length, Sepal.Width, Petal.Width"
    And Method input should have value "XGBoost"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: A classifier's results; a weaker model lists the rows it gets wrong (5.13)
    Then "Predicted Species vs Actual" results panel should be visible
    And "Confusions" results panel should be visible
    And "Wrong predictions" results panel should be absent
    And the "text of cell 1 of Train" reading of grid viewer in "Performance" results panel should be "1.000"
    And the "rows shown" reading of scatter plot viewer in "Predicted Species vs Actual" results panel should be 150
    And the "bars" reading of bar chart viewer in "Class balance" results panel should be 3
    And the "axes" reading of pc plot viewer in "Features" results panel should be 4
    And "Residuals" results panel should be absent
    And "Statistics" results panel should be absent
    When user enters "1" into Iterations input
    Then a new training's result should be shown
    And the "text of cell 1 of Train" reading of grid viewer in "Performance" results panel should be "0.980"
    And "Wrong predictions" results panel should be visible
    And the "rows" reading of grid viewer in "Wrong predictions" results panel should be 3
    And the "row order" reading of grid viewer in "Wrong predictions" results panel should be "71, 78, 84"
    And the "text of cell 71 of Sepal.Length" reading of grid viewer in "Wrong predictions" results panel should be "5.9"
    And the "text of cell 71 of Sepal.Width" reading of grid viewer in "Wrong predictions" results panel should be "3.2"
    And the "bars" reading of bar chart viewer in "Class balance" results panel should be 3
    When user enters "20" into Iterations input
    Then a new training's result should be shown
    And "Wrong predictions" results panel should be absent
    And no error or warning balloon should have been shown
    And no errors should have been logged

  @known-failure
  Scenario: A disabled OK of Save says why when hovered (2.2 step 3)
    A platform regression since 2026-10-02/03, confirmed by hand on 2026-10-15: no tooltip shows. The reason is a
    ui.setDisabled tooltip, whose hover catcher sits under the dialog (z-index 3000): the pointer reaches the dialog's
    command bar. Drop the tag once the platform is fixed. The scenario stops at its failing step, so the CANCEL below
    keeps the journey clean only once the tooltip shows; until then the harness closes the dialog with Escape.
    When user clicks on Save button
    And user clears Name input in "Save model" dialog
    Then OK button in "Save model" dialog should be disabled
    When user hovers over OK button in "Save model" dialog
    Then tooltip should contain text "Enter a name for the model."
    When user clicks on CANCEL button in "Save model" dialog
    Then the "Save model" dialog should close

  Scenario: Save keeps the model without its data and links the run; the catalog lists it (2.2, 2.3, 2.8, 4.7)
    When user clicks on Save button
    Then "Save model" dialog should be visible
    And Name input in "Save model" dialog should have value "Predict Species by Sepal.Length, Petal.Length, Se..."
    And Description input in "Save model" dialog should have value ""
    And there should be 0 visible tag chips in "Save model" dialog
    And the placeholder of tags box in "Save model" dialog should start with "Type a tag and press Enter"
    And "Save model" dialog should contain text "A link to System:DemoFiles/iris.csv is saved; the data stays where it is."
    And OK button in "Save model" dialog should be enabled
    When user clears Name input in "Save model" dialog
    Then Name input in "Save model" dialog should be invalid
    And OK button in "Save model" dialog should be disabled
    When user types "   " into Name input in "Save model" dialog
    Then OK button in "Save model" dialog should be disabled
    When user presses Enter in Name input in "Save model" dialog
    Then "Save model" dialog should be visible
    When user enters "forge-test-iris-{time}" into Name input in "Save model" dialog
    Then Name input in "Save model" dialog should be valid
    When user types "demo" into tags box in "Save model" dialog
    And user presses Enter in tags box in "Save model" dialog
    Then there should be 1 visible tag chips in "Save model" dialog
    When user selects "None" in "Data storage" input in "Save model" dialog
    Then "Save model" dialog should contain text "Only a summary of the data is saved."
    When user clicks on OK button in "Save model" dialog
    Then the "Save model" dialog should close
    And an info balloon containing 'Model "forge-test-iris-{time}" saved.' should have been shown
    And Save button should be disabled
    And 1 Forge model named "forge-test-iris-{time}" should be on the server
    And the Forge model "forge-test-iris-{time}" on the server should have:
      | engine_name                     | XGBoost        |
      | engine_namespace                | Eda            |
      | engine_kind                     | function       |
      | task                            | classification |
      | target_name                     | Species        |
      | feature_count                   | 4              |
      | row_count                       | 150            |
      | has_training_rows               | false          |
      | storage_mode                    | none           |
      | dataset_table_id                | empty          |
      | dataset_ref                     | empty          |
      | tags                            | demo           |
      | dataset_fingerprint.rowCount    | 150            |
      | dataset_fingerprint.columnCount | 5              |
    And the fingerprint of the Forge model "forge-test-iris-{time}" should describe the columns "Sepal.Length, Petal.Length, Sepal.Width, Petal.Width, Species" and hold no row values
    And the model file folder of "forge-test-iris-{time}" should hold exactly "model.bin"
    And no table named "forge-test-iris-{time}" should be on the server
    And the newest training run of the table "forge-test-iris-{time}" should be linked to the model "forge-test-iris-{time}"
    And the shown seed should be the seed of the Forge model "forge-test-iris-{time}" on the server
    When user switches to the "forge-test-iris-{time}" table view
    And user picks "ML > Forge > Models" from the top menu
    Then the "Forge" view should be current
    And second grid viewer should list the Forge model "forge-test-iris-{time}"
    Given the context panel is open
    When user clicks on the Forge model "forge-test-iris-{time}" in second grid viewer
    Then "Apply model" icon should be enabled
    And "Delete model" icon should be enabled
    And the context panel should show the Forge model "forge-test-iris-{time}"
    When user switches to the "Predictive model" view
    Then the same training's result should be shown
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Data storage at Save: a reference links the file, a copy uploads the training columns (5.5)
    When user enters "0.1" into Alpha input
    Then a new training's result should be shown
    And the shown training should be the newest training run of the table "forge-test-iris-{time}", with:
      | engine_name           | XGBoost |
      | hyperparameters.alpha | 0.1     |
    When user clicks on Save button
    Then "Save model" dialog should contain text "A link to System:DemoFiles/iris.csv is saved; the data stays where it is."
    When user selects "Copy" in "Data storage" input in "Save model" dialog
    Then "Save model" dialog should contain text "The training columns (150 rows) will be uploaded to the server."
    When user selects "Reference" in "Data storage" input in "Save model" dialog
    And user enters "forge-test-iris-ref-{time}" into Name input in "Save model" dialog
    And user clicks on OK button in "Save model" dialog
    Then the "Save model" dialog should close
    And the Forge model "forge-test-iris-ref-{time}" on the server should have:
      | storage_mode     | reference                 |
      | dataset_ref.kind | file                      |
      | dataset_ref.path | System:DemoFiles/iris.csv |
      | dataset_table_id | empty                     |
    When user enters "0" into Alpha input
    Then a new training's result should be shown
    When user clicks on Save button
    And user selects "Copy" in "Data storage" input in "Save model" dialog
    And user enters "forge-test-iris-copy-{time}" into Name input in "Save model" dialog
    And user clicks on OK button in "Save model" dialog
    Then the "Save model" dialog should close
    And the Forge model "forge-test-iris-copy-{time}" on the server should have:
      | storage_mode     | copy    |
      | dataset_table_id | present |
    And the training copy of the Forge model "forge-test-iris-copy-{time}" should be the table "forge-test-iris-copy-{time} (training data)" of 150 rows and the columns "Sepal.Length, Petal.Length, Sepal.Width, Petal.Width, Species"
    When user switches to the "forge-test-iris-{time}" table view
    And user picks "ML > Forge > Models" from the top menu
    Then the "Forge" view should be current
    When user clicks on the Forge model "forge-test-iris-copy-{time}" in second grid viewer
    Then the context panel should show the Forge model "forge-test-iris-copy-{time}"
    And "Details" pane in context panel should contain text "forge-test-iris-copy-{time} (training data)"
    When user clicks on the Forge model "forge-test-iris-ref-{time}" in second grid viewer
    Then the context panel should show the Forge model "forge-test-iris-ref-{time}"
    And "Details" pane in context panel should contain text "System:DemoFiles/iris.csv"
    When user switches to the "Predictive model" view
    And the Forge model "forge-test-iris-copy-{time}" is deleted through the API
    Then no table named "forge-test-iris-copy-{time} (training data)" should be on the server
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: A method whose model file holds training rows is recorded so, without a warning (5.6)
    When user selects "SVM" in Method input
    Then a new training's result should be shown
    And the shown training should be the newest training run of the table "forge-test-iris-{time}", with:
      | status      | completed |
      | engine_name | SVM       |
    When user clicks on Save button
    Then "Save model" dialog should contain text "A link to System:DemoFiles/iris.csv is saved; the data stays where it is."
    And "Save model" dialog should not contain text "training rows"
    When user selects "None" in "Data storage" input in "Save model" dialog
    And user enters "forge-test-iris-svm-{time}" into Name input in "Save model" dialog
    And user clicks on OK button in "Save model" dialog
    Then the "Save model" dialog should close
    And the Forge model "forge-test-iris-svm-{time}" on the server should have:
      | engine_name       | SVM  |
      | has_training_rows | true |
    When user selects "XGBoost" in Method input
    Then a new training's result should be shown
    And the shown training should be the newest training run of the table "forge-test-iris-{time}", with:
      | status      | completed |
      | engine_name | XGBoost   |
    When user clicks on Save button
    Then "Save model" dialog should contain text "A link to System:DemoFiles/iris.csv is saved; the data stays where it is."
    And "Save model" dialog should not contain text "training rows"
    When user clicks on CANCEL button in "Save model" dialog
    Then the "Save model" dialog should close
    And Save button should be enabled
    And "Method" group should not contain text "server"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: The balloon after Save links to the Apply dialog and to the catalog (5.7)
    When user enters "25" into Iterations input
    Then a new training's result should be shown
    Given no balloon is shown
    When user clicks on Save button
    And user selects "None" in "Data storage" input in "Save model" dialog
    And user enters "forge-test-iris-link-{time}" into Name input in "Save model" dialog
    And user clicks on OK button in "Save model" dialog
    Then the "Save model" dialog should close
    When user enters "20" into Iterations input
    Then a new training's result should be shown
    When user clicks on Save button
    And user selects "None" in "Data storage" input in "Save model" dialog
    And user enters "forge-test-iris-link2-{time}" into Name input in "Save model" dialog
    And user clicks on OK button in "Save model" dialog
    Then the "Save model" dialog should close
    When user clicks on "Apply..." link in the balloon of the Forge model "forge-test-iris-link-{time}"
    Then "Apply predictive model" dialog should be visible
    And Model input in "Apply predictive model" dialog should have value "forge-test-iris-link-{time}"
    And editor of Table input in "Apply predictive model" dialog should have text "forge-test-iris-{time}"
    When user clicks on CANCEL button in "Apply predictive model" dialog
    Then the "Apply predictive model" dialog should close
    When user clicks on "Show in the catalog" link in the balloon of the Forge model "forge-test-iris-link2-{time}"
    Then the "Forge" view should be current
    And second grid viewer should list the Forge model "forge-test-iris-link2-{time}"
    And no error or warning balloon should have been shown
    And no errors should have been logged
