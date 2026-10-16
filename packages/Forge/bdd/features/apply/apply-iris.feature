@journey @serial
Feature: Applying an iris classifier: the Apply dialog, column matching, missing values and the API
  ML | Forge | Apply..., the catalog's Apply model icon and Forge:applyModel over iris opened from the demo files, with
  an iris classifier the Background trains and saves through Forge:trainModel (XGBoost, no data stored) in place of
  the manual's Train view and Save (train-iris.feature walks those). Translated from the manual scenarios 3.1 (steps
  3-5), 3.3, 3.4, 3.5 (steps 1-2), 3.7 (steps 1, 2, 4), 3.8, 3.9, 3.13, 3.14 (steps 2-5) and 4.10 of
  core/docs/features/platform/predictive-modeling/forge/manual-checks.md. The predictions are claimed against the
  training: the share of rows where Species and its prediction agree (an empty prediction fails it) equals the stored
  train Accuracy. An application's end is its task bar entry gone (followed as an element: the entry shows the
  progress instead of "Predicting Species" while it runs), and each scenario counts the entries it caused. Every
  application is counted on the server for this feature's model, so a count is exact; 3.13's Impute is claimed by
  what it predicts (copies of the features imputed with Eda:knnImpute and the dialog's settings, predicted with the
  same model through Forge:applyModel), a check that records one more application after 3.13's counts, so the counts
  of the later scenarios include it. The second iris copy (3.3), the copy with empty cells (3.13) and demog (3.8, 3.14)
  are opened by their scenario and closed by it.

  Not claimed here: 3.1 steps 1-2 (the Background trains through the API: the model the view would save) and "the
  other columns' values unchanged" (the package's Apply tests); 3.3 step 2's "listed below the fitting models" (the
  model is chosen by name, and the order of the stand's other models is not the feature's); 3.5 step 2's "the previous
  row is scenario 3.3's" (the newest application is claimed after each scenario's own application instead) and step 3
  (the Domains grid); 3.7 step 1 (the row is made with grok s, not through the Domains view's New model...) and step 3
  (catalog.feature deletes its own hand-made row); 3.8's console output and "no dialog opens" (the refusal is the
  call's own error); 3.9's "no Predicted By entry in the tags table" (a text claim is case-insensitive, so Forge's
  Predicted by would match it); 3.13's tooltip of the Sepal.Width row beyond its missing-values sentence, Impute's
  other settings, and the values Impute filled in: the two rows' predicted classes do not depend on them (the check
  with 1 neighbour instead of 4 passes as well), so here the check proves only that Impute predicted those rows as the
  imputed copy does (apply-demog.feature's regressor tells the imputations apart); 3.14 step 1 (the Train view on demog, M1c) and "newest first" in step 2 (the model is chosen by
  name); 4.10 step 1's card lines beyond the target, the method and the name, and the tags table; 4.10's "no Predicted
  by pane" on Species is read once the panel's title is Species (no signal says the panel has added its last pane).
  The reason the disabled OK gives on hover (3.3 step 4, 3.7 step 2, 3.14 steps 2 and 4): the platform shows no tooltip
  of a disabled button inside a dialog (a regression, train-iris.feature's @known-failure scenario); OK is claimed
  disabled and the reason is claimed on the input it concerns.

  Where the product moved on from the manual text: the catalog is ML | Forge | Models (it focuses an open catalog);
  3.4's Apply model icon presets the table the model fits, so no table is chosen by hand. 3.8 step 2's refusal is
  claimed on demog, not winequality: any table without the model's columns is refused with the same message, and
  demog is the table the feature opens anyway (3.14).

  Background:
    Given user is logged in
    And no Forge models named "forge-test-apply-iris-{time}, forge-test-apply-hand-{time}" are on the server
    And no Forge training runs of the table "forge-test-iris-{time}" are on the server
    And user opens the demo file "iris.csv" as "forge-test-iris-{time}"
    And a Forge model "forge-test-apply-iris-{time}" is trained on the current table through the API with:
      | target | Species |
      | method | XGBoost |

  Scenario: The newest fitting model is preset; the predictions agree with the training as its Accuracy says (3.1, 3.5)
    Given user counts the task bar entries
    When user picks "ML > Forge > Apply..." from the top menu
    Then there should be 1 visible "Apply predictive model" dialog
    And Table input in "Apply predictive model" dialog should have value "forge-test-iris-{time}"
    And Model input in "Apply predictive model" dialog should have value "forge-test-apply-iris-{time}"
    And summary of "Columns" group should have the text "4 of 4 matched"
    And summary of "Columns" group should be valid
    And "Columns" group should be collapsed
    And "Missing values" input in "Apply predictive model" dialog should be hidden
    And "More options" group should be collapsed
    And OK button in "Apply predictive model" dialog should be enabled
    When user expands "Columns" group
    Then editor of "Sepal.Length" input in "Apply predictive model" dialog should have text "Sepal.Length"
    And editor of "Sepal.Width" input in "Apply predictive model" dialog should have text "Sepal.Width"
    And editor of "Petal.Length" input in "Apply predictive model" dialog should have text "Petal.Length"
    And editor of "Petal.Width" input in "Apply predictive model" dialog should have text "Petal.Width"
    When user expands "More options" group
    Then "Batch size" input in "Apply predictive model" dialog should have value "10000"
    When user clicks on OK button in "Apply predictive model" dialog
    Then the "Apply predictive model" dialog should close
    And the task bar entry that showed "Predicting Species" should have ended
    And the task bar should have shown 1 entry "Predicting Species"
    And a new column "Species (predicted)" should have been added
    And "Species (predicted)" column should have no missing values
    And an info balloon containing 'Added the column "Species (predicted)" to forge-test-iris-{time}.' should have been shown
    And the table should have the columns "col 1, Sepal.Length, Sepal.Width, Petal.Length, Petal.Width, Species, Species (predicted)"
    And the share of rows where "Species" and "Species (predicted)" columns of the table "forge-test-iris-{time}" agree should equal the Train Accuracy of the Forge model "forge-test-apply-iris-{time}" within 0.001
    And the model "forge-test-apply-iris-{time}" should have 1 application on the server
    And the newest application of the model "forge-test-apply-iris-{time}" should have:
      | table_name   | forge-test-iris-{time} |
      | row_count    | 150                    |
      | skipped_rows | 0                      |
      | column_name  | Species (predicted)    |
      | source       | ui                     |
      | status       | completed              |
      | duration_ms  | number                 |
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: The catalog's Apply model icon applies under a new name; the column points at its model (3.4, 4.10)
    Given user counts the task bar entries
    When user picks "ML > Forge > Models" from the top menu
    Then the "Forge" view should be current
    And "Apply model" icon in forge app should be disabled
    Given the context panel is open
    When user clicks on the Forge model "forge-test-apply-iris-{time}" in second grid viewer
    Then "Apply model" icon in forge app should be enabled
    When user hovers over "Apply model" icon in forge app
    Then tooltip should contain text "Apply model"
    When user clicks on "Apply model" icon in forge app
    Then there should be 1 visible "Apply predictive model" dialog
    And Model input in "Apply predictive model" dialog should have value "forge-test-apply-iris-{time}"
    And Table input in "Apply predictive model" dialog should have value "forge-test-iris-{time}"
    When user clicks on OK button in "Apply predictive model" dialog
    Then the "Apply predictive model" dialog should close
    And the task bar entry that showed "Predicting Species" should have ended
    And the task bar should have shown 1 entry "Predicting Species"
    And the "forge-test-iris-{time}" view should be current
    And the table should have a column "Species (predicted 2)"
    And an info balloon containing 'Added the column "Species (predicted 2)" to forge-test-iris-{time}.' should have been shown
    And the "Species (predicted 2)" column should be tagged with the id of the Forge model "forge-test-apply-iris-{time}"
    And the "Species (predicted 2)" column should not have the tag "predictive.model"
    And the model "forge-test-apply-iris-{time}" should have 2 applications on the server
    Given the "Species (predicted 2)" column is the current object
    Then the context panel title should be "Species (predicted 2)"
    When user expands "Predicted by" pane in context panel
    Then "Predicted by" pane in context panel should contain text "forge-test-apply-iris-{time}"
    And "Predicted by" pane in context panel should contain text "Predict Species"
    And "Predicted by" pane in context panel should contain text "using XGBoost"
    When user hovers over "Predict Species" label in "Predicted by" pane in context panel
    Then tooltip should contain text "Predicts Species from 4 features with XGBoost"
    Given the "Species" column is the current object
    Then the context panel title should be "Species"
    And "Predicted by" pane in context panel should be absent
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Close names are matched, a far one is not guessed, a column is never used twice, the table is not renamed (3.3)
    Given user counts the task bar entries
    And user opens the demo file "iris.csv" as "forge-test-iris2-{time}"
    When user renames "Sepal.Width" column to "sepal width"
    And user removes "Petal.Width" column
    And user picks "ML > Forge > Apply..." from the top menu
    Then there should be 1 visible "Apply predictive model" dialog
    When user selects "forge-test-apply-iris-{time}" in Model input in "Apply predictive model" dialog
    Then summary of "Columns" group should have the text "3 of 4 matched"
    And summary of "Columns" group should be invalid
    And "Columns" group should be expanded
    And editor of "Sepal.Length" input in "Apply predictive model" dialog should have text "Sepal.Length"
    And editor of "Sepal.Width" input in "Apply predictive model" dialog should have text "sepal width"
    And editor of "Petal.Length" input in "Apply predictive model" dialog should have text "Petal.Length"
    And editor of "Petal.Width" input in "Apply predictive model" dialog should have text ""
    And editor of "Petal.Width" input in "Apply predictive model" dialog should be invalid
    And OK button in "Apply predictive model" dialog should be disabled
    When user hovers over label of "Petal.Width" input in "Apply predictive model" dialog
    Then tooltip should contain text "Table column used as the feature 'Petal.Width' (numbers)."
    And tooltip should contain text "Choose a column for 'Petal.Width'."
    And the column picker of "Petal.Width" input in "Apply predictive model" dialog should offer exactly "col 1, Sepal.Length, sepal width, Petal.Length"
    And editor of "Petal.Width" input in "Apply predictive model" dialog should have text ""
    When user selects "col 1" in "Petal.Width" input in "Apply predictive model" dialog
    Then editor of "Petal.Width" input in "Apply predictive model" dialog should be valid
    And summary of "Columns" group should have the text "4 of 4 matched"
    And summary of "Columns" group should be valid
    And "Columns" group should be expanded
    And OK button in "Apply predictive model" dialog should be enabled
    When user selects "col 1" in "Petal.Length" input in "Apply predictive model" dialog
    Then editor of "Petal.Width" input in "Apply predictive model" dialog should be invalid
    And summary of "Columns" group should have the text "3 of 4 matched"
    And summary of "Columns" group should be invalid
    And OK button in "Apply predictive model" dialog should be disabled
    When user hovers over label of "Petal.Width" input in "Apply predictive model" dialog
    Then tooltip should contain text "'col 1' is also used for 'Petal.Length'. Choose a different column."
    When user selects "Petal.Length" in "Petal.Length" input in "Apply predictive model" dialog
    Then editor of "Petal.Width" input in "Apply predictive model" dialog should be valid
    And OK button in "Apply predictive model" dialog should be enabled
    When user clicks on OK button in "Apply predictive model" dialog
    Then the "Apply predictive model" dialog should close
    And the task bar entry that showed "Predicting Species" should have ended
    And the task bar should have shown 1 entry "Predicting Species"
    And a new column "Species (predicted)" should have been added
    And the table should have the columns "col 1, Sepal.Length, sepal width, Petal.Length, Species, Species (predicted)"
    And the model "forge-test-apply-iris-{time}" should have 3 applications on the server
    And the newest application of the model "forge-test-apply-iris-{time}" should have:
      | table_name  | forge-test-iris2-{time} |
      | column_name | Species (predicted)     |
      | status      | completed               |
    When user closes the current view
    Then no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Rows with missing values are skipped and counted by default; Impute predicts every row (3.13)
    Given user counts the task bar entries
    And user opens the demo file "iris.csv" as "forge-test-irisgaps-{time}"
    When user clears "Sepal.Width" column in rows "3, 10"
    And user clears "Species" column in rows "20"
    And user picks "ML > Forge > Apply..." from the top menu
    Then there should be 1 visible "Apply predictive model" dialog
    And Model input in "Apply predictive model" dialog should have value "forge-test-apply-iris-{time}"
    And summary of "Columns" group should have the text "4 of 4 matched"
    And summary of "Columns" group should be valid
    And "Columns" group should be collapsed
    And "Missing values" input in "Apply predictive model" dialog should be visible
    And "Missing values" input in "Apply predictive model" dialog should have "Skip rows" chosen
    And "Neighbors" input in "Apply predictive model" dialog should be hidden
    And OK button in "Apply predictive model" dialog should be enabled
    When user expands "Columns" group
    And user hovers over label of "Sepal.Width" input in "Apply predictive model" dialog
    Then tooltip should contain text "'Sepal.Width' has 2 missing values."
    When user hovers over label of "Missing values" input in "Apply predictive model" dialog
    Then tooltip should contain text "Rows with missing values in: Sepal.Width: 2."
    When user clicks on OK button in "Apply predictive model" dialog
    Then the "Apply predictive model" dialog should close
    And the task bar entry that showed "Predicting Species" should have ended
    And the task bar should have shown 1 entry "Predicting Species"
    And a new column "Species (predicted)" should have been added
    And an info balloon containing 'Added the column "Species (predicted)" to forge-test-irisgaps-{time}; 2 rows skipped (missing values).' should have been shown
    And the "Species (predicted)" column should be empty exactly in the rows where any of "Sepal.Length, Sepal.Width, Petal.Length, Petal.Width" columns is empty
    And the value of "Species (predicted)" column in row 3 should be ""
    And the value of "Species (predicted)" column in row 10 should be ""
    And the value of "Species (predicted)" column in row 20 should match "^(setosa|versicolor|virginica)$"
    And the model "forge-test-apply-iris-{time}" should have 4 applications on the server
    And the newest application of the model "forge-test-apply-iris-{time}" should have:
      | table_name   | forge-test-irisgaps-{time} |
      | row_count    | 150                        |
      | skipped_rows | 2                          |
      | column_name  | Species (predicted)        |
      | status       | completed                  |
    When user picks "ML > Forge > Apply..." from the top menu
    Then there should be 1 visible "Apply predictive model" dialog
    When user selects "Impute" in "Missing values" input in "Apply predictive model" dialog
    Then "Neighbors" input in "Apply predictive model" dialog should be visible
    And "Neighbors" input in "Apply predictive model" dialog should have value "4"
    And "Distance" input in "Apply predictive model" dialog should be visible
    And "Distance" input in "Apply predictive model" dialog should have value "Euclidean"
    When user clicks on OK button in "Apply predictive model" dialog
    Then the "Apply predictive model" dialog should close
    And the task bar entry that showed "Predicting Species" should have ended
    And the task bar should have shown 2 entries "Predicting Species"
    And a new column "Species (predicted 2)" should have been added
    And an info balloon containing 'Added the column "Species (predicted 2)" to forge-test-irisgaps-{time}.' should have been shown
    And "Species (predicted 2)" column should have no missing values
    And "Sepal.Width" column should have missing values
    And the value of "Sepal.Width" column in row 3 should be ""
    And the value of "Sepal.Width" column in row 10 should be ""
    And the model "forge-test-apply-iris-{time}" should have 5 applications on the server
    And the newest application of the model "forge-test-apply-iris-{time}" should have:
      | table_name   | forge-test-irisgaps-{time} |
      | skipped_rows | 0                          |
      | column_name  | Species (predicted 2)      |
      | status       | completed                  |
    And the "Species (predicted 2)" column should hold the predictions of the Forge model "forge-test-apply-iris-{time}" on "Sepal.Length, Sepal.Width, Petal.Length, Petal.Width" imputed with 4 neighbors and "Euclidean" distance
    When user closes the current view
    Then no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: A model row without a model file is refused with the reason, and is no application (3.7)
    Given user switches to the "forge-test-iris-{time}" table view
    And a Forge model row "forge-test-apply-hand-{time}" without a model file is on the server
    When user picks "ML > Forge > Apply..." from the top menu
    Then there should be 1 visible "Apply predictive model" dialog
    When user selects "forge-test-apply-hand-{time}" in Model input in "Apply predictive model" dialog
    Then "Apply predictive model" dialog should contain text "The model 'forge-test-apply-hand-{time}' has no model file, so it cannot be applied. It was not trained and saved in Forge."
    And Model input in "Apply predictive model" dialog should be invalid
    And "Columns" group in "Apply predictive model" dialog should be absent
    And OK button in "Apply predictive model" dialog should be disabled
    When user hovers over label of Model input in "Apply predictive model" dialog
    Then tooltip should contain text "The model 'forge-test-apply-hand-{time}' has no model file, so it cannot be applied."
    When user clicks on CANCEL button in "Apply predictive model" dialog
    Then the "Apply predictive model" dialog should close
    And the model "forge-test-apply-hand-{time}" should have 0 applications on the server
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: A script applies a model through Forge:applyModel; a table without the features is refused (3.8)
    Given user counts the task bar entries
    And user switches to the "forge-test-iris-{time}" table view
    When user applies the Forge model "forge-test-apply-iris-{time}" to the current table through the API
    Then the task bar entry that showed "Predicting Species" should have ended
    And the task bar should have shown 1 entry "Predicting Species"
    And the table should have a column "Species (predicted 3)"
    And the newest application of the model "forge-test-apply-iris-{time}" should have:
      | table_name  | forge-test-iris-{time} |
      | column_name | Species (predicted 3)  |
      | source      | api                    |
      | status      | completed              |
    Given user opens the demo file "demog.csv" as "forge-test-demog-{time}"
    Then applying the Forge model "forge-test-apply-iris-{time}" to the table "forge-test-demog-{time}" through the API should fail with "The table has no columns 'Sepal.Length', 'Sepal.Width', 'Petal.Length', 'Petal.Width' the model needs. Map them in columnNamesMap."
    And the table "forge-test-demog-{time}" should have no column matching "\(predicted"
    And the model "forge-test-apply-iris-{time}" should have 7 applications on the server
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Every problem of the mapping marks its row with a plain text and keeps OK disabled (3.14)
    Given user switches to the "forge-test-demog-{time}" table view
    When user picks "ML > Forge > Apply..." from the top menu
    Then there should be 1 visible "Apply predictive model" dialog
    When user selects "forge-test-apply-iris-{time}" in Model input in "Apply predictive model" dialog
    Then summary of "Columns" group should have the text "0 of 4 matched"
    And summary of "Columns" group should be invalid
    And "Columns" group should be expanded
    And editor of "Sepal.Length" input in "Apply predictive model" dialog should have text ""
    And editor of "Sepal.Width" input in "Apply predictive model" dialog should have text ""
    And editor of "Petal.Length" input in "Apply predictive model" dialog should have text ""
    And editor of "Petal.Width" input in "Apply predictive model" dialog should have text ""
    And editor of "Sepal.Length" input in "Apply predictive model" dialog should be invalid
    And editor of "Sepal.Width" input in "Apply predictive model" dialog should be invalid
    And editor of "Petal.Length" input in "Apply predictive model" dialog should be invalid
    And editor of "Petal.Width" input in "Apply predictive model" dialog should be invalid
    And "Missing values" input in "Apply predictive model" dialog should be hidden
    And OK button in "Apply predictive model" dialog should be disabled
    When user hovers over label of "Sepal.Length" input in "Apply predictive model" dialog
    Then tooltip should contain text "Choose a column for 'Sepal.Length'."
    And the column picker of "Sepal.Length" input in "Apply predictive model" dialog should offer exactly "AGE, HEIGHT, WEIGHT"
    When user selects "AGE" in "Sepal.Length" input in "Apply predictive model" dialog
    Then summary of "Columns" group should have the text "1 of 4 matched"
    When user selects "HEIGHT" in "Sepal.Width" input in "Apply predictive model" dialog
    Then summary of "Columns" group should have the text "2 of 4 matched"
    When user selects "WEIGHT" in "Petal.Length" input in "Apply predictive model" dialog
    Then summary of "Columns" group should have the text "3 of 4 matched"
    When user selects "AGE" in "Petal.Width" input in "Apply predictive model" dialog
    Then summary of "Columns" group should have the text "3 of 4 matched"
    And summary of "Columns" group should be invalid
    And editor of "Sepal.Length" input in "Apply predictive model" dialog should be valid
    And editor of "Petal.Width" input in "Apply predictive model" dialog should be invalid
    And OK button in "Apply predictive model" dialog should be disabled
    And "Missing values" input in "Apply predictive model" dialog should be visible
    And "Missing values" input in "Apply predictive model" dialog should have "Skip rows" chosen
    When user hovers over label of "Petal.Width" input in "Apply predictive model" dialog
    Then tooltip should contain text "'AGE' is also used for 'Sepal.Length'. Choose a different column."
    When user hovers over label of "Missing values" input in "Apply predictive model" dialog
    Then tooltip should contain text "HEIGHT: "
    And tooltip should contain text "WEIGHT: "
    When user removes "WEIGHT" column
    Then editor of "Petal.Length" input in "Apply predictive model" dialog should be invalid
    And summary of "Columns" group should have the text "2 of 4 matched"
    When user hovers over label of "Petal.Length" input in "Apply predictive model" dialog
    Then tooltip should contain text "The column 'WEIGHT' is no longer in the table. Choose another column for 'Petal.Length'."
    When user clicks on CANCEL button in "Apply predictive model" dialog
    Then the "Apply predictive model" dialog should close
    When user closes the current view
    Then no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: The old tool's menu and dialog know nothing of Forge; the column has Forge's tag and pane (3.9)
    Given user switches to the "forge-test-iris-{time}" table view
    Then the "ML > Models" top menu group should list exactly "Train Model..., Apply Model..."
    And the "ML > Forge" top menu group should list exactly "Models, Train..., Apply..."
    When user picks "ML > Models > Apply Model..." from the top menu
    Then there should be 1 visible "Apply predictive model" dialog
    And Inputs input in "Apply predictive model" dialog should be visible
    And "Columns" group in "Apply predictive model" dialog should be absent
    And Model input in "Apply predictive model" dialog should not contain text "forge-test-"
    When user clicks on CANCEL button in "Apply predictive model" dialog
    Then the "Apply predictive model" dialog should close
    And the "Species (predicted 2)" column should be tagged with the id of the Forge model "forge-test-apply-iris-{time}"
    And the "Species (predicted 2)" column should not have the tag "predictive.model"
    Given the "Species (predicted 2)" column is the current object
    Then the context panel title should be "Species (predicted 2)"
    When user expands "Predicted by" pane in context panel
    Then "Predicted by" pane in context panel should contain text "forge-test-apply-iris-{time}"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Deleted with the trash icon, the model is gone and its prediction column says so (4.10)
    When user picks "ML > Forge > Models" from the top menu
    Then the "Forge" view should be current
    When user remembers the catalog loads
    And user clicks on "Refresh" icon in forge app
    Then the shown catalog should have reloaded
    When user clicks on the Forge model "forge-test-apply-hand-{time}" in second grid viewer
    Then the context panel should show the Forge model "forge-test-apply-hand-{time}"
    When user clicks on the Forge model "forge-test-apply-iris-{time}" in second grid viewer
    Then the context panel should show the Forge model "forge-test-apply-iris-{time}"
    When user remembers the catalog loads
    And user clicks on "Delete model" icon in forge app
    Then there should be 1 visible "Delete model" dialog
    And "Delete model" dialog should contain text 'Delete the model "forge-test-apply-iris-{time}"?'
    When user clicks on OK button in "Delete model" dialog
    Then every catalog should have reloaded
    And second grid viewer should not list the Forge model "forge-test-apply-iris-{time}"
    And the context panel should not show the Forge model "forge-test-apply-iris-{time}"
    And 0 Forge models named "forge-test-apply-iris-{time}" should be on the server
    And the model "forge-test-apply-iris-{time}" should have 0 applications on the server
    Given user switches to the "forge-test-iris-{time}" table view
    And the "Species (predicted)" column is the current object
    Then the context panel title should be "Species (predicted)"
    When user expands "Predicted by" pane in context panel
    Then "Predicted by" pane in context panel should contain text "The model that predicted this column is no longer available to you."
    And no error or warning balloon should have been shown
    And no errors should have been logged
