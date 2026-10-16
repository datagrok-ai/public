@journey @serial
Feature: Applying the quality regressor to winequality, and cancelling an application
  ML | Forge | Apply... over winequality opened from the demo files, with a quality regressor the Background trains and
  saves through Forge:trainModel (XGBoost, the method the manual chooses; no data stored). Translated from the manual
  scenarios 3.2, 3.6 and 4.4 (step 3), with 3.5's records of these applications, of
  core/docs/features/platform/predictive-modeling/forge/manual-checks.md. The predictions are claimed against the
  training: the mean absolute difference of quality and its prediction, computed over the table (an empty prediction
  fails it), equals the stored train MAE. The cancel is made reliable by a batch of one row: the application becomes
  thousands of steps, and the cross of its task bar entry is clicked once the entry is there; the entry's first
  progress text shows that it ran in batches of one row. That it stopped is the application's status, cancelled, and no
  prediction column; the entry never showing the last batch's text only adds to it (a text shown between two reads of
  the watch would be missed). An application's end is its task bar entry
  gone (followed as an element: the entry shows the progress instead of "Predicting quality" while it runs). Server
  claims read the stand through grok s, scoped to this feature's model.

  Not claimed here: 3.2 step 1 as a separate read (the MAE claim reads the stored metric itself) and "the iris models
  follow in the list" (the order of models that do not fit is the stand's, not the feature's); the rows of the
  Columns group between fixed.acidity and alcohol; 4.4 step 4 as a separate query (the newest application is claimed
  instead), the Activity grid's When, and step 3's play icon (the feature opens the Forge view through ML | Forge |
  Models, the catalog's menu path, not the icon).

  Background:
    Given user is logged in
    And no Forge models named "forge-test-apply-quality-{time}" are on the server
    And no Forge training runs of the table "forge-test-wine-{time}" are on the server
    And user opens the demo file "winequality.csv" as "forge-test-wine-{time}"
    And a Forge model "forge-test-apply-quality-{time}" is trained on the current table through the API with:
      | target | quality |
      | method | XGBoost |

  Scenario: The quality model is preset and matches every column; its predictions are whole numbers with the training's MAE (3.2, 3.5)
    Given user counts the task bar entries
    When user picks "ML > Forge > Apply..." from the top menu
    Then there should be 1 visible "Apply predictive model" dialog
    And Table input in "Apply predictive model" dialog should have value "forge-test-wine-{time}"
    And Model input in "Apply predictive model" dialog should have value "forge-test-apply-quality-{time}"
    And summary of "Columns" group should have the text "11 of 11 matched"
    And summary of "Columns" group should be valid
    And "Columns" group should be collapsed
    And "Missing values" input in "Apply predictive model" dialog should be hidden
    And "More options" group should be collapsed
    And OK button in "Apply predictive model" dialog should be enabled
    When user expands "Columns" group
    Then editor of "fixed.acidity" input in "Apply predictive model" dialog should have text "fixed.acidity"
    And editor of "alcohol" input in "Apply predictive model" dialog should have text "alcohol"
    When user expands "More options" group
    Then "Batch size" input in "Apply predictive model" dialog should have value "10000"
    When user clicks on OK button in "Apply predictive model" dialog
    Then the "Apply predictive model" dialog should close
    And the task bar entry that showed "Predicting quality" should have ended
    And the task bar should have shown 1 entry "Predicting quality"
    And a new column "quality (predicted)" should have been added
    And "quality (predicted)" column should have no missing values
    And an info balloon containing 'Added the column "quality (predicted)" to forge-test-wine-{time}.' should have been shown
    And every value of "quality (predicted)" column should match "^-?\d+$"
    And the mean absolute difference of "quality" and "quality (predicted)" columns of the table "forge-test-wine-{time}" should equal the Train MAE of the Forge model "forge-test-apply-quality-{time}" within 0.001
    And the model "forge-test-apply-quality-{time}" should have 1 application on the server
    And the newest application of the model "forge-test-apply-quality-{time}" should have:
      | table_name   | forge-test-wine-{time} |
      | row_count    | 6497                   |
      | skipped_rows | 0                      |
      | column_name  | quality (predicted)    |
      | source       | ui                     |
      | status       | completed              |
      | duration_ms  | number                 |
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: A cancelled application adds no column and is recorded; the model's Activity and Last run show it (3.6, 4.4)
    Given user counts the task bar entries
    When user picks "ML > Forge > Apply..." from the top menu
    Then there should be 1 visible "Apply predictive model" dialog
    And Model input in "Apply predictive model" dialog should have value "forge-test-apply-quality-{time}"
    When user expands "More options" group
    And user enters "1" into "Batch size" input in "Apply predictive model" dialog
    And user clicks on OK button in "Apply predictive model" dialog
    Then the "Apply predictive model" dialog should close
    When user cancels the "Predicting quality" task bar entry
    Then a warning balloon containing "Application was cancelled." should have been shown
    And the task bar entry that showed "Predicting quality" should have ended
    And the task bar should have shown 1 entry "Predicting quality"
    And the task bar entry that showed "Predicting quality" should have shown "Rows 1-1 of 6497"
    And the task bar entry that showed "Predicting quality" should not have shown "Rows 6497-6497 of 6497"
    And the table "forge-test-wine-{time}" should have no column matching "^quality \(predicted \d+\)$"
    And the newest application of the model "forge-test-apply-quality-{time}" should have:
      | table_name  | forge-test-wine-{time} |
      | status      | cancelled              |
      | column_name | empty                  |
    And the model "forge-test-apply-quality-{time}" should have 2 applications on the server
    When user picks "ML > Forge > Models" from the top menu
    Then the "Forge" view should be current
    Given the context panel is open
    When user clicks on the Forge model "forge-test-apply-quality-{time}" in second grid viewer
    Then the context panel should show the Forge model "forge-test-apply-quality-{time}"
    And value of "Applications" detail in context panel should have the text "2"
    And value of "Last run" detail in context panel should contain text "(cancelled)"
    And value of "Last run" detail in context panel should start with the time of the newest application of the model "forge-test-apply-quality-{time}"
    When user expands "Activity" pane in context panel
    Then "Activity" pane in context panel should contain text "2 applications"
    When user scrolls grid viewer in "Activity" pane in context panel to the "Prediction column" column
    Then the "text of cell 1 of Prediction column" reading of grid viewer in "Activity" pane in context panel should be ""
    And the "text of cell 2 of Prediction column" reading of grid viewer in "Activity" pane in context panel should be "quality (predicted)"
    When user scrolls grid viewer in "Activity" pane in context panel to the "Status" column
    Then the "text of cell 1 of Status" reading of grid viewer in "Activity" pane in context panel should be "cancelled"
    And the "text of cell 2 of Status" reading of grid viewer in "Activity" pane in context panel should be "completed"
    And no error or warning balloon should have been shown
    And no errors should have been logged
