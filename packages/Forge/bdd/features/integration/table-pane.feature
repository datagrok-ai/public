@journey @serial
Feature: The table pane Forge models: the models that fit a table, their tooltip, menu and run icon, and Train
  The context panel of a table holds Forge's pane Forge models after the built-in panes: one line per saved model
  that fits the table, newest first, then Train. Translated from the manual scenarios 6.1 (steps 1-7) and 4.13 (step
  3) of core/docs/features/platform/predictive-modeling/forge/manual-checks.md. The Background trains four models
  through Forge:trainModel (no data stored): two iris classifiers (XGBoost; the second, saved after the first, has the
  name that sorts first, so "newest first" is not the names' order), a wine-quality regressor and a regressor of
  demog's AGE by HEIGHT and WEIGHT (Linear Regression), each on its own open table. The user's models fit the same
  tables, so the lines are claimed among this run's own (the forge-test- names ending in the run's time), in their
  order. The pane is built anew each time it is shown, after the models are read, so every claim about it polls; the
  run icon's application is counted on the server for its model. Train is clicked in the pane of the demog table
  while the iris view is current, so the Train view opens on the pane's table, not on the current view's.

  Not claimed here: the corner click of steps 1 and 4 (a click on the grid's top-left corner does not make the table
  current again, a platform issue: a step sets the table as the current object instead); step 4's cell click back to
  the table (a click on a cell of the iris grid makes the cell the current object, not the table: the same step sets
  the table again); "newest first" beyond the two iris models; the run icon's colour and "no model icon" (visual);
  the built-in Models pane's content (the old tool's list; the pane is claimed present); the values of the tooltip's
  Created by, Created and Updated (the account's name and the times of the run; the Created by caption is claimed); the
  menu items Edit..., Clone, Delete and Watch
  of step 3 (the row menu of a line lists Open, Share..., History, Copy link, Apply... and Download); step 8 (the
  package test UI: Table panel "tableModelsPanel is a panel of a table"). "No model fits this table." reads every
  model the account sees: a model of the user's whose features all have a close column in country-languages
  (Country, Languages) would make it fail.

  Where the product moved on from the manual text: the corner click is replaced by the table set as the current
  object; the new column of step 5 is claimed by its name, not as "a new column" (that count is armed by a top-menu
  command only). Step 6's cars is not a table no model fits: its height and curb.weight are close columns of HEIGHT
  and WEIGHT, so the demog regressor is listed there (claimed), and "No model fits this table." is claimed on
  country-languages instead.

  Background:
    Given user is logged in
    And no Forge models named "forge-test-pane-iris-{time}, forge-test-pane-flower-{time}, forge-test-pane-quality-{time}, forge-test-pane-demog-{time}" are on the server
    And no Forge training runs of the tables "forge-test-iris-{time}, forge-test-wine-{time}, forge-test-demog-{time}" are on the server
    And user opens the demo file "winequality.csv" as "forge-test-wine-{time}"
    And a Forge model "forge-test-pane-quality-{time}" is trained on the current table through the API with:
      | target | quality           |
      | method | Linear Regression |
    And user opens the demo file "demog.csv" as "forge-test-demog-{time}"
    And a Forge model "forge-test-pane-demog-{time}" is trained on the current table through the API with:
      | target   | AGE               |
      | features | HEIGHT, WEIGHT    |
      | method   | Linear Regression |
    And user opens the demo file "iris.csv" as "forge-test-iris-{time}"
    And a Forge model "forge-test-pane-iris-{time}" is trained on the current table through the API with:
      | target | Species |
      | method | XGBoost |
    And a Forge model "forge-test-pane-flower-{time}" is trained on the current table through the API with:
      | target | Species |
      | method | XGBoost |
    And the context panel is open
    Then no balloon of any colour should have been shown
    And no errors should have been logged

  Scenario: A table's context panel lists the models that fit it, newest first, after the built-in panes (6.1, 4.13)
    Given the table "forge-test-iris-{time}" is the current object
    Then "Models" pane in context panel should be visible
    And the "Forge models" pane in context panel should come after the "Script" pane
    When user expands "Forge models" pane in context panel
    Then the model lines of "Forge models" pane in context panel should list exactly "forge-test-pane-flower-{time}, forge-test-pane-iris-{time}" among those ending in "-{time}"
    And "forge-test-pane-quality-{time}" model line in "Forge models" pane in context panel should be absent
    And "forge-test-pane-demog-{time}" model line in "Forge models" pane in context panel should be absent
    And menu arrow of "forge-test-pane-iris-{time}" model line in "Forge models" pane in context panel should be visible
    And "Train" button in "Forge models" pane in context panel should be visible
    When user hovers over run icon of "forge-test-pane-iris-{time}" model line in "Forge models" pane in context panel
    Then tooltip should contain text "Run prediction"
    When user hovers over "Train" button in "Forge models" pane in context panel
    Then tooltip should contain text "Train a new model on this table."
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: A line's name shows the model's tooltip and the row's menu (6.1)
    When user hovers over name of "forge-test-pane-iris-{time}" model line in "Forge models" pane in context panel
    Then tooltip should contain text "forge-test-pane-iris-{time}"
    And tooltip should contain text "Created by"
    And value of "Trained on" detail in tooltip should have the text "forge-test-iris-{time} (150 rows)"
    And value of "Last run" detail in tooltip should have the text "Never"
    And value of "Tags" detail in tooltip should have the text ""
    And value of "Inputs" detail in tooltip should have the text "Sepal.Length, Sepal.Width, Petal.Length, Petal.Width"
    And value of "Outputs" detail in tooltip should have the text "Species"
    And value of "Method" detail in tooltip should have the text "XGBoost"
    And value of "Applicable to" detail in tooltip should have the text "forge-test-iris-{time}"
    When user right-clicks on name of "forge-test-pane-iris-{time}" model line in "Forge models" pane in context panel
    Then the open menu should list "Open"
    And the open menu should list "Share..."
    And the open menu should list "History"
    And the open menu should list "Copy link"
    And the open menu should list "Apply..."
    And the open menu should list "Download"
    When user presses Escape
    Then context menu should be absent
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: The name opens the model in the context panel (6.1)
    When user clicks on name of "forge-test-pane-iris-{time}" model line in "Forge models" pane in context panel
    Then the context panel should show the Forge model "forge-test-pane-iris-{time}"
    And "Details" pane in context panel should be visible
    Given the table "forge-test-iris-{time}" is the current object
    Then "Forge models" pane in context panel should be visible
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: The run icon applies the model to the table (6.1)
    Given user counts the task bar entries
    When user expands "Forge models" pane in context panel
    And user clicks on run icon of "forge-test-pane-iris-{time}" model line in "Forge models" pane in context panel
    Then there should be 1 visible "Apply predictive model" dialog
    And Table input in "Apply predictive model" dialog should have value "forge-test-iris-{time}"
    And Model input in "Apply predictive model" dialog should have value "forge-test-pane-iris-{time}"
    And summary of "Columns" group should have the text "4 of 4 matched"
    When user clicks on OK button in "Apply predictive model" dialog
    Then the "Apply predictive model" dialog should close
    And the task bar entry that showed "Predicting Species" should have ended
    And the table should have a column "Species (predicted)"
    And an info balloon containing 'Added the column "Species (predicted)" to forge-test-iris-{time}.' should have been shown
    And the model "forge-test-pane-iris-{time}" should have 1 application on the server
    And the newest application of the model "forge-test-pane-iris-{time}" should have:
      | source      | ui                  |
      | status      | completed           |
      | column_name | Species (predicted) |
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Each table lists only the models that fit it; a table no model fits says so (6.1)
    When user switches to the "forge-test-wine-{time}" table view
    Given the table "forge-test-wine-{time}" is the current object
    When user expands "Forge models" pane in context panel
    Then the model lines of "Forge models" pane in context panel should list exactly "forge-test-pane-quality-{time}" among those ending in "-{time}"
    When user switches to the "forge-test-demog-{time}" table view
    Given the table "forge-test-demog-{time}" is the current object
    When user expands "Forge models" pane in context panel
    Then the model lines of "Forge models" pane in context panel should list exactly "forge-test-pane-demog-{time}" among those ending in "-{time}"
    Given user opens the demo file "cars.csv" as "forge-test-cars-{time}"
    And the table "forge-test-cars-{time}" is the current object
    When user expands "Forge models" pane in context panel
    Then the model lines of "Forge models" pane in context panel should list exactly "forge-test-pane-demog-{time}" among those ending in "-{time}"
    Given user opens the demo file "country-languages.csv" as "forge-test-langs-{time}"
    And the table "forge-test-langs-{time}" is the current object
    When user expands "Forge models" pane in context panel
    Then "Forge models" pane in context panel should contain text "No model fits this table."
    And the model lines of "Forge models" pane in context panel should list exactly "" among those ending in "-{time}"
    And "Train" button in "Forge models" pane in context panel should be visible
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Train in the pane opens the Train view on the pane's table, not the current view's (6.1)
    When user switches to the "forge-test-iris-{time}" table view
    Given the table "forge-test-demog-{time}" is the current object
    And user remembers the shown training
    When user expands "Forge models" pane in context panel
    And user clicks on "Train" button in "Forge models" pane in context panel
    Then the "Predictive model" view should be current
    And Table input should have value "forge-test-demog-{time}"
    And a new training's result should be shown
    And the shown training should be the newest training run of the table "forge-test-demog-{time}", with:
      | status      | completed |
      | engine_name | XGBoost   |
      | target_name | SEVERITY  |
      | source      | ui        |
    When user closes the current view
    Then no error or warning balloon should have been shown
    And no errors should have been logged
