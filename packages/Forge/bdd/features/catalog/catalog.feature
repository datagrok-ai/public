@journey @serial
Feature: The model catalog: a model's panes, its commands, comparing, sharing, evaluating and deleting
  ML | Forge | Models with three iris and quality models the Background trains and saves through Forge:trainModel
  (XGBoost; the iris classifier tagged "demo", a second iris classifier keeping a reference to iris.csv, the quality
  regressor saved last, so it is the newest), and a model row made by hand without a model file. Translated from the
  manual scenarios 2.6, 3.7 (step 3), 4.0 (step 3), 4.1, 4.2 (steps 1-2), 4.3, 4.4 (steps 1-2), 4.5, 4.6, 4.7 (steps
  2-3), 4.8, 4.9, 4.11, 4.12 (steps 1-3) and 5.15 of core/docs/features/platform/predictive-modeling/forge/manual-checks.md.
  A reload of the catalog is claimed by the load it shows (data-forge-catalog-load on the models pane of every open
  Forge view, a hidden one too): a delete or an edit is followed by a later load in every open catalog before the grid
  is read. Every claim names this feature's models; the stand's other models are listed too and never named. The
  context panel is rebuilt only when another row is clicked (an application of the shown model is followed in place
  by its Details and Activity), so a scenario that needs a model's panel built anew reopens the model (another row
  shown first). The Performance pane is claimed against the model row's stored metrics, and Last run
  against the newest application's time; an application's or an evaluation's end is its task bar entry gone, and
  each scenario counts the entries it caused.

  Not claimed here: 4.1's title (the star and the context actions; the package's UI: Catalog tests), Created and
  Updated, and the Data source of the Reference model; 4.2's Domains card and card mode (manual), the tooltip's
  "Created <date>" line and "another cell keeps the grid's own tooltip"; 4.3's other header and cell tooltips and the
  seed's copy icon; 4.4's When (the time) beyond Last run; 4.5's "no card cut, no scroll bar" (visual), the
  comparison's title icon and card fields, "Compare opens the same view" (the Compare models view is claimed by its
  table, not against the panel), and "the values equal the Performance panes" (claimed as ranges and empty cells: the
  comparison's float columns have no format of three decimals); 4.6's Domains row menu (manual); 4.7 step 1
  (train-iris.feature), step 3's Description input, step 4 (the Domains view, manual) and the browser's autofill
  (headless never shows it); 4.8's
  "a few dozen KB, the size of the folder" (the file is claimed by its size against the server's model file) and the
  Domains row menu (manual); 4.9 step 2 (the Domains view, manual) and "a comparison that included it compares the
  rest" (no comparison is open when the model is deleted); 4.9's and 2.6's "the panel shows the current row's
  model or nothing" (the panel is claimed to have left the deleted model); 4.11's "every model" beyond this feature's
  three, and Applicable to's exact list (claimed by the table it holds and the one it lost: the platform's table input
  lists every open table twice, the shell's tables and its context's); 4.12 step 1's All users row, its "View and use" and the empty Can edit; 4.12 step 4 (the second user on the
  release build, manual) and the dialog's group "All users" (the feature shares with the second account
  DATAGROK_SHARING_LOGIN and takes the grant back: a deleted model's entity keeps its grants, so a grant to everyone
  would outlive the run); 5.15's tab icon, the Insights & Tips text, step 2's "the Target list offers the text columns
  only" and the divider drag; 2.6's Refresh tooltip. The reason the disabled OK of Edit model gives on hover (4.7 step
  3, "Enter a name for the model."): the platform shows no tooltip of a disabled button inside a dialog
  (train-iris.feature's @known-failure scenario); OK is claimed disabled.

  Where the product moved on from the manual text: ML | Forge | Models focuses an open catalog, so 2.6's second
  catalog is the app opened again; 4.4's "click the catalog row again" is not needed: the shown model's Details and
  Activity show the application by themselves.
  4.6 opens a second iris copy, so the iris model fits two tables, and then the catalog again: the demo opener names
  the copy after the platform has listed it ("iris (2)"), and the open catalog's Applicable to keeps a table under the
  name it was listed with until a table closes. 4.6's "make iris the current table (click its view, then back to the
  catalog)": no table is current while the catalog is in front (it is not a table view), so with Applicable to empty
  Apply... opens on the first open table the model fits, the first iris copy, though the second one's view was shown
  last.

  Background:
    Given user is logged in
    And no Forge models named "forge-test-cat-iris-{time}, forge-test-cat-ref-{time}, forge-test-cat-quality-{time}, forge-test-cat-hand-{time}, forge-test-cat-edited-{time}" are on the server
    And no Forge training runs of the tables "forge-test-iris-{time}, forge-test-wine-{time}" are on the server
    And user opens the demo file "iris.csv" as "forge-test-iris-{time}"
    And a Forge model "forge-test-cat-ref-{time}" is trained on the current table through the API with:
      | target  | Species   |
      | method  | XGBoost   |
      | storage | reference |
    And a Forge model "forge-test-cat-iris-{time}" is trained on the current table through the API with:
      | target | Species |
      | method | XGBoost |
      | tags   | demo    |
    And user opens the demo file "winequality.csv" as "forge-test-wine-{time}"
    And a Forge model "forge-test-cat-quality-{time}" is trained on the current table through the API with:
      | target | quality |
      | method | XGBoost |
    When user picks "ML > Forge > Models" from the top menu
    Then the "Forge" view should be current
    Given the context panel is open
    Then second grid viewer should list the Forge model "forge-test-cat-iris-{time}"
    And second grid viewer should list the Forge model "forge-test-cat-ref-{time}"
    And second grid viewer should list the Forge model "forge-test-cat-quality-{time}"
    And no balloon of any colour should have been shown
    And no errors should have been logged

  Scenario: A saved model is a platform entity with its author's grants; its Details (4.0, 4.1)
    Then the Forge model "forge-test-cat-iris-{time}" should be an entity on the server with the author's grants
    When user clicks on the Forge model "forge-test-cat-iris-{time}" in second grid viewer
    Then the context panel should show the Forge model "forge-test-cat-iris-{time}"
    And "Details" pane in context panel should be expanded
    And "Performance" pane in context panel should be visible
    And "Activity" pane in context panel should be visible
    And "Sharing" pane in context panel should be visible
    And "History" pane in context panel should be visible
    And the context panel should not contain text "Chats"
    And value of "Author" detail in context panel should contain text "admin"
    And value of "Table" detail in context panel should have the text "forge-test-iris-{time} (150 rows)"
    And value of "Data storage" detail in context panel should have the text "None"
    And value of "Last run" detail in context panel should have the text "Never"
    And value of "Applications" detail in context panel should have the text "0"
    And value of "Features" detail in context panel should have the text "Sepal.Length, Sepal.Width, Petal.Length, Petal.Width"
    And value of "Target" detail in context panel should have the text "Species"
    And value of "Method" detail in context panel should have the text "XGBoost"
    And value of "Task" detail in context panel should have the text "classification"
    And value of "Applicable to" detail in context panel should have the text "forge-test-iris-{time}"
    And there should be 1 visible tag chips in "Details" pane in context panel
    When user switches to the "forge-test-iris-{time}" table view
    And user closes the current view
    And user switches to the "Forge" view
    And user clicks on the Forge model "forge-test-cat-quality-{time}" in second grid viewer
    Then the context panel should show the Forge model "forge-test-cat-quality-{time}"
    And value of "Table" detail in context panel should have the text "forge-test-wine-{time} (6497 rows)"
    And value of "Task" detail in context panel should have the text "regression"
    And value of "Applicable to" detail in context panel should have the text "forge-test-wine-{time}"
    When user clicks on the Forge model "forge-test-cat-iris-{time}" in second grid viewer
    Then the context panel should show the Forge model "forge-test-cat-iris-{time}"
    And value of "Method" detail in context panel should have the text "XGBoost"
    And "Details" pane in context panel should not contain text "Applicable to"
    Given user opens the demo file "iris.csv" as "forge-test-iris-{time}"
    And user switches to the "Forge" view
    Then no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: A model's tooltip in the catalog says what it predicts, from what, how, and whose it is (4.2)
    When user hovers over the Forge model "forge-test-cat-iris-{time}" in second grid viewer
    Then tooltip should contain text "forge-test-cat-iris-{time}"
    And tooltip should contain text "Predicts Species from 4 features with XGBoost · classification · 150 rows"
    And tooltip should contain text "Sepal.Length, Sepal.Width, Petal.Length, Petal.Width"
    And tooltip should contain text "by admin"
    And tooltip should contain text "Data storage: None · Tags: demo"
    When user hovers over the Forge model "forge-test-cat-quality-{time}" in second grid viewer
    Then tooltip should contain text "Predicts quality from 11 features with XGBoost · regression · 6497 rows"
    And tooltip should contain text "fixed.acidity, volatile.acidity, citric.acid, residual.sugar, chlorides +6 more"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Performance shows the stored metrics and offers Run Evaluation; a model without metrics says so (4.3)
    When user reopens the Forge model "forge-test-cat-iris-{time}" in second grid viewer
    Then the context panel should show the Forge model "forge-test-cat-iris-{time}"
    When user expands "Performance" pane in context panel
    Then the "text of cell 1 of Metric" reading of grid viewer in "Performance" pane in context panel should be "Accuracy"
    And the "text of cell 2 of Metric" reading of grid viewer in "Performance" pane in context panel should be "F1"
    And the "rows" reading of grid viewer in "Performance" pane in context panel should be 2
    And the "text of cell 1 of Train" reading of grid viewer in "Performance" pane in context panel should be the "metrics.train.accuracy" of the Forge model "forge-test-cat-iris-{time}" on the server to three decimals
    And the "text of cell 1 of Validation" reading of grid viewer in "Performance" pane in context panel should be the "metrics.validation.accuracy" of the Forge model "forge-test-cat-iris-{time}" on the server to three decimals
    And the "text of cell 2 of Train" reading of grid viewer in "Performance" pane in context panel should be the "metrics.train.f1" of the Forge model "forge-test-cat-iris-{time}" on the server to three decimals
    And the "text of cell 2 of Validation" reading of grid viewer in "Performance" pane in context panel should be the "metrics.validation.f1" of the Forge model "forge-test-cat-iris-{time}" on the server to three decimals
    And "Performance" pane in context panel should contain text "Validation: 5-fold cross-validation on 150 rows"
    And the seed in the context panel should be the seed of the Forge model "forge-test-cat-iris-{time}" on the server
    When user hovers over the "header Validation" area of grid viewer in "Performance" pane in context panel
    Then tooltip should contain text "Value on rows the model did not see"
    And "Run Evaluation" button in "Performance" pane in context panel should be visible
    When user hovers over "Run Evaluation" button in "Performance" pane in context panel
    Then tooltip should contain text "Apply the model to its training data, or to an open table, and show how well it predicts."
    When user clicks on the Forge model "forge-test-cat-quality-{time}" in second grid viewer
    Then the context panel should show the Forge model "forge-test-cat-quality-{time}"
    When user expands "Performance" pane in context panel
    Then the "text of cell 1 of Metric" reading of grid viewer in "Performance" pane in context panel should be "MSE"
    And the "text of cell 2 of Metric" reading of grid viewer in "Performance" pane in context panel should be "RMSE"
    And the "text of cell 3 of Metric" reading of grid viewer in "Performance" pane in context panel should be "MAE"
    And the "text of cell 4 of Metric" reading of grid viewer in "Performance" pane in context panel should be "R2"
    Given a Forge model row "forge-test-cat-hand-{time}" without a model file is on the server
    When user remembers the catalog loads
    And user clicks on "Refresh" icon in forge app
    Then every catalog should have reloaded
    When user clicks on the Forge model "forge-test-cat-hand-{time}" in second grid viewer
    Then the context panel should show the Forge model "forge-test-cat-hand-{time}"
    When user expands "Performance" pane in context panel
    Then "Performance" pane in context panel should contain text "No metrics were recorded for this model."
    And "Run Evaluation" button in "Performance" pane in context panel should be absent
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Applied from the catalog, the model's Activity, Last run and Applications show the application (4.4)
    Given user counts the task bar entries
    When user clicks on the Forge model "forge-test-cat-iris-{time}" in second grid viewer
    Then the context panel should show the Forge model "forge-test-cat-iris-{time}"
    And value of "Last run" detail in context panel should have the text "Never"
    And value of "Applications" detail in context panel should have the text "0"
    When user expands "Activity" pane in context panel
    Then "Activity" pane in context panel should contain text "Not applied yet."
    When user clicks on "Apply model" icon in forge app
    Then there should be 1 visible "Apply predictive model" dialog
    And Model input in "Apply predictive model" dialog should have value "forge-test-cat-iris-{time}"
    And Table input in "Apply predictive model" dialog should have value "forge-test-iris-{time}"
    When user clicks on OK button in "Apply predictive model" dialog
    Then the "Apply predictive model" dialog should close
    And the "forge-test-iris-{time}" view should be current
    And the task bar entry that showed "Predicting Species" should have ended
    And the task bar should have shown 1 entry "Predicting Species"
    And an info balloon containing 'Added the column "Species (predicted)" to forge-test-iris-{time}.' should have been shown
    And the table should have a column "Species (predicted)"
    And the model "forge-test-cat-iris-{time}" should have 1 application on the server
    Given user switches to the "Forge" view
    Then the context panel should show the Forge model "forge-test-cat-iris-{time}"
    And value of "Applications" detail in context panel should have the text "1"
    And value of "Last run" detail in context panel should start with the time of the newest application of the model "forge-test-cat-iris-{time}"
    When user expands "Activity" pane in context panel
    Then "Activity" pane in context panel should contain text "1 application"
    And the "text of cell 1 of Who" reading of grid viewer in "Activity" pane in context panel should be "admin"
    And the "text of cell 1 of Table" reading of grid viewer in "Activity" pane in context panel should be "forge-test-iris-{time}"
    And the "text of cell 1 of Rows" reading of grid viewer in "Activity" pane in context panel should be "150"
    When user scrolls grid viewer in "Activity" pane in context panel to the "Prediction column" column
    Then the "text of cell 1 of Prediction column" reading of grid viewer in "Activity" pane in context panel should be "Species (predicted)"
    When user scrolls grid viewer in "Activity" pane in context panel to the "Status" column
    Then the "text of cell 1 of Status" reading of grid viewer in "Activity" pane in context panel should be "completed"
    When user hovers over the "header Status" area of grid viewer in "Activity" pane in context panel
    Then tooltip should contain text "Outcome of the application."
    When user hovers over the "cell 1 of Status" area of grid viewer in "Activity" pane in context panel
    Then tooltip should contain text "Completed: the column was added."
    When user scrolls grid viewer in "Activity" pane in context panel to the "Source" column
    Then the "text of cell 1 of Source" reading of grid viewer in "Activity" pane in context panel should be "ui"
    When user hovers over the "cell 1 of Source" area of grid viewer in "Activity" pane in context panel
    Then tooltip should contain text "The Apply dialog or the catalog."
    When user scrolls grid viewer in "Activity" pane in context panel to the "Duration (ms)" column
    Then the "text of cell 1 of Duration (ms)" reading of grid viewer in "Activity" pane in context panel should be between 0 and 600000
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Two selected models compare in the context panel and in a view of their own (4.5)
    Then "Compare in a new view" icon in forge app should be disabled
    When user selects the Forge models "forge-test-cat-quality-{time}, forge-test-cat-iris-{time}" in second grid viewer
    Then "Compare in a new view" icon in forge app should be enabled
    And the context panel should contain text "Compare 2 models"
    And forms viewer in context panel should be visible
    When user hovers over "Compare in a new view" icon in forge app
    Then tooltip should contain text "Compare in a new view"
    When user clicks on "Compare in a new view" icon in forge app
    Then the "Compare models" view should be current
    And table "Compare models" should have columns "Name, Description, Method, Task, Target, Training rows, Created, MSE (train), MSE (validation), RMSE (train), RMSE (validation), MAE (train), MAE (validation), R2 (train), R2 (validation), Accuracy (train), Accuracy (validation), F1 (train), F1 (validation)"
    And table "Compare models" should have 2 rows
    And forms viewer should be visible
    And the "text of cell 1 of Name" reading of grid viewer should be "forge-test-cat-quality-{time}"
    And the "text of cell 2 of Name" reading of grid viewer should be "forge-test-cat-iris-{time}"
    When user scrolls the grid to the "R2 (validation)" column
    Then the "text of cell 1 of R2 (validation)" reading of grid viewer should be between 0.3 and 0.7
    When user scrolls the grid to the "MSE (train)" column
    Then the "text of cell 2 of MSE (train)" reading of grid viewer should be ""
    When user scrolls the grid to the "Accuracy (validation)" column
    Then the "text of cell 1 of Accuracy (validation)" reading of grid viewer should be ""
    And the "text of cell 2 of Accuracy (validation)" reading of grid viewer should be between 0.85 and 1
    When user closes the current view
    And user switches to the "Forge" view
    And user right-clicks on the Forge model "forge-test-cat-quality-{time}" in second grid viewer
    Then the open menu should list "Compare"
    When user closes the context menu
    And user toggles the selection of the Forge model "forge-test-cat-quality-{time}" in second grid viewer
    Then "Compare in a new view" icon in forge app should be disabled
    And the context panel should show the Forge model "forge-test-cat-quality-{time}"
    When user toggles the selection of the Forge model "forge-test-cat-iris-{time}" in second grid viewer
    Then "Compare in a new view" icon in forge app should be disabled
    And second grid viewer should have 0 selected rows
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: A model's row menu has Forge's commands; Apply... opens on the table Applicable to holds, else the first open one it fits (4.6)
    When user right-clicks on the Forge model "forge-test-cat-iris-{time}" in second grid viewer
    Then the open menu should list "Apply..."
    And the open menu should list "Edit model..."
    And the open menu should list "Download"
    And the open menu should list "Delete model"
    And the open menu should not list "Apply to"
    And the open menu should not list "Compare"
    When user closes the context menu
    Given user opens the demo file "iris.csv" as "forge-test-iris2-{time}"
    And user switches to the "Forge" view
    When user closes the current view
    And user switches to the "forge-test-wine-{time}" table view
    And user picks "ML > Forge > Models" from the top menu
    Then the "Forge" view should be current
    And there should be 1 open view named "Forge"
    When user selects "forge-test-iris2-{time}" in "Applicable to" input in forge app
    Then "Applicable to" input in forge app should have value "forge-test-iris2-{time}"
    When user right-clicks on the Forge model "forge-test-cat-iris-{time}" in second grid viewer
    And user picks "Apply..." from the open menu
    Then there should be 1 visible "Apply predictive model" dialog
    And Table input in "Apply predictive model" dialog should have value "forge-test-iris2-{time}"
    And Model input in "Apply predictive model" dialog should have value "forge-test-cat-iris-{time}"
    When user clicks on CANCEL button in "Apply predictive model" dialog
    Then the "Apply predictive model" dialog should close
    When user selects "" in "Applicable to" input in forge app
    Then "Applicable to" input in forge app should have value ""
    Given user switches to the "forge-test-iris2-{time}" table view
    And user switches to the "Forge" view
    When user right-clicks on the Forge model "forge-test-cat-iris-{time}" in second grid viewer
    And user picks "Apply..." from the open menu
    Then there should be 1 visible "Apply predictive model" dialog
    And Table input in "Apply predictive model" dialog should have value "forge-test-iris-{time}"
    And Model input in "Apply predictive model" dialog should have value "forge-test-cat-iris-{time}"
    When user clicks on CANCEL button in "Apply predictive model" dialog
    Then the "Apply predictive model" dialog should close
    Given user switches to the "forge-test-iris2-{time}" table view
    When user closes the current view
    And user switches to the "Forge" view
    Then no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Tags are edited in Details and in Edit model..., which renames the model too; the name is never empty (4.7)
    Given the table "forge-test-iris2-{time}" is closed
    And user switches to the "Forge" view
    When user reopens the Forge model "forge-test-cat-iris-{time}" in second grid viewer
    Then the context panel should show the Forge model "forge-test-cat-iris-{time}"
    And there should be 1 visible tag chips in "Details" pane in context panel
    When user remembers the catalog loads
    And user types "second" into tags box in "Details" pane in context panel
    And user presses Enter in tags box in "Details" pane in context panel
    Then there should be 2 visible tag chips in "Details" pane in context panel
    And every catalog should have reloaded
    And the "tags" cell of the Forge model "forge-test-cat-iris-{time}" in second grid viewer should read "demo, second"
    And the Forge model "forge-test-cat-iris-{time}" on the server should have:
      | tags | demo, second |
    When user right-clicks on the Forge model "forge-test-cat-iris-{time}" in second grid viewer
    And user picks "Edit model..." from the open menu
    Then there should be 1 visible "Edit model" dialog
    And Name input in "Edit model" dialog should have value "forge-test-cat-iris-{time}"
    And there should be 2 visible tag chips in "Edit model" dialog
    When user clears Name input in "Edit model" dialog
    Then Name input in "Edit model" dialog should be invalid
    And OK button in "Edit model" dialog should be disabled
    When user enters "forge-test-cat-edited-{time}" into Name input in "Edit model" dialog
    And user types "third" into tags box in "Edit model" dialog
    And user presses Enter in tags box in "Edit model" dialog
    Then there should be 3 visible tag chips in "Edit model" dialog
    When user remembers the catalog loads
    And user clicks on OK button in "Edit model" dialog
    Then the "Edit model" dialog should close
    And an info balloon containing 'Model "forge-test-cat-edited-{time}" updated.' should have been shown
    And every catalog should have reloaded
    And second grid viewer should list the Forge model "forge-test-cat-edited-{time}"
    And second grid viewer should not list the Forge model "forge-test-cat-iris-{time}"
    And the context panel should show the Forge model "forge-test-cat-edited-{time}"
    And the Forge model "forge-test-cat-edited-{time}" on the server should have:
      | tags | demo, second, third |
    And 0 Forge models named "forge-test-cat-iris-{time}" should be on the server
    When user right-clicks on the Forge model "forge-test-cat-edited-{time}" in second grid viewer
    And user picks "Edit model..." from the open menu
    Then there should be 1 visible "Edit model" dialog
    When user enters "forge-test-cat-iris-{time}" into Name input in "Edit model" dialog
    And user remembers the catalog loads
    And user clicks on OK button in "Edit model" dialog
    Then the "Edit model" dialog should close
    And every catalog should have reloaded
    And second grid viewer should list the Forge model "forge-test-cat-iris-{time}"
    And 1 Forge model named "forge-test-cat-iris-{time}" should be on the server
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Download takes the model file out; a row without one offers no Download and is deleted with the trash icon (4.8, 3.7)
    Given user watches downloads
    When user right-clicks on the Forge model "forge-test-cat-iris-{time}" in second grid viewer
    And user picks "Download" from the open menu
    Then a file "forge-test-cat-iris-{time}.bin" should have been downloaded
    And the downloaded file "forge-test-cat-iris-{time}.bin" should have the size of the model file of the Forge model "forge-test-cat-iris-{time}"
    When user right-clicks on the Forge model "forge-test-cat-hand-{time}" in second grid viewer
    Then the open menu should list "Delete model"
    And the open menu should not list "Download"
    When user closes the context menu
    And user clicks on the Forge model "forge-test-cat-hand-{time}" in second grid viewer
    Then the context panel should show the Forge model "forge-test-cat-hand-{time}"
    When user remembers the catalog loads
    And user clicks on "Delete model" icon in forge app
    Then there should be 1 visible "Delete model" dialog
    And "Delete model" dialog should contain text 'Delete the model "forge-test-cat-hand-{time}"?'
    When user clicks on OK button in "Delete model" dialog
    Then every catalog should have reloaded
    And second grid viewer should not list the Forge model "forge-test-cat-hand-{time}"
    And 0 Forge models named "forge-test-cat-hand-{time}" should be on the server
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Applicable to keeps the models an open table fits, in the browser; a closed table leaves it (4.11)
    When user selects "forge-test-iris-{time}" in "Applicable to" input in forge app
    Then second grid viewer should list the Forge model "forge-test-cat-iris-{time}"
    And second grid viewer should list the Forge model "forge-test-cat-ref-{time}"
    And second grid viewer should not list the Forge model "forge-test-cat-quality-{time}"
    When user selects "forge-test-wine-{time}" in "Applicable to" input in forge app
    Then second grid viewer should list the Forge model "forge-test-cat-quality-{time}"
    And second grid viewer should not list the Forge model "forge-test-cat-iris-{time}"
    And second grid viewer should not list the Forge model "forge-test-cat-ref-{time}"
    When user selects "" in "Applicable to" input in forge app
    Then second grid viewer should list the Forge model "forge-test-cat-iris-{time}"
    And second grid viewer should list the Forge model "forge-test-cat-ref-{time}"
    And second grid viewer should list the Forge model "forge-test-cat-quality-{time}"
    When user selects "forge-test-wine-{time}" in "Applicable to" input in forge app
    And user switches to the "forge-test-wine-{time}" table view
    And user closes the current view
    And user switches to the "Forge" view
    Then "Applicable to" input in forge app should have value ""
    And second grid viewer should list the Forge model "forge-test-cat-iris-{time}"
    And "Applicable to" input in forge app should contain text "forge-test-iris-{time}"
    And "Applicable to" input in forge app should not contain text "forge-test-wine-{time}"
    And no table named "forge-test-iris-{time}" should be on the server
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Shared with a user through the platform's dialog, the model shows the grant; the grant is taken back (4.12)
    When user reopens the Forge model "forge-test-cat-quality-{time}" in second grid viewer
    Then the context panel should show the Forge model "forge-test-cat-quality-{time}"
    When user expands "Sharing" pane in context panel
    Then "Sharing" pane in context panel should contain text "Not shared yet. Only its author and administrators can see this model."
    When user clicks on "Share..." button in "Sharing" pane in context panel
    Then there should be 1 visible "Share forge-test-cat-quality-{time}" dialog
    And "Share forge-test-cat-quality-{time}" dialog should contain text "Full access"
    When user picks the sharing user by name in "User, group, or email" input in "Share forge-test-cat-quality-{time}" dialog
    And user unchecks "Send notifications" input in "Share forge-test-cat-quality-{time}" dialog
    And user clicks on OK button in "Share forge-test-cat-quality-{time}" dialog
    Then the "Share forge-test-cat-quality-{time}" dialog should close
    And an info balloon containing "Shared" should have been shown
    And value of "Can view" detail in context panel should name the sharing user
    And the Forge model "forge-test-cat-quality-{time}" should be shared with the sharing user on the server
    When user clicks on "Share..." button in "Sharing" pane in context panel
    Then there should be 1 visible "Share forge-test-cat-quality-{time}" dialog
    And "Share forge-test-cat-quality-{time}" dialog should contain text "Full access"
    When user removes the sharing user by name from "Share forge-test-cat-quality-{time}" dialog
    And user clicks on OK button in "Share forge-test-cat-quality-{time}" dialog
    Then the "Share forge-test-cat-quality-{time}" dialog should close
    And "Sharing" pane in context panel should contain text "Not shared yet. Only its author and administrators can see this model."
    And the Forge model "forge-test-cat-quality-{time}" should not be shared with the sharing user on the server
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Run Evaluation measures a model on its data or an open table, without a record or a column (5.15)
    Given user counts the task bar entries
    When user clicks on the Forge model "forge-test-cat-ref-{time}" in second grid viewer
    Then the context panel should show the Forge model "forge-test-cat-ref-{time}"
    When user expands "Performance" pane in context panel
    And user clicks on "Run Evaluation" button in "Performance" pane in context panel
    Then the "Evaluation: forge-test-cat-ref-{time}" view should be current
    And results area should become ready
    And the task bar entry that showed "Evaluating forge-test-cat-ref-{time}" should have ended
    And the task bar should have shown 1 entry "Evaluating forge-test-cat-ref-{time}"
    And "Performance" results panel should be visible
    And the "text of cell 1 of Evaluation" reading of grid viewer in "Performance" results panel should equal the "text of cell 1 of Train" reading of grid viewer in "Performance" results panel
    And the "text of cell 2 of Evaluation" reading of grid viewer in "Performance" results panel should equal the "text of cell 2 of Train" reading of grid viewer in "Performance" results panel
    And "Performance" results panel should contain text "Evaluation: iris (150 rows)"
    And "Predicted Species vs Actual" results panel should be visible
    And the open table views should be exactly "forge-test-iris-{time}"
    And table "forge-test-iris-{time}" should have columns "col 1, Sepal.Length, Sepal.Width, Petal.Length, Petal.Width, Species, Species (predicted)"
    And the model "forge-test-cat-ref-{time}" should have 0 applications on the server
    When user closes the current view
    And user switches to the "Forge" view
    And user clicks on the Forge model "forge-test-cat-iris-{time}" in second grid viewer
    Then the context panel should show the Forge model "forge-test-cat-iris-{time}"
    When user expands "Performance" pane in context panel
    And user clicks on "Run Evaluation" button in "Performance" pane in context panel
    Then there should be 1 visible "Evaluate model" dialog
    And Table input in "Evaluate model" dialog should have value "forge-test-iris-{time}"
    And editor of Target input in "Evaluate model" dialog should have text "Species"
    When user clicks on OK button in "Evaluate model" dialog
    Then the "Evaluation: forge-test-cat-iris-{time}" view should be current
    And results area should become ready
    And the task bar entry that showed "Evaluating forge-test-cat-iris-{time}" should have ended
    And the task bar should have shown 1 entry "Evaluating forge-test-cat-iris-{time}"
    And "Performance" results panel should contain text "Evaluation: forge-test-iris-{time} (150 rows)"
    And table "forge-test-iris-{time}" should have columns "col 1, Sepal.Length, Sepal.Width, Petal.Length, Petal.Width, Species, Species (predicted)"
    And the table "forge-test-iris-{time}" should have no mouse-over row
    When user hovers over the "marker of row 5" area of scatter plot viewer in "Predicted Species vs Actual" results panel
    Then the mouse-over row of the table "forge-test-iris-{time}" should be the "hovered row" reading of scatter plot viewer in "Predicted Species vs Actual" results panel
    And the model "forge-test-cat-iris-{time}" should have 1 application on the server
    When user closes the current view
    And user switches to the "forge-test-iris-{time}" table view
    And user closes the current view
    And user switches to the "Forge" view
    And user reopens the Forge model "forge-test-cat-iris-{time}" in second grid viewer
    Then the context panel should show the Forge model "forge-test-cat-iris-{time}"
    When user expands "Performance" pane in context panel
    And user clicks on "Run Evaluation" button in "Performance" pane in context panel
    Then a warning balloon containing "No open table has the columns the model needs. Open its training data or a similar table." should have been shown
    When user clicks on the Forge model "forge-test-cat-ref-{time}" in second grid viewer
    Then the context panel should show the Forge model "forge-test-cat-ref-{time}"
    When user expands "Performance" pane in context panel
    And user clicks on "Run Evaluation" button in "Performance" pane in context panel
    Then the "Evaluation: forge-test-cat-ref-{time}" view should be current
    And results area should become ready
    And no table should be left in the workspace
    When user closes the current view
    Then no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Deleted with the trash icon in a second catalog, the model leaves both catalogs, its file and the panel; its run stays (2.6)
    Given user opens the "Forge" app
    Then "Delete model" icon in forge app should be disabled
    And the newest training run of the table "forge-test-wine-{time}" should be linked to the model "forge-test-cat-quality-{time}"
    When user clicks on the Forge model "forge-test-cat-quality-{time}" in second grid viewer
    Then "Delete model" icon in forge app should be enabled
    And the context panel should show the Forge model "forge-test-cat-quality-{time}"
    When user hovers over "Delete model" icon in forge app
    Then tooltip should contain text "Delete model"
    When user remembers the catalog loads
    And user clicks on "Delete model" icon in forge app
    Then there should be 1 visible "Delete model" dialog
    And "Delete model" dialog should contain text 'Delete the model "forge-test-cat-quality-{time}"?'
    When user clicks on OK button in "Delete model" dialog
    Then there should be 2 open views named "Forge"
    And every catalog should have reloaded
    And second grid viewer should not list the Forge model "forge-test-cat-quality-{time}"
    And the context panel should not show the Forge model "forge-test-cat-quality-{time}"
    Given user switches to the "Forge" view
    Then the first Forge view should be shown
    And the shown catalog should have reloaded
    And second grid viewer should not list the Forge model "forge-test-cat-quality-{time}"
    And second grid viewer should list the Forge model "forge-test-cat-iris-{time}"
    When user remembers the catalog loads
    And user clicks on "Refresh" icon in forge app
    Then the shown catalog should have reloaded
    And second grid viewer should not list the Forge model "forge-test-cat-quality-{time}"
    And second grid viewer should list the Forge model "forge-test-cat-iris-{time}"
    And the model file folder of the deleted Forge model "forge-test-cat-quality-{time}" should be gone
    And the newest training run of the table "forge-test-wine-{time}" should have:
      | status   | completed |
      | model_id | empty     |
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Deleted from the row menu, the model goes with its file and its applications (4.9)
    Then the model "forge-test-cat-iris-{time}" should have 1 application on the server
    When user reopens the Forge model "forge-test-cat-iris-{time}" in second grid viewer
    Then the context panel should show the Forge model "forge-test-cat-iris-{time}"
    When user right-clicks on the Forge model "forge-test-cat-iris-{time}" in second grid viewer
    And user remembers the catalog loads
    And user picks "Delete model" from the open menu
    Then there should be 1 visible "Delete model" dialog
    And "Delete model" dialog should contain text 'Delete the model "forge-test-cat-iris-{time}"?'
    When user clicks on OK button in "Delete model" dialog
    Then every catalog should have reloaded
    And second grid viewer should not list the Forge model "forge-test-cat-iris-{time}"
    And the context panel should not show the Forge model "forge-test-cat-iris-{time}"
    And 0 Forge models named "forge-test-cat-iris-{time}" should be on the server
    And the model file folder of the deleted Forge model "forge-test-cat-iris-{time}" should be gone
    And the model "forge-test-cat-iris-{time}" should have 0 applications on the server
    And no error or warning balloon should have been shown
    And no errors should have been logged
