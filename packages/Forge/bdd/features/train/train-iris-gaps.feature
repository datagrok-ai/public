@journey @serial
Feature: Training on iris with missing values: Skip rows, Impute, and a cell cleared later
  ML | Forge | Train... over iris opened from the demo files, with two cells of Sepal.Width (rows 3 and 10) and one of
  Species (row 20) cleared before the view opens. Translated from the manual scenarios 3.11, 3.12, 5.11 (steps 1-3)
  and 5.13 (step 3) of core/docs/features/platform/predictive-modeling/forge/manual-checks.md. A training is claimed
  by the run its results show (data-forge-run-id): a new run's result, then that run as the table's newest on the
  server with the rows it used and the preparation it recorded (options.missingValues as a whole,
  options.preprocessingInfo), and the seed the Performance panel shows as that run's seed.
  The feature saves no model and its runs are under its own table name; it is serial all the same, to keep a run to
  two lanes (three sessions starting at once have stalled the stand). The hovered point of 5.13 is claimed by its
  values: the plot draws the 146 rows the training kept, so with rows 3, 5, 10 and 20 skipped its row 17 is the
  table's row 21, and the table's mouse-over row holds the point's four measurements (the manual's "the skipped row
  shifts nothing").

  Not claimed here: the two options of Missing values on one line and "Method follows Data directly" with no gap
  (visual: the order of the shown groups is claimed); the iris grid's two cells "still empty" (read from the table's
  column, not from the grid's drawing); the grid highlight of a hovered point (the table's mouse-over row is read
  through the API); 5.11 step 4 (the deletes of the demog feature's models: its sweep). Neighbors after a Skip rows /
  Impute round trip (a hidden setting goes back to its default, so 3 is claimed only before step 5).

  Where the product moved on from the manual text: the cells are cleared before the view opens (3.11's own order) and
  in 5.13 step 3 while it is open: the Train view does not react to a cell edit (no data-change subscription, as in the
  old tool, which reacts to added and removed columns only), so a cleared cell counts at the next input change. The
  scenario clears a Petal.Length cell while no checked feature has a gap (Preparation hidden), and the next change,
  Petal.Width unchecked, brings Preparation back with the gap counted; with every feature checked again the run uses
  146 rows. The feature trains with SVM: XGBoost, the default method, handles missing values itself (`mlnativenans`)
  and hides Missing values, so its kept gaps are claimed once, in the first scenario.

  Background:
    Given user is logged in
    And no Forge training runs of the table "forge-test-irisgaps-{time}" are on the server
    And user opens the demo file "iris.csv" as "forge-test-irisgaps-{time}"
    When user clears "Sepal.Width" column in rows "3, 10"
    And user clears "Species" column in rows "20"
    And user remembers the shown training
    And user picks "ML > Forge > Train..." from the top menu
    Then the "Predictive model" view should be current
    And a new training's result should be shown
    When user selects "SVM" in Method input
    Then a new training's result should be shown
    And no balloon of any colour should have been shown
    And no errors should have been logged

  Scenario: Rows with missing values are skipped and counted; Missing values lives in Preparation (3.11, 5.11)
    Then the groups of train form should be "Data, Preparation, Method"
    And the "Preparation" group should show exactly the inputs "Missing values"
    And "Missing values" input in "Data" group should be absent
    And "Missing values" input in "Preparation" group should have "Skip rows" chosen
    And "Missing values" input should contain text "Impute"
    And editor of "Missing values" input should be visible
    And editor of "Missing values" input should be valid
    And Train button should be hidden
    When user hovers over label of "Missing values" input
    Then tooltip should contain text "Rows with missing values in: Sepal.Width: 2."
    When user hovers over label of Target input
    Then tooltip should contain text "Species: 1 row without a value is skipped."
    And no balloon of any colour should have been shown
    And "Performance" results panel should contain text "Rows: 147 used, 3 skipped (missing values)"
    And "Performance" results panel should contain text "Validation: 5-fold cross-validation on 147 rows"
    And the shown seed should be the seed of the shown training on the server
    And the "text of cell 1 of Metric" reading of grid viewer in "Performance" results panel should be "Accuracy"
    And "Predicted Species vs Actual" results panel should be visible
    And "Confusions" results panel should be visible
    And "Class balance" results panel should be visible
    And the shown training should be the newest training run of the table "forge-test-irisgaps-{time}", with:
      | status                    | completed                        |
      | row_count                 | 147                              |
      | options.preprocessingInfo | ["ignore-missing"]               |
      | options.missingValues     | {"mode":"skip","skippedRows":3} |
    And the "text of cell 1 of Preparation" reading of grid viewer in "Runs" results panel should be ""
    When user selects "XGBoost" in Method input
    Then a new training's result should be shown
    And the groups of train form should be "Data, Method"
    And "Performance" results panel should contain text "Missing values: kept (the method handles them)"
    And "Performance" results panel should contain text "Validation: 5-fold cross-validation on 149 rows"
    And the shown training should be the newest training run of the table "forge-test-irisgaps-{time}", with:
      | status                    | completed                        |
      | row_count                 | 149                              |
      | options.preprocessingInfo | []                               |
      | options.missingValues     | {"mode":"keep","skippedRows":1} |
    When user selects "SVM" in Method input
    Then a new training's result should be shown
    And the groups of train form should be "Data, Preparation, Method"
    When user clicks on editor of Features input
    And user toggles the "Sepal.Width" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then a new training's result should be shown
    And the groups of train form should be "Data, Method"
    And "Performance" results panel should contain text "Rows: 149 used, 1 skipped (missing values)"
    When user clicks on editor of Features input
    And user toggles the "Sepal.Width" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then a new training's result should be shown
    And the groups of train form should be "Data, Preparation, Method"
    And "Missing values" input in "Preparation" group should have "Skip rows" chosen
    And "Performance" results panel should contain text "Rows: 147 used, 3 skipped (missing values)"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Impute fills the features' gaps inside Preparation; an invalid setting trains nothing (3.12, 5.11)
    When user selects "Impute" in "Missing values" input
    Then there should be 0 visible dialog
    And Neighbors input should have value "4"
    And Distance input should have value "Euclidean"
    And Distance input should offer "Euclidean, Manhattan"
    And the "Preparation" group should show exactly the inputs "Missing values, Neighbors, Distance"
    And a new training's result should be shown
    And "Performance" results panel should contain text "Rows: 149 used, 1 skipped (missing values)"
    And the shown training should be the newest training run of the table "forge-test-irisgaps-{time}", with:
      | row_count             | 149                                                                       |
      | options.missingValues | {"mode":"impute","distance":"Euclidean","neighbors":4,"skippedRows":1} |
    When user hovers over label of Neighbors input
    Then tooltip should contain text "Nearest rows used to fill a value."
    When user remembers the training run count of the table "forge-test-irisgaps-{time}"
    And user enters "0" into Neighbors input
    Then editor of Neighbors input should be invalid
    And results area should be invalid
    And "Performance" results panel should contain text "Fix the settings."
    And Train button should be hidden
    When user enters "3" into Neighbors input
    Then a new training's result should be shown
    And the training run count of the table "forge-test-irisgaps-{time}" should be one more than remembered
    And "Performance" results panel should contain text "Rows: 149 used, 1 skipped (missing values)"
    And the "text of cell 1 of Metric" reading of grid viewer in "Performance" results panel should be "Accuracy"
    And "Predicted Species vs Actual" results panel should be visible
    And "Confusions" results panel should be visible
    And "Class balance" results panel should be visible
    And the shown training should be the newest training run of the table "forge-test-irisgaps-{time}", with:
      | status                    | completed                                                                 |
      | row_count                 | 149                                                                       |
      | options.preprocessingInfo | ["impute-missing"]                                                        |
      | options.missingValues     | {"mode":"impute","distance":"Euclidean","neighbors":3,"skippedRows":1} |
    And the "text of cell 1 of Preparation" reading of grid viewer in "Runs" results panel should be "impute (3 neighbors, Euclidean)"
    When user switches to the "forge-test-irisgaps-{time}" table view
    Then "Sepal.Width" column should have 2 missing values
    And the value of "Sepal.Width" column in row 3 should be ""
    And the value of "Sepal.Width" column in row 10 should be ""
    When user switches to the "Predictive model" view
    And user selects "Skip rows" in "Missing values" input
    Then the "Preparation" group should show exactly the inputs "Missing values"
    And a new training's result should be shown
    And "Performance" results panel should contain text "Rows: 147 used, 3 skipped (missing values)"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: A cell cleared while the view is open counts at the next change; a hovered point is still the table's row (5.13)
    When user clicks on editor of Features input
    And user toggles the "Sepal.Width" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then a new training's result should be shown
    And the groups of train form should be "Data, Method"
    When user switches to the "forge-test-irisgaps-{time}" table view
    And user clears "Petal.Length" column in rows "5"
    And user switches to the "Predictive model" view
    And user clicks on editor of Features input
    And user toggles the "Petal.Width" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then a new training's result should be shown
    And editor of Features input should have text "(2) Sepal.Length, Petal.Length"
    And the groups of train form should be "Data, Preparation, Method"
    And "Performance" results panel should contain text "Rows: 148 used, 2 skipped (missing values)"
    When user hovers over label of "Missing values" input
    Then tooltip should contain text "Rows with missing values in: Petal.Length: 1."
    When user clicks on editor of Features input
    And user toggles the "Sepal.Width" column in the column list of "Select columns..." dialog
    And user toggles the "Petal.Width" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then a new training's result should be shown
    And editor of Features input should have text "(4) Sepal.Length, Petal.Length, Sepal.Width, Petal.Width"
    And "Performance" results panel should contain text "Rows: 146 used, 4 skipped (missing values)"
    And the shown training should be the newest training run of the table "forge-test-irisgaps-{time}", with:
      | row_count             | 146                              |
      | options.missingValues | {"mode":"skip","skippedRows":4} |
    When user hovers over label of "Missing values" input
    Then tooltip should contain text "Petal.Length: 1"
    And tooltip should contain text "Sepal.Width: 2"
    And the "rows shown" reading of scatter plot viewer in "Predicted Species vs Actual" results panel should be 146
    And the table "forge-test-irisgaps-{time}" should have no mouse-over row
    When user hovers over the "marker of row 17" area of scatter plot viewer in "Predicted Species vs Actual" results panel
    Then the mouse-over row of the table "forge-test-irisgaps-{time}" should hold the "Sepal.Length, Sepal.Width, Petal.Length, Petal.Width" values of the hovered row of scatter plot viewer in "Predicted Species vs Actual" results panel
    And no error or warning balloon should have been shown
    And no errors should have been logged
