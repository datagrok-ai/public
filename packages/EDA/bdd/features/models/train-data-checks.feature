@journey @eda @realizes:ml.menu.models.train-model
Feature: The Train Model view checks the data before it trains
  ML | Models | Train Model... runs its data checks as soon as the target and the features are set and
  lists what they find under Insights & Tips in the preview: class imbalance, categorical features,
  identifier-like categories, highly correlated columns and missing values. Translated from
  files/TestTrack/Models/models-validators-edge.md and models-bug-grok-3525.md (GROK-3525) and their
  playwright-public specs, which built 30–200-row tables through the JS API and checked the tip text
  and a Save button that was only attached.

  Each scenario opens a table of its own, written into the feature, small enough that the numbers the
  checks print can be stated: eight A and two B make A 1.600 and B 0.400 of the two classes' even share.
  A check that only warns leaves the training going (the preview becomes ready and Save enabled). A
  missing value is an error that stops the training; a categorical or identifier-like feature leaves no
  engine that takes the data. Either way the preview reports an invalid result, Save stays disabled,
  and the view offers its remedy (Ignore missing, One-hot encoding, Skip unique categories); it trains
  once One-hot encoding or Ignore missing is ticked. Each scenario closes its Train Model view.

  Not translated, and why: the old cases' claim that a categorical feature "does not block training"
  — no engine of the stand takes a text feature, so the view shows no Model Engine until One-hot
  encoding is on (the scenario says so); the "dapi.ml.save was not called" checks of GROK-3525, since
  nothing reaches Save while it is disabled; the categorical-target variant on demog, whose RACE has no
  missing values in demog-1000 (a table with an empty label stands in). The tip for an identifier-like
  column reads "contains contain too many unique categories" (MISSING.md, to check by hand), so the
  claim names only its last words. Ticking Skip unique categories offers the engines but trains
  nothing — the preview stays empty and "Invalid argument (namesOrColumns): Not supported type: null"
  is logged from PredictiveModelingEngine.apply — so the scenario stops before the tick; the rest waits
  in MISSING.md for its ticket.

  Background:
    Given user is logged in

  Scenario: Class imbalance is reported and the model still trains
    Given user opens a table "imbalanced" with:
      | f1 | f2 | target |
      | 1  | 7  | A      |
      | 2  | 3  | A      |
      | 3  | 9  | A      |
      | 4  | 1  | A      |
      | 5  | 8  | A      |
      | 6  | 2  | A      |
      | 7  | 10 | A      |
      | 8  | 4  | A      |
      | 9  | 6  | B      |
      | 10 | 5  | B      |
    When user picks "ML > Models > Train Model..." from the top menu
    Then the "Predictive model" view should be current
    When user selects "target" in Predict input
    And user clicks on editor of Features input
    And user clicks on None label in "Select columns..." dialog
    And user clicks on the "cell 1 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on the "cell 2 of x" area of grid viewer in "Select columns..." dialog
    Then "Select columns..." dialog should contain text "2 checked"
    When user clicks on OK button in "Select columns..." dialog
    Then editor of Features input should contain text "(2)"
    And model preview should be ready
    And model preview should contain text "Some columns contain class imbalance"
    And model preview should contain text "A (1.600)"
    And model preview should contain text "B (0.400)"
    And "Model Engine" input should be visible
    And Save button should be enabled
    And no error or warning balloon should have been shown
    And no errors should have been logged
    When user closes the current view

  Scenario: A categorical feature is reported, and One-hot encoding lets the model train
    Given user opens a table "colors" with:
      | cat   | num | target |
      | red   | 1   | 3.1    |
      | green | 2   | 1.2    |
      | blue  | 3   | 4.4    |
      | red   | 4   | 2.5    |
      | green | 5   | 5.3    |
      | blue  | 6   | 1.7    |
      | red   | 7   | 3.9    |
      | green | 8   | 2.2    |
      | blue  | 9   | 4.8    |
    When user picks "ML > Models > Train Model..." from the top menu
    And user selects "target" in Predict input
    And user clicks on editor of Features input
    And user clicks on None label in "Select columns..." dialog
    And user clicks on the "cell 1 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on the "cell 2 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then editor of Features input should contain text "(2)"
    And model preview should contain text "Column 'cat' is categorical. Most models require converting them to numerical."
    And model preview should be invalid
    And "One-hot encoding" input should be visible
    And "Model Engine" input should be absent
    And Save button should be disabled
    When user checks "One-hot encoding" input
    Then model preview should be ready
    And "Model Engine" input should be visible
    And model preview should not contain text "is categorical"
    And "R squared" table row should be visible
    And Save button should be enabled
    And no error or warning balloon should have been shown
    And no errors should have been logged
    When user closes the current view

  Scenario: An identifier-like column is reported, and the view offers to skip it
    Given user opens a table "identifiers" with:
      | id_like | feature1 | target |
      | row_1   | 1        | 3.1    |
      | row_2   | 2        | 1.2    |
      | row_3   | 3        | 4.4    |
      | row_4   | 4        | 2.5    |
      | row_5   | 5        | 5.3    |
      | row_6   | 6        | 1.7    |
    When user picks "ML > Models > Train Model..." from the top menu
    And user selects "target" in Predict input
    And user clicks on editor of Features input
    And user clicks on None label in "Select columns..." dialog
    And user clicks on the "cell 1 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on the "cell 2 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then editor of Features input should contain text "(2)"
    And model preview should contain text "Column 'id_like' contains"
    And model preview should contain text "too many unique categories."
    And model preview should be invalid
    And "Skip unique categories" input should be visible
    And "Model Engine" input should be absent
    And Save button should be disabled
    And no error or warning balloon should have been shown
    And no errors should have been logged
    When user closes the current view

  Scenario: Highly correlated features are named in pairs and the model still trains
    Given user opens a table "correlated" with:
      | feat_a | feat_b | target |
      | 1      | 1.1    | 3.1    |
      | 2      | 2.3    | 1.2    |
      | 3      | 2.9    | 4.4    |
      | 4      | 4.2    | 2.5    |
      | 5      | 5.1    | 5.3    |
      | 6      | 6.4    | 1.7    |
      | 7      | 6.8    | 3.9    |
      | 8      | 8.3    | 2.2    |
    When user picks "ML > Models > Train Model..." from the top menu
    And user selects "target" in Predict input
    And user clicks on editor of Features input
    And user clicks on None label in "Select columns..." dialog
    And user clicks on the "cell 1 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on the "cell 2 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then editor of Features input should contain text "(2)"
    And model preview should be ready
    And model preview should contain text "Columns are highly correlated"
    And model preview should contain text "feat_a <> feat_b"
    And Save button should be enabled
    And no error or warning balloon should have been shown
    And no errors should have been logged
    When user closes the current view

  # GROK-3525: a target with missing values is reported, and nothing trains until the rows are dropped
  Scenario: Missing values in a numeric target block training until Ignore missing is ticked
    Given user opens a table "gaps" with:
      | X  | Y  |
      | 1  | 2  |
      | 2  | 4  |
      | 3  |    |
      | 4  | 8  |
      | 5  | 10 |
      | 6  |    |
      | 7  | 14 |
      | 8  | 16 |
      | 9  | 18 |
      | 10 | 20 |
    When user picks "ML > Models > Train Model..." from the top menu
    And user selects "Y" in Predict input
    And user clicks on editor of Features input
    And user clicks on None label in "Select columns..." dialog
    And user clicks on the "cell 1 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then editor of Features input should contain text "(1)"
    And model preview should contain text "Column 'Y' contains missing values."
    And model preview should be invalid
    And "Ignore missing" input should be visible
    And "Impute missing" input should be visible
    And "Model Engine" input should be absent
    And Save button should be disabled
    When user checks "Ignore missing" input
    Then model preview should be ready
    And "Impute missing" input should be hidden
    And "Model Engine" input should be visible
    And model preview should not contain text "contains missing values"
    And "R squared" table row should be visible
    And Save button should be enabled
    And no error or warning balloon should have been shown
    And no errors should have been logged
    When user closes the current view

  # GROK-3525, the categorical branch
  Scenario: Missing labels in a categorical target block training too
    Given user opens a table "unlabelled" with:
      | f1 | f2 | label |
      | 1  | 7  | A     |
      | 2  | 3  | B     |
      | 3  | 9  |       |
      | 4  | 1  | A     |
      | 5  | 8  | B     |
      | 6  | 2  | A     |
      | 7  | 10 | B     |
      | 8  | 4  |       |
    When user picks "ML > Models > Train Model..." from the top menu
    And user selects "label" in Predict input
    And user clicks on editor of Features input
    And user clicks on None label in "Select columns..." dialog
    And user clicks on the "cell 1 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on the "cell 2 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then editor of Features input should contain text "(2)"
    And model preview should contain text "Column 'label' contains missing values."
    And model preview should be invalid
    And "Model Engine" input should be absent
    And Save button should be disabled
    When user checks "Ignore missing" input
    Then model preview should be ready
    And "Model Engine" input should be visible
    And model preview should not contain text "contains missing values"
    And Save button should be enabled
    And no error or warning balloon should have been shown
    And no errors should have been logged
    When user closes the current view
