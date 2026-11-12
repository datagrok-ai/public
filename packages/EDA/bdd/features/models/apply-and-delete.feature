@journey @eda @realizes:ml.menu.models.train-model @realizes:ml.menu.models.apply-model @realizes:eda.model.pls-regression @realizes:eda.model.linear-regression
Feature: Saved models applied to new data and deleted from the gallery
  Two regression models — Petal.Width predicted from Sepal.Length, Sepal.Width and Petal.Length, by PLS
  Regression and by Linear Regression — trained on iris, saved, applied through ML | Models | Apply
  Model... to a fresh copy of the table and to another table that has the input columns, then deleted
  from Browse > Platform > Predictive models. Translated from TestTrack General/predictive-models-spec.ts.

  The Apply dialog offers the models suggested for the table, each shown as "<saved at>: <name>" cut at
  40 characters and keyed by its id, so the model is chosen by name through its id; a change of the
  choice makes the chosen model the current object, which the context panel shows.

  Not translated, and why: General/chemprop-spec.ts — the Chemprop engine trains and predicts in the
  chem-chemprop Docker container (the lead's rule); its in-browser half (the Train Model view, Predict
  and Features) is claimed here and in train-on-cars. The old spec trained on
  sensors/accelerometer.csv (accel_x from accel_y, accel_z and time_offset), which is not a registered
  dataset, and applied a model to a grok.data.testData random walk, which no step opens and whose
  columns are none of the model's inputs; iris and a table with the input columns stand in (MISSING.md).
  The view asks every engine whether it applies, the Samples package's Docker engine (PyKNN) included: a
  stand where Samples is installed without its container logs "Container is not started". The view remembers hyperparameters, so the
  Components claim reads the value the view offers, which is the default only on a fresh account.

  PLS is trained with two components: with three on three features it spans them all and predicts as the
  linear regression does, so the two applied columns could not be told apart. Which model OK applied is
  claimed by the two prediction columns differing; the model a column names in its predictive.model tag
  is not readable by any step. How well a model predicts is the EDA package tests' claim, not this
  feature's. Both models are removed through the gallery in the last scenario, and by name now and at
  feature end.

  Background:
    Given user is logged in
    And no predictive model named "BDD-Iris-PLS-{run}" is on the server
    And no predictive model named "BDD-Iris-LR-{run}" is on the server

  Scenario: A PLS model predicts Petal.Width from three other columns, and is saved
    Given user opens iris dataset
    Then the table should have 6 columns
    When user picks "ML > Models > Train Model..." from the top menu
    Then the "Predictive model" view should be current
    When user selects "Petal.Width" in Predict input
    And user clicks on editor of Features input
    And user clicks on None label in "Select columns..." dialog
    Then the "text of cell 2 of __name" reading of grid viewer in "Select columns..." dialog should be "Sepal.Length"
    And the "text of cell 4 of __name" reading of grid viewer in "Select columns..." dialog should be "Petal.Length"
    When user clicks on the "cell 2 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on the "cell 3 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on the "cell 4 of x" area of grid viewer in "Select columns..." dialog
    Then "Select columns..." dialog should contain text "3 checked"
    When user clicks on OK button in "Select columns..." dialog
    Then editor of Features input should contain text "(3)"
    When user selects "Eda: PLS Regression" in "Model Engine" input
    Then Components input should have value "3"
    And model preview should be ready
    # three components on three features would span them all and predict as the linear regression does
    When user enters "2" into Components input
    Then Components input should have value "2"
    And model preview should be ready
    And "Eda: PLS Regression" heading should be visible
    And "R squared" table row should be visible
    When user clicks on Save button
    And user enters "BDD-Iris-PLS-{run}" into Name input in dialog
    And user clicks on OK button in dialog
    Then dialog should be absent
    And 1 predictive model named "BDD-Iris-PLS-{run}" should be on the server
    And no errors should have been logged

  Scenario: A linear regression on the same inputs is saved too
    When user selects "Eda: Linear Regression" in "Model Engine" input
    Then model preview should be ready
    And "Eda: Linear Regression" heading should be visible
    And "R squared" table row should be visible
    When user clicks on Save button
    And user enters "BDD-Iris-LR-{run}" into Name input in dialog
    And user clicks on OK button in dialog
    Then dialog should be absent
    And 1 predictive model named "BDD-Iris-LR-{run}" should be on the server
    And no errors should have been logged

  Scenario: Apply Model offers the saved models and adds the chosen one's prediction
    Given user opens iris dataset
    And the context panel is open
    When user picks "ML > Models > Apply Model..." from the top menu
    Then "Apply predictive model" dialog should be visible
    When user selects the predictive model "BDD-Iris-PLS-{run}" in Model input in "Apply predictive model" dialog
    Then the context panel should show "BDD-Iris-PLS-{run}"
    When user selects the predictive model "BDD-Iris-LR-{run}" in Model input in "Apply predictive model" dialog
    Then the context panel should show "BDD-Iris-LR-{run}"
    And Inputs input in "Apply predictive model" dialog should contain text "(3/3)"
    When user clicks on OK button in "Apply predictive model" dialog
    Then the "Apply predictive model" dialog should close
    And 1 new column should have been added
    And the table should have 7 columns
    And the newest column matching "^Petal\.Width" should have no missing values
    And no errors should have been logged

  Scenario: The other model adds its own prediction beside it
    When user picks "ML > Models > Apply Model..." from the top menu
    Then "Apply predictive model" dialog should be visible
    When user selects the predictive model "BDD-Iris-PLS-{run}" in Model input in "Apply predictive model" dialog
    Then the context panel should show "BDD-Iris-PLS-{run}"
    And Inputs input in "Apply predictive model" dialog should contain text "(3/3)"
    When user clicks on OK button in "Apply predictive model" dialog
    Then the "Apply predictive model" dialog should close
    And 1 new column should have been added
    And the table should have 8 columns
    And the newest column matching "^Petal\.Width" should have no missing values
    # the linear regression wrote "Petal.Width (2)" one scenario earlier; with two components the PLS predicts otherwise
    And some value of "Petal.Width (3)" column should differ from "Petal.Width (2)" column in the same row
    And no errors should have been logged

  Scenario: A model applies to another table that has its input columns
    Given user opens a table "new readings" with:
      | Sepal.Length | Sepal.Width | Petal.Length |
      | 5.1          | 3.5         | 1.4          |
      | 6.4          | 3.2         | 4.5          |
      | 6.3          | 3.3         | 6.0          |
    When user picks "ML > Models > Apply Model..." from the top menu
    Then "Apply predictive model" dialog should be visible
    When user selects the predictive model "BDD-Iris-LR-{run}" in Model input in "Apply predictive model" dialog
    Then the context panel should show "BDD-Iris-LR-{run}"
    And Inputs input in "Apply predictive model" dialog should contain text "(3/3)"
    When user clicks on OK button in "Apply predictive model" dialog
    Then the "Apply predictive model" dialog should close
    And 1 new column should have been added
    And the table should have 4 columns
    And the newest column matching "^Petal\.Width" should have no missing values
    And no errors should have been logged

  Scenario: The gallery shows the model's details and performance, and deletes it
    Given the browse panel is open
    And Platform tree node inside browse tree is expanded
    When user clicks on "Predictive models" tree node inside browse tree
    Then the "Models" view should be current
    When user clicks on "BDD-Iris-LR-{run}" label in gallery
    Then the context panel should show "BDD-Iris-LR-{run}"
    And "Details" pane in context panel should be visible
    And "Performance" pane in context panel should be visible
    And "Sharing" pane in context panel should be visible
    When user picks "Delete" from the context menu of "BDD-Iris-LR-{run}" label in gallery
    Then "Are you sure?" dialog should be visible
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And "BDD-Iris-LR-{run}" label in gallery should be absent
    And 0 predictive models named "BDD-Iris-LR-{run}" should be on the server
    And 1 predictive model named "BDD-Iris-PLS-{run}" should be on the server
    When user picks "Delete" from the context menu of "BDD-Iris-PLS-{run}" label in gallery
    And user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And "BDD-Iris-PLS-{run}" label in gallery should be absent
    And 0 predictive models named "BDD-Iris-PLS-{run}" should be on the server
    And no errors should have been logged
