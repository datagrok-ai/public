@journey
Feature: The Forge app, the platform's entry points that open it
  The Forge app from Browse > Apps, from Browse > Platform > Predictive models and from Forge:forgeModels, its Methods
  grid, its catalog's columns and icons, and ML | Models opening Forge's Train view and Apply dialog. Translated from
  the manual scenarios 1.1, 1.2, 1.3, 1.4, 1.6, 2.9, 4.13 (steps 1-2) and 6.4 of
  core/docs/features/platform/predictive-modeling/forge/manual-checks.md. Only the Train view's training run of the
  feature's own table is saved, swept at the feature's start and end, so the feature runs in the parallel lane; it opens iris
  from the demo files because the top menu shows ML only over a table view. The grids'
  captions are read as the grids draw them (the header cell under each header), and the word "Engine" is looked for in
  every drawn cell and header of both grids as well as in the view's own text.

  Not claimed here: 1.1's place of the node next to the other apps; 1.4's visual parts (the grid as wide as the view
  without a horizontal scroll bar, the header on one line, the icon colours), the missing tooltip of a disabled icon,
  the tooltips of the enabled Apply model, Delete model and Compare icons (apply-iris.feature 3.4 claims Apply model's,
  catalog.feature the other two), the Applicable
  to input's own icons, the Methods grid's tooltips other than the Roles header, a Method type cell and the Roles cells
  of XGBoost and PLS Regression, and "No models yet." (the stand holds models); 1.3's "Engine appears nowhere" beyond
  the view's text and the two grids' drawn cells and headers (the catalog's rows below the drawn ones, the tooltips
  not hovered); 1.4 step 2's tab cross (the view is closed through the API and claimed gone); 1.6 step 1's table,
  target and features inputs (the Train view is told apart from the old one by its form, and by the old one's model
  preview missing); 1.6 step 3 and 2.9 step 3 (models of the earlier tool are not named); 4.13 step 2's Predicted By
  tag (apply-iris.feature claims the Forge side), step 3 (the table pane, M1c); 6.4's MLFlow node hidden under
  Predictive models; that ML has no Forge group (the package test UI "catalog menu item").

  Where the product moved on from the manual text: ML | Models opens Forge, and the catalog has no menu item; the
  feature opens it with Forge:forgeModels, which focuses an open catalog rather than opening a second one (1.2's "a
  second instance"), so the feature claims that it makes the Forge view current from a table view and that one Forge
  view is open; Methods lists the methods by name (stage 5), not in 1.3's order, and
  lists every installed package's methods, so the feature finds each EDA method by its name and claims neither the
  row count nor a row's position (a stand with Chem lists Chemprop too).

  Background:
    Given user is logged in
    And no Forge training runs of the table "forge-test-iris-{time}" are on the server
    And user opens the demo file "iris.csv" as "forge-test-iris-{time}"
    Then no balloon of any colour should have been shown
    And no errors should have been logged

  Scenario: Forge opens from Apps, from the Predictive models node and from Forge:forgeModels (1.1, 1.2, 6.4, 4.13)
    Then the top menu should list:
      | ML > Models > Train Model... |
      | ML > Models > Apply Model... |
    And the "ML > Models" top menu group should list exactly "Train Model..., Apply Model..."
    Given the browse panel is open
    When user expands "Apps" tree node inside browse tree
    And user clicks on "Apps > Forge" tree node inside browse tree
    Then the "Forge" view should be current
    And "Methods" heading in forge app should be visible
    And "Models" heading in forge app should be visible
    And first grid viewer should be visible
    And second grid viewer should be visible
    When user expands "Platform" tree node inside browse tree
    And user clicks on "Predictive models" tree node inside browse tree
    Then the "Forge" view should be current
    And forge app should be visible
    When user switches to the "forge-test-iris-{time}" table view
    And user calls "Forge:forgeModels" function
    Then the "Forge" view should be current
    And there should be 1 open view named "Forge"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Methods lists the installed methods with their package, type, roles and settings; Engine appears nowhere (1.3)
    Given user switches to the "Forge" view
    Then first grid viewer should show the columns "Method, Package, Method type, Roles, Hyperparameters"
    And the method "CatBoost" should be listed above the method "Decision Tree" in first grid viewer
    And the "Package" cell of the method "CatBoost" in first grid viewer should read "Eda"
    And the "Package" cell of the method "Decision Tree" in first grid viewer should read "Eda"
    And the "Package" cell of the method "Extra Trees" in first grid viewer should read "Eda"
    And the "Package" cell of the method "Hist Gradient Boosting" in first grid viewer should read "Eda"
    And the "Package" cell of the method "k-Nearest Neighbors" in first grid viewer should read "Eda"
    And the "Package" cell of the method "LDA" in first grid viewer should read "Eda"
    And the "Package" cell of the method "Linear Regression" in first grid viewer should read "Eda"
    And the "Package" cell of the method "Naive Bayes" in first grid viewer should read "Eda"
    And the "Package" cell of the method "Neural Network" in first grid viewer should read "Eda"
    And the "Package" cell of the method "PLS Regression" in first grid viewer should read "Eda"
    And the "Package" cell of the method "Random Forest" in first grid viewer should read "Eda"
    And the "Package" cell of the method "RealMLP" in first grid viewer should read "Eda"
    And the "Package" cell of the method "Softmax" in first grid viewer should read "Eda"
    And the "Package" cell of the method "SVM" in first grid viewer should read "Eda"
    And the "Package" cell of the method "TabM" in first grid viewer should read "Eda"
    And the "Package" cell of the method "XGBoost" in first grid viewer should read "Eda"
    And the "Method type" cell of the method "CatBoost" in first grid viewer should read "function"
    And the "Method type" cell of the method "Decision Tree" in first grid viewer should read "function"
    And the "Method type" cell of the method "Extra Trees" in first grid viewer should read "function"
    And the "Method type" cell of the method "Hist Gradient Boosting" in first grid viewer should read "function"
    And the "Method type" cell of the method "k-Nearest Neighbors" in first grid viewer should read "function"
    And the "Method type" cell of the method "LDA" in first grid viewer should read "function"
    And the "Method type" cell of the method "Linear Regression" in first grid viewer should read "function"
    And the "Method type" cell of the method "Naive Bayes" in first grid viewer should read "function"
    And the "Method type" cell of the method "Neural Network" in first grid viewer should read "function"
    And the "Method type" cell of the method "PLS Regression" in first grid viewer should read "function"
    And the "Method type" cell of the method "Random Forest" in first grid viewer should read "function"
    And the "Method type" cell of the method "RealMLP" in first grid viewer should read "function"
    And the "Method type" cell of the method "Softmax" in first grid viewer should read "function"
    And the "Method type" cell of the method "SVM" in first grid viewer should read "function"
    And the "Method type" cell of the method "TabM" in first grid viewer should read "function"
    And the "Method type" cell of the method "XGBoost" in first grid viewer should read "function"
    And the "Roles" cell of the method "XGBoost" in first grid viewer should read "train, apply, isApplicable, isInteractive"
    And the "Hyperparameters" cell of the method "XGBoost" in first grid viewer should read "iterations, eta, maxDepth, lambda, alpha, subsample, colsample, seed"
    And the "Roles" cell of the method "PLS Regression" in first grid viewer should read "train, apply, isApplicable, isInteractive, visualize"
    And the "Hyperparameters" cell of the method "PLS Regression" in first grid viewer should read "components"
    And forge app should not contain text "Engine"
    And no drawn cell or header of first grid viewer should contain "Engine"
    And no drawn cell or header of second grid viewer should contain "Engine"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: The catalog shows its captions and icons; the Methods grid explains itself in tooltips (1.4)
    Given user switches to the "Forge" view
    Then second grid viewer should show the columns "Name, Method, Task, Target, Data storage, Training rows, Tags, Created"
    And "Refresh" icon in forge app should be enabled
    And "Apply model" icon in forge app should be disabled
    And "Delete model" icon in forge app should be disabled
    And "Compare in a new view" icon in forge app should be disabled
    And "Applicable to" input in forge app should be visible
    And "Applicable to" input in forge app should have value ""
    When user hovers over "Refresh" icon in forge app
    Then tooltip should contain text "Refresh"
    When user hovers over the "header Roles" area of first grid viewer
    Then tooltip should contain text "What the method can do."
    When user hovers over the "Method type" cell of the method "Decision Tree" in first grid viewer
    Then tooltip should contain text "A package function."
    When user hovers over the "Roles" cell of the method "XGBoost" in first grid viewer
    Then tooltip should contain text "train: Trains a model on a table and a target column, returns the model"
    And tooltip should contain text "apply: Applies a trained model to a table, returns the predictions"
    And tooltip should contain text "isApplicable: Tells whether the method can learn from the given data"
    And tooltip should contain text "isInteractive: Tells whether training is fast enough to rerun on every change"
    When user hovers over the "Roles" cell of the method "PLS Regression" in first grid viewer
    Then tooltip should contain text "isInteractive: Tells whether training is fast enough to rerun on every change"
    And tooltip should contain text "visualize: Shows method-specific views of a trained model"
    When user closes the current view
    Then forge app should be absent
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: ML | Models opens Forge's Train view and Apply dialog, not the old tool's (1.6, 2.9)
    Given user switches to the "forge-test-iris-{time}" table view
    And user remembers the shown training
    When user picks "ML > Models > Train Model..." from the top menu
    Then the "Predictive model" view should be current
    And a new training's result should be shown
    And train form should be visible
    And model preview should be absent
    When user closes the current view
    And user switches to the "forge-test-iris-{time}" table view
    And user picks "ML > Models > Apply Model..." from the top menu
    Then there should be 1 visible "Apply predictive model" dialog
    And Inputs input in "Apply predictive model" dialog should be absent
    When user clicks on CANCEL button in "Apply predictive model" dialog
    Then the "Apply predictive model" dialog should close
    When user calls "Forge:forgeModels" function
    Then the "Forge" view should be current
    And no error or warning balloon should have been shown
    And no errors should have been logged
