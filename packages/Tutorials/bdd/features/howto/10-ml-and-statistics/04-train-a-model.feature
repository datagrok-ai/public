@guide @tutorial
Feature: How-to · ML
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). Train, save, then apply. ML > Models > Train Model... opens the
  Forge package's Train view: on iris it takes Species (the last column) as the target and the four measurements as
  features, suggests XGBoost and trains at once in the browser (EDA's Crux methods, no server script), with 5-fold
  cross-validation. The clip points at Performance (5-fold cross-validation scores) and Confusions (the model trained
  on all rows, predicting those same rows: all 150 right, while Performance's validation accuracy is 0.927, so the
  line says which rows it shows rather than "which species it mixes up"); another method (Random Forest) retrains by
  itself and Runs lists both trainings. Save names the model "Iris species" and keeps no data (Data storage None);
  Browse > Platform > Machine Learning (shown when Forge is installed) opens the Models view with its card. Then the
  view tab of "new flowers" (nine iris rows without Species: rows 1, 51, 101, 11, 61, 111, 21, 81, 121, three of each
  species, interleaved) and ML > Models > Apply Model...: the dialog presets the newest model that fits the table and
  matches each of its inputs to the column of the same name; OK predicts in the browser and adds "Species
  (predicted)". Both tables are opened off camera, "new flowers" first so that "iris flowers" is current. Every
  training is recorded as a Forge training run of the table, so iris is opened as "iris flowers": the model, its
  applications and the runs of that table are deleted at feature start and end (howto-ml.ts, through
  Forge:deleteModel and grok.dapi.domains). Merged from the former clips "Train a model" and "Apply a model to new
  data".

  Scenario: Train a model and apply it
    Train a model that predicts one column from the others, save it, and fill in its predictions for new rows.

    Given user is logged in
    And simple mode is off
    And the Forge model "Iris species" and the training runs of the table "iris flowers" are removed now and when the feature ends
    And user opens a table "new flowers" with:
      | Sepal.Length | Sepal.Width | Petal.Length | Petal.Width |
      | 5.1          | 3.5         | 1.4          | 0.2         |
      | 7.0          | 3.2         | 4.7          | 1.4         |
      | 6.3          | 3.3         | 6.0          | 2.5         |
      | 5.4          | 3.7         | 1.5          | 0.2         |
      | 5.0          | 2.0         | 3.5          | 1.0         |
      | 6.5          | 3.2         | 5.1          | 2.0         |
      | 5.4          | 3.4         | 1.7          | 0.2         |
      | 5.5          | 2.4         | 3.8          | 1.1         |
      | 6.9          | 3.2         | 5.7          | 2.3         |
    And user opens iris dataset keeping the first 150 rows as "iris flowers"
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    When the video says "150 flowers: four measurements and the species" at the top
    # caption: Train Model is under ML > Models
    And user picks "ML > Models > Train Model..." from the top menu
    Then the Train view should show a finished training
    # caption: Target: the column to predict
    When user hovers over editor of Target input
    # caption: Features: the columns it learns from
    And user hovers over editor of Features input
    # caption: Method: XGBoost here; training runs in your browser
    And user hovers over editor of Method input
    # caption: It trains right away. Performance: scores on held-out rows
    And user hovers over title of "Performance" results panel
    Then "Performance" results panel should contain text "5-fold cross-validation"
    # caption: Confusions: predicted against actual, on the rows it learned from
    When user hovers over title of "Confusions" results panel
    # caption: Pick another method and it trains again
    And user remembers the training the Train view shows
    And user selects "Random Forest" in Method input
    Then the Train view should show a new training
    # caption: Runs: every training so far, to compare
    When user hovers over title of "Runs" results panel
    # caption: Save keeps the model under a name
    When user clicks on Save button
    Then "Save model" dialog should be visible
    When user enters "Iris species" into Name input in "Save model" dialog
    And user clicks on OK button in "Save model" dialog
    Then the "Save model" dialog should close
    # caption: Saved models are in Browse > Platform > Machine Learning
    When user clicks on browse tab
    And user expands "Platform" tree node inside browse tree
    And user clicks on "Machine Learning" tree node inside browse tree
    Then the current view should show the text "Iris species"
    # caption: Your model's card: what it predicts, from what, and how
    When user hovers over the text "Iris species" in gallery
    # caption: Now apply it: a table of nine new flowers
    And user clicks on the tab of the "new flowers" view
    Then the "new flowers" view should be current
    When the video says "Measured, but the species is not known yet" at the bottom
    # caption: Apply Model is under ML > Models
    And user picks "ML > Models > Apply Model..." from the top menu
    Then "Apply predictive model" dialog should be visible
    # caption: Model: the one just saved is already picked
    When user hovers over Model input in "Apply predictive model" dialog
    Then Model input in "Apply predictive model" dialog should have value "Iris species"
    # caption: Columns: each input of the model gets the column of the same name
    When user expands Columns collapsible group in "Apply predictive model" dialog
    And user hovers over editor of "Petal.Length" input in "Apply predictive model" dialog
    # caption: OK adds the prediction as a new column
    And user clicks on OK button in "Apply predictive model" dialog
    Then "Apply predictive model" dialog should be hidden
    And the table should have a column "Species (predicted)"
    And "Species (predicted)" column should have no missing values
    And "Species (predicted)" column should have at least 3 distinct values
    # caption: Each flower gets the species the model expects from its measurements
    When user hovers over the "cell 2 of Species (predicted)" area of grid
    Then no error or warning balloon should have been shown
    And no errors should have been logged
