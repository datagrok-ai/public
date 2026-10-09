@journey @eda @realizes:views.models @realizes:ml.menu.models.train-model @realizes:ml.menu.models.apply-model
Feature: The Predictive models gallery
  Browse > Platform > Predictive models: finding a model by its name and through the quick filters,
  editing its description, re-running its evaluation, applying it to an open table from its card,
  and saving it as a zip. Translated from files/TestTrack/Models/models-testdemog-lifecycle-smoke.md
  block 4 and models-lifecycle-csv-table.md scenarios 3 and 4, and their playwright-public specs,
  which reached the gallery by its route, checked the Run Evaluation click only for "no new error"
  and the Filter templates icon for "some section exists".

  Two models of the feature's own, each trained on a table written into it: one predicts the class of
  "readings" from f1 and f2, the other the score of "levels" from g1 and g2, so each is applicable to
  its own table only, which its card says ("Applicable to levels"). The gallery search is fuzzy:
  "Kestrel" keeps one of the two models and not the other, and the counter drops; Created by me writes
  "author = @current" into an empty search (a typed search is kept and combined with it). Run Evaluation applies the model to the table it was trained
  on and draws its result charts and metrics (Accuracy, Confusions, a scatter plot) in the pane. Apply to on a card opens the Apply dialog on that table and that model, and makes the table's view
  current once the prediction is added. The dialog makes its model the current object too, as a walk
  by hand shows; under the test the context panel kept the model clicked before (a guard of the panel
  dropping the change, not instrumented), so the scenario claims the dialog's own Model input instead.

  Not translated, and why: picking two cards with Ctrl+click and comparing them (Actions > Compare
  opens "Compare models" with Name, Description, Method and Source) — no step clicks an element with a
  key held (MISSING.md); the Is applicable to... quick filter, which saves the chosen open table to
  the server (no step removes such a table) and, on this stand, kept both models for a table only one of
  them fits (MISSING.md, to check by hand); the Context Panel "walk every tab" of the smoke case, which is
  apply-and-delete's Details / Performance / Sharing claim; Share..., which share-model walks; Delete,
  which apply-and-delete walks — a model deleted through the gallery leaves the table its save uploaded
  and its wrapper project on the server (MISSING.md), so the models here go by name, now and at feature
  end, through the sweep that removes those too.

  Background:
    Given user is logged in
    And no predictive model named "BDD-Kestrel-{run}" is on the server
    And no predictive model named "BDD-Osprey-{run}" is on the server

  Scenario: Two models are trained and saved, each on a table of its own
    Given user opens a table "readings" with:
      | f1 | f2 | target |
      | 1  | 7  | A      |
      | 2  | 3  | A      |
      | 3  | 9  | B      |
      | 4  | 1  | A      |
      | 5  | 8  | B      |
      | 6  | 2  | A      |
      | 7  | 10 | B      |
      | 8  | 4  | B      |
      | 9  | 6  | A      |
      | 10 | 5  | B      |
    When user picks "ML > Models > Train Model..." from the top menu
    And user selects "target" in Predict input
    And user clicks on editor of Features input
    And user clicks on None label in "Select columns..." dialog
    And user clicks on the "cell 1 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on the "cell 2 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then model preview should be ready
    When user clicks on Save button
    And user enters "BDD-Kestrel-{run}" into Name input in dialog
    And user clicks on OK button in dialog
    Then dialog should be absent
    And 1 predictive model named "BDD-Kestrel-{run}" should be on the server
    Given user opens a table "levels" with:
      | g1 | g2  | score |
      | 1  | 1.1 | 3.1   |
      | 2  | 2.3 | 1.2   |
      | 3  | 2.9 | 4.4   |
      | 4  | 4.2 | 2.5   |
      | 5  | 5.1 | 5.3   |
      | 6  | 6.4 | 1.7   |
      | 7  | 6.8 | 3.9   |
      | 8  | 8.3 | 2.2   |
    When user picks "ML > Models > Train Model..." from the top menu
    And user selects "score" in Predict input
    And user clicks on editor of Features input
    And user clicks on None label in "Select columns..." dialog
    And user clicks on the "cell 1 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on the "cell 2 of x" area of grid viewer in "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then model preview should be ready
    When user clicks on Save button
    And user enters "BDD-Osprey-{run}" into Name input in dialog
    And user clicks on OK button in dialog
    Then dialog should be absent
    And 1 predictive model named "BDD-Osprey-{run}" should be on the server
    And no errors should have been logged

  Scenario: The gallery search finds a model by a word of its name
    Given the context panel is open
    And the browse panel is open
    And Platform tree node inside browse tree is expanded
    When user clicks on "Predictive models" tree node inside browse tree
    Then the "Models" view should be current
    And "BDD-Kestrel-{run}" label in gallery should be visible
    And "BDD-Osprey-{run}" label in gallery should be visible
    When user remembers the gallery counter
    And user types "Kestrel" into gallery search
    Then the gallery counter should be lower than remembered
    And "BDD-Kestrel-{run}" label in gallery should be visible
    And "BDD-Osprey-{run}" label in gallery should be absent
    When user types "Osprey" into gallery search
    Then "BDD-Osprey-{run}" label in gallery should be visible
    And the gallery counter should be lower than remembered
    And "BDD-Kestrel-{run}" label in gallery should be absent
    # other features save and delete models at the same time, so the counter is claimed only under a search
    When user clears gallery search
    Then "BDD-Kestrel-{run}" label in gallery should be visible
    And no errors should have been logged

  Scenario: The Created by me quick filter writes its query into the search
    When user clicks on "Toggle filters" icon
    Then "Is applicable to..." tag should be visible
    When user clicks on "Created by me" tag
    Then gallery search should have value "author = @current"
    When user clicks on "All" tag
    Then gallery search should have value ""
    When user clicks on "Toggle filters" icon
    Then "Is applicable to..." tag should be hidden
    And no errors should have been logged

  Scenario: Edit... changes the description the Details pane shows
    When user picks "Edit..." from the context menu of "BDD-Kestrel-{run}" label in gallery
    Then "Predictive model" dialog should be visible
    And Name input in "Predictive model" dialog should have value "BDD-Kestrel-{run}"
    When user enters "Classifies the readings, {run}" into Description input in "Predictive model" dialog
    And user clicks on OK button in "Predictive model" dialog
    Then the "Predictive model" dialog should close
    When user clicks on "BDD-Kestrel-{run}" label in gallery
    Then the context panel should show "BDD-Kestrel-{run}"
    And "Details" pane in context panel should contain text "Classifies the readings, {run}"
    And "Details" pane in context panel should contain text "f1, f2"
    And no errors should have been logged

  Scenario: Run Evaluation draws the model's charts and metrics on its training table
    When user clicks on "Performance" pane in context panel
    Then "Run Evaluation" button in context panel should be visible
    And "Performance" pane in context panel should not contain text "Accuracy"
    When user clicks on "Run Evaluation" button in context panel
    Then "Performance" pane in context panel should contain text "Accuracy"
    And "Performance" pane in context panel should contain text "Confusions"
    And scatter plot viewer in "Performance" pane in context panel should be visible
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: A card applies its model to an open table it fits, and only to that one
    Then "BDD-Osprey-{run}" gallery card should contain text "Applicable to levels"
    And "BDD-Osprey-{run}" gallery card should not contain text "readings"
    And "BDD-Kestrel-{run}" gallery card should contain text "Applicable to readings"
    And "BDD-Kestrel-{run}" gallery card should not contain text "levels"
    When user picks "Apply to > levels (8 rows, 3 columns)" from the context menu of "BDD-Osprey-{run}" label in gallery
    Then "Apply predictive model" dialog should be visible
    And Model input in "Apply predictive model" dialog should contain text "BDD-Osprey-"
    And Model input in "Apply predictive model" dialog should not contain text "BDD-Kestrel-"
    And Inputs input in "Apply predictive model" dialog should contain text "(2/2)"
    When user clicks on OK button in "Apply predictive model" dialog
    Then the "Apply predictive model" dialog should close
    And the "levels" view should be current
    And the table should have 4 columns
    And the newest column matching "^score \(" should have no missing values
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Save as Zip downloads the model
    Given user watches downloads
    And user switches to the "Models" view
    When user picks "Save as Zip" from the context menu of "BDD-Kestrel-{run}" label in gallery
    Then a file matching "\.zip$" should have been downloaded
    And no errors should have been logged
