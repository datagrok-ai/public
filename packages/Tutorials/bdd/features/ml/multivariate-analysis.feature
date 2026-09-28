@tutorials @serial @realizes:tutorials.multivariate-analysis
Feature: The Multivariate Analysis tutorial
  Walks Machine Learning > Multivariate Analysis from its card to the end: the PLS dialog, the
  response, fifteen predictors picked through "All" minus price, three components, the names
  column, the run, and the tour of the five charts the analysis docks. Each step is claimed as
  ticked and as done — the analysis adds its columns and charts, and each stop of the tour describes
  the next chart in order.
  Translated from playwright-tests/e2e/tutorials/multivariate-analysis.test.ts, which picked the
  predictors by probing the column picker's canvas pixel by pixel and claimed nothing about the
  analysis. PLS runs in the browser (a WASM web worker), not on the server.

  Fixed in the tutorial for this translation: the counter never reached its last step (7 declared
  for 7 actions); the tour took the charts two seconds after RUN by their position, which a slower
  first run or another docking left short or wrong — it now waits for the charts and finds them by
  title; and its texts were one off from the charts (the Scores text was split in two, so every
  chart after Scores showed the text of the one before).
  Not claimed: which chart a tour stop points at — ui.hints.addHint marks nothing on its anchor.

  Serial, because a finished tutorial writes its completion record into the account's settings,
  which every page syncs whole.

  Background:
    Given user is logged in
    And the "tutorials" user settings are put back at feature end
    And the "achievement-badges" user settings are put back at feature end
    And the "Multivariate Analysis" tutorial is not completed yet
    And the Tutorials app is open

  Scenario: A learner completes the Multivariate Analysis tutorial
    When user starts the "Multivariate Analysis" tutorial
    Then the tutorial progress should be 1 of 8
    When user picks "ML > Analyze > Multivariate Analysis..." from the top menu
    Then the tutorial step "Click on \"ML | Analyze | Multivariate Analysis...\"" should be done
    And "Multivariate Analysis (PLS)" dialog should be visible

    When user selects "price" in Predict input in "Multivariate Analysis (PLS)" dialog
    Then the tutorial step "Set \"Predict\" to \"price\"" should be done
    When user clicks on editor of Using input in "Multivariate Analysis (PLS)" dialog
    Then "Select columns..." dialog should be visible
    When user clicks on "All" link in "Select columns..." dialog
    Then "Select columns..." dialog should contain text "16 checked"
    # price is below the list's fold: the search narrows the list to it
    When user types "price" into "Search" input in "Select columns..." dialog
    And user toggles the "price" column in the column list of "Select columns..." dialog
    Then "Select columns..." dialog should contain text "15 checked"
    When user clicks on OK button in "Select columns..." dialog
    Then the tutorial step "Select all columns, except \"price\", as \"Using\"" should be done
    And editor of Using input in "Multivariate Analysis (PLS)" dialog should contain text "(15)"
    When user enters "3" into Components input in "Multivariate Analysis (PLS)" dialog
    Then the tutorial step "Set the number of components to \"3\"" should be done
    When user selects "model" in Names input in "Multivariate Analysis (PLS)" dialog
    Then the tutorial step "Set \"Names\" to \"model\"" should be done

    When user clicks on RUN button in "Multivariate Analysis (PLS)" dialog
    Then the tutorial step "Click \"RUN\" and wait for the analysis to complete" should be done
    And the table should have a column "price (predicted)"
    And the table should have a column "x.score.t3"
    And table "cars(Features Analysis)" should have 15 rows
    And the open tableview should have 3 scatter plot viewers
    And the open tableview should have 3 bar chart viewers

    Then hint popup should contain text "Observed vs. Predicted"
    When user clicks on "next" button in hint popup
    Then hint popup should contain text "Scores"
    When user clicks on "next" button in hint popup
    Then hint popup should contain text "Loadings"
    When user clicks on "next" button in hint popup
    Then hint popup should contain text "Variable Importance"
    When user clicks on "next" button in hint popup
    Then hint popup should contain text "Explained Variance"
    When user clicks on "done" button in hint popup
    Then the tutorial step "Explore each viewer" should be done

    And the "Multivariate Analysis" tutorial should be completed
    And the tutorial should have listed 7 steps
    And the tutorial progress should be 8 of 8
    And no hint should be shown
    And no errors should have been logged
