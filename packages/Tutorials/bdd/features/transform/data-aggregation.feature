@tutorials @serial @realizes:tutorials.data-aggregation
Feature: The Data Aggregation tutorial
  Walks Data Transformation > Data Aggregation from its card to the end: the Aggregation Editor, two
  group-by columns and a pivot, the aggregation narrowed, repointed and switched to a median, the
  aggregated table selecting and filtering its source, and the parameters saved to the history.
  Each step is claimed as ticked and as done — the editor's configuration as the pivot computed it,
  37 selected rows for "Asian, F"
  and 75 filtered rows for "Other, M", the saved entry in the history menu.
  Translated from playwright-tests/e2e/tutorials/data-aggregation.test.ts, which reached the
  aggregated rows by pixel sweeps over a guessed canvas and claimed no configuration at all.

  The saved parameters live in the browser's storage and are cleared before the walk, so the entry
  the last step saves is the only one.

  Serial, because a finished tutorial writes its completion record into the account's settings,
  which every page syncs whole.

  Background:
    Given user is logged in
    And the "tutorials" user settings are put back at feature end
    And the "achievement-badges" user settings are put back at feature end
    And the "Data Aggregation" tutorial is not completed yet
    And user clears the saved pivot table parameters
    And the Tutorials app is open

  Scenario: A learner completes the Data Aggregation tutorial
    When user starts the "Data Aggregation" tutorial
    Then the tutorial progress should be 1 of 11
    When user picks "Data > Aggregate Rows..." from the top menu
    Then the tutorial step "Open Aggregation Editor" should be done
    And the open tableview should have 1 pivot table viewer
    # the tutorial resets the editor to avg(AGE) and avg(HEIGHT) with nothing grouped
    And the "aggregate" reading of pivot table viewer should be "avg(AGE), avg(HEIGHT)"

    When user adds "RACE" to the "group by" row of pivot table viewer
    Then the tutorial step "Group rows by column \"RACE\"" should be done
    When user adds "SEX" to the "group by" row of pivot table viewer
    Then the tutorial step "Group rows by column \"SEX\"" should be done
    And the "group by" reading of pivot table viewer should be "RACE, SEX"
    When user adds "DIS_POP" to the "pivot" row of pivot table viewer
    Then the tutorial step "Pivot data by column \"DIS_POP\"" should be done
    And the "pivot" reading of pivot table viewer should be "DIS_POP"

    When user picks "Remove others" from the context menu of the "aggregate chip avg(AGE)" area of pivot table viewer
    Then the tutorial step "Leave only the \"avg(AGE)\" aggregation" should be done
    And the "aggregate" reading of pivot table viewer should be "avg(AGE)"
    When user picks "Column > WEIGHT" from the context menu of the "aggregate chip avg(AGE)" area of pivot table viewer
    And user closes the context menu
    Then the tutorial step "Change a column to \"WEIGHT\"" should be done
    And the "aggregate" reading of pivot table viewer should be "avg(WEIGHT)"
    When user picks "Aggregation > med" from the context menu of the "aggregate chip avg(WEIGHT)" area of pivot table viewer
    And user closes the context menu
    Then the tutorial step "Change the aggregation function to \"med\"" should be done
    And the "aggregate" reading of pivot table viewer should be "med(WEIGHT)"

    # the first aggregated row is "Asian, F": 37 patients
    When user clicks on the "grid row header 1" area of pivot table viewer holding Shift
    Then the tutorial step "Select rows in the source table with values of the first aggregated row" should be done
    And 37 rows should be selected
    And every selected row should pass the filter
    And no rows where "RACE" is "Caucasian" should be selected
    And no rows where "RACE" is "Other" should be selected
    And no rows where "SEX" is "M" should be selected
    When user presses Escape
    Then the tutorial step "Remove selection by pressing \"Esc\"" should be done
    And no rows should be selected

    # the last of the 8 aggregated rows is "Other, M": 75 patients
    When user clicks on the "grid cell 8 of RACE" area of pivot table viewer
    Then the tutorial step "Click on the last row in the aggregated table to filter by it" should be done
    And 75 rows should pass the filter
    And no rows where "RACE" is "Asian" should pass the filter
    And no rows where "SEX" is "F" should pass the filter

    When user picks "Save parameters" from the history menu of pivot table viewer
    Then the tutorial step "Save parameters" should be done
    And the "history entries" reading of pivot table viewer should contain "med(WEIGHT)"

    And the "Data Aggregation" tutorial should be completed
    And the tutorial should have listed 11 steps
    And the tutorial progress should be 11 of 11
    And no hint should be shown
    And no errors should have been logged
