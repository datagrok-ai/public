@journey @full-stand @serial @realizes:views.scripts
Feature: Running a script with data from every source
  An R script counts the cells of the table it gets; its Run... dialog takes the table from an open
  table, from a local file, from Datagrok Files and from a database query, and the console runs it
  by its qualified name. A run started from the Scripts view shows nothing in the view: the console
  logs the call and its outputs. Translated from
  files/TestTrack/Scripts/run.md and playwright-public/scripts/scripts-run-debugged.test.ts.

  The script is this feature's own ({time} in its name; the md's shared "testRscript" chain is gone)
  and is deleted with its chats at the end. R runs in a container, so the feature is @full-stand.

  Not translated, and why: nothing of the md is left out. The Activity pane's count of runs is
  claimed in scripts-browser.feature: the context panel renders the current object once, and the
  same script clicked again after a run keeps the count it showed before. The old spec's local-file test passed when
  no file chooser opened at all; here the upload answers the chooser the icon opens, and the count
  proves the file's table was the one run.

  Serial: every scenario here works in the Scripts view, whose search text and view mode are the
  account's own settings — two features searching it at the same time would see each other's text.

  Background:
    Given user is logged in
    And a script "BddScriptRun{time}" is on the server:
      """
      #language: r
      #input: dataframe table
      #output: int count
      #output: string newParam
      count <- nrow(table) * ncol(table)
      newParam <- "test"
      """
    And the context panel is open

  Scenario: Run... with the open cars table counts its 510 cells
    Given user opens cars dataset
    And user opens the Scripts view
    When user clears gallery search
    And user types "BddScriptRun{time}" into gallery search
    And user clicks on "BddScriptRun{time}" link in gallery
    Then the context panel should show "BddScriptRun{time}"
    When user notes the console output
    And user picks "Run..." from the context menu of "BddScriptRun{time}" link in gallery
    Then "BddScriptRun{time}" dialog should be visible
    When user selects "cars" in Table input in "BddScriptRun{time}" dialog
    And user clicks on OK button in "BddScriptRun{time}" dialog
    Then the "BddScriptRun{time}" dialog should close
    And the console should show "count: 510"
    And the console should show "newParam: \"test\""
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A local file is the table
    When user notes the console output
    And user picks "Run..." from the context menu of "BddScriptRun{time}" link in gallery
    Then "BddScriptRun{time}" dialog should be visible
    When user uploads "fixtures/cars-small.csv" through "Open file" icon in "BddScriptRun{time}" dialog
    Then Table input in "BddScriptRun{time}" dialog should have value "cars-small"
    When user clicks on OK button in "BddScriptRun{time}" dialog
    Then the "BddScriptRun{time}" dialog should close
    And the console should show "count: 18"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A file from Datagrok Files is the table
    When user notes the console output
    And user picks "Run..." from the context menu of "BddScriptRun{time}" link in gallery
    Then "BddScriptRun{time}" dialog should be visible
    When user clicks on "Add file from Files" icon in "BddScriptRun{time}" dialog
    Then "Select a file" dialog should be visible
    When user expands "Files > Demo" tree node inside "Select a file" dialog
    And user clicks on "Files > Demo > cars.csv" tree node inside "Select a file" dialog
    And user clicks on OK button in "Select a file" dialog
    Then the "Select a file" dialog should close
    And Table input in "BddScriptRun{time}" dialog should have value "cars"
    When user clicks on OK button in "BddScriptRun{time}" dialog
    Then the "BddScriptRun{time}" dialog should close
    And the console should show "count: 510"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A database query's result is the table
    When user notes the console output
    And user picks "Run..." from the context menu of "BddScriptRun{time}" link in gallery
    Then "BddScriptRun{time}" dialog should be visible
    When user clicks on "Query database" icon in "BddScriptRun{time}" dialog
    Then "Select a database query" dialog should be visible
    When user expands "Postgres" tree node inside "Select a database query" dialog
    And user expands "Postgres > NorthwindTest" tree node inside "Select a database query" dialog
    And user double-clicks on "Postgres > NorthwindTest > PostgresAll" tree node inside "Select a database query" dialog
    Then the "Select a database query" dialog should close
    And Table input in "BddScriptRun{time}" dialog should have value "PostgresAll"
    When user clicks on OK button in "BddScriptRun{time}" dialog
    Then the "BddScriptRun{time}" dialog should close
    And the console should show "count: 11620"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The console runs the script by its qualified name
    When user notes the console output
    And user calls the script "BddScriptRun{time}" from the console with '"cars"'
    Then the console should show "count: 510"
    And the console should show "newParam: \"test\""
    And no errors should have been logged

  Scenario: CANCEL runs nothing
    Given user opens the Scripts view
    When user clears gallery search
    And user types "BddScriptRun{time}" into gallery search
    And user notes the console output
    And user picks "Run..." from the context menu of "BddScriptRun{time}" link in gallery
    Then "BddScriptRun{time}" dialog should be visible
    When user clicks on CANCEL button in "BddScriptRun{time}" dialog
    Then the "BddScriptRun{time}" dialog should close
    # the run after it is the positive the claim needs: one call logged, not two
    When user picks "Run..." from the context menu of "BddScriptRun{time}" link in gallery
    And user selects "cars" in Table input in "BddScriptRun{time}" dialog
    And user clicks on OK button in "BddScriptRun{time}" dialog
    Then the "BddScriptRun{time}" dialog should close
    And the console should show "count: 510" 1 time
    And no errors should have been logged
    # the gallery keeps its search text for the next visit: leave it empty
    When user clears gallery search
