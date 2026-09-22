@journey @full-stand @serial @realizes:views.scripts
Feature: Creating a script
  A new R script from the Scripts view: the New menu offers every language, the template opens in
  the editor unsaved, the sample icon brings its sample table, the Signature editor renames the
  script and adds a parameter that lands in the code, a run with cars answers count = 510 in the
  results under the editor, and Save stores the script with the parameters the server parses from
  its header. Translated from files/TestTrack/Scripts/create.md (1-12) and
  playwright-public/scripts/scripts-create-debugged.test.ts (the R test).

  R runs in a container, so the feature is @full-stand (dev has it). Every name carries {time}: the
  old chain of files sharing one "testRscript" is gone, and the create feature renames the template
  before it saves, so it never adds another of the old "Template_N" scripts.

  Not translated, and why: the md's new parameter is an output string; the Signature editor's
  parameter grid adds one (as an input bool), but its Direction and Type cells open no editor on a
  click, a double-click or Enter (probed on dev 22 Sep), so the feature claims what the "+" icon
  does and the header it writes. Step 12's "x" is the view's own close, which the harness does
  through the shell — several "Close view" icons are on the page at once.

  Serial: every scenario here works in the Scripts view, whose search text and view mode are the
  account's own settings — two features searching it at the same time would see each other's text.

  Background:
    Given user is logged in
    And no script named "BddScriptCreate{time}" is on the server
    And the browse panel is open

  Scenario: The Scripts view opens from Platform > Functions
    When user expands "Platform" tree node inside browse tree
    And user expands "Platform > Functions" tree node inside browse tree
    And user clicks on "Platform > Functions > Scripts" tree node inside browse tree
    Then the "Scripts" view should be current
    And gallery should be visible
    And no errors should have been logged

  Scenario: New offers every language and opens the R template unsaved
    When user clicks on New button
    Then the open menu should list "R Script..."
    And the open menu should list "Python Script..."
    And the open menu should list "Octave Script..."
    And the open menu should list "NodeJS Script..."
    And the open menu should list "Julia Script..."
    And the open menu should list "JavaScript Script..."
    And the open menu should list "Grok Script..."
    And the open menu should list "Pyodide Script..."
    When user picks "R Script..." from the open menu
    Then the "Template" view should be current
    And code editor should contain the text "#language: r"
    And code editor should be visible
    And no errors should have been logged

  Scenario: The sample icon opens the sample table
    When user clicks on "Open script sample table" icon
    Then table "cars" should be open
    And table "cars" should have 30 rows
    And the "Template" view should be current
    And no errors should have been logged

  Scenario: The Signature editor names the script and adds a parameter to its header
    Given user switches to the "Template" view
    When user clicks on "Open Signature Editor" icon
    Then PARAMETERS tab should be visible
    When user enters "BddScriptCreate{time}" into Name input
    And user clicks on PARAMETERS tab
    Then there should be 2 visible "Add the param" icon
    When user clicks on first "Add the param" icon
    Then there should be 3 visible "Add the param" icon
    Then "Open function editor" icon should be visible
    When user clicks on "Open function editor" icon
    Then "Run script (F5)" icon should be visible
    And code editor should contain the text "#name: BddScriptCreate{time}"
    And code editor should contain the text "#input: bool newParam"
    And no errors should have been logged

  Scenario: Running with cars answers the number of cells
    # the run dialog still carries the template's name: the editor re-reads the header on Save
    Given user switches to the "Template" view
    When user clicks on "Run script (F5)" icon
    Then "Template" dialog should be visible
    When user selects "cars" in Table input in "Template" dialog
    And user clicks on OK button in "Template" dialog
    Then the "Template" dialog should close
    And the script results should show "count" as "510"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Save stores the script with the parameters of its header
    When user saves the script
    Then an info balloon containing "Script saved." should have been shown
    And 1 script named "BddScriptCreate{time}" should be on the server
    And the script "BddScriptCreate{time}" on the server should have an output "count" of type "int"
    And the script "BddScriptCreate{time}" on the server should have an input "newParam" of type "bool"
    And the script "BddScriptCreate{time}" on the server should have an input "table" of type "dataframe"
    And the "BddScriptCreate{time}" view should be current
    And no errors should have been logged

  Scenario: Closing the editor returns to Scripts, where the new script is found
    When user closes the current view
    Then the "Scripts" view should be current
    When user remembers the gallery counter
    And user types "BddScriptCreate{time}" into gallery search
    Then the gallery counter should be lower than remembered
    And "BddScriptCreate{time}" link in gallery should be visible
    And no errors should have been logged
