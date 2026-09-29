@tutorials @serial @realizes:tutorials.scripting
Feature: The Scripting tutorial
  Walks Machine learning > Scripting from its card to the end: a Python script editor opened from the
  Browse tree, its sample table opened, the script run for that table, the console opened, a second
  output added to the script and the script run again. Each step is claimed as ticked and as done —
  the editor, the table, the runs and their outputs.
  Translated from playwright-tests/e2e/tutorials/scripting.test.ts.

  The script runs in Jupyter on the server, and the tutorial will not start without it (its
  prerequisite): where the stand does not run it, the walk is skipped, not failed.
  Fixed in the tutorial for this translation: the counter was 11 for 10 steps.
  Green on dev (Jupyter running). The outputs are shown under the editor, in "Results"
  (count 510), later a "Clone" grid; the console the tutorial sends the learner to shows none
  of them, and no table view opens for the dataframe output, although the tutorial text promises
  both — the claims are on what the editor shows.

  Serial, because a finished tutorial writes its completion record into the account's settings,
  which every page syncs whole.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the "tutorials" user settings are put back at feature end
    And the "achievement-badges" user settings are put back at feature end
    And the "Scripting" tutorial is not completed yet
    And the Tutorials app is open

  Scenario: A learner completes the Scripting tutorial
    # the tutorial will not start without Jupyter (its prerequisite), and its script runs there
    Given the stand runs the "Jupyter" service
    When user starts the "Scripting" tutorial
    Then the tutorial progress should be 1 of 10
    Given the tutorial step "In the Browse Panel, click Platform > Functions > Scripts > New > Python Script. This opens a script editor." should not be done yet
    When user clicks on Platform---Functions---Scripts tree node inside browse tree
    And user clicks on NEW button
    And user picks "Python Script" from the open menu
    Then the tutorial step "In the Browse Panel, click Platform > Functions > Scripts > New > Python Script. This opens a script editor." should be done
    And code editor should be visible
    Given the tutorial step "Open a sample table for the script" should not be done yet
    When user clicks on asterisk icon
    Then the tutorial step "Open a sample table for the script" should be done
    Given the tutorial step "Run the script" should not be done yet
    When user clicks on play icon
    Then the tutorial step "Run the script" should be done
    And "Template" dialog should be visible
    When user selects "cars" in "Table" input in "Template" dialog
    Then the tutorial step "Set \"Table\" to cars" should be done
    When user clicks on OK button in "Template" dialog
    # the dialog stays open while the script runs in Jupyter, past the 15 s of a claim on a busy stand
    Then the "Template" dialog should close
    And the tutorial step "Click \"OK\"" should be done
    # the scalar output is shown under the editor, not in the console the next step opens (see the description)
    And "Results" dock panel should contain text "510"
    Given the tutorial step "Find the results in the console" should not be done yet
    # "You can use the ~ key to control the console visibility"
    When user opens the console
    Then the tutorial step "Find the results in the console" should be done

    Given the tutorial step "Add the second output value to the script" should not be done yet
    When user types "#output: dataframe clone" into the first empty line of code editor
    And user types "clone = table" into the last line of code editor
    Then the tutorial step "Add the second output value to the script" should be done
    Given the tutorial step "Run the script" should be listed 2 times
    When user clicks on play icon
    Then the tutorial step "Run the script" should be done 2 times
    When user selects "cars" in "Table" input in "Template" dialog
    Then the tutorial step "Set \"Table\" to cars" should be done 2 times
    When user clicks on OK button in "Template" dialog
    Then the "Template" dialog should close
    And the tutorial step "Click \"OK\"" should be done 2 times
    And "Clone" dock panel should be visible

    And the "Scripting" tutorial should be completed
    And the tutorial should have listed 10 steps
    And the tutorial progress should be 10 of 10
    And no hint should be shown
    And no errors should have been logged
