@tutorials @serial @realizes:tutorials.scripting
Feature: The Scripting tutorial
  Walks Machine learning > Scripting from its card to its first run: a Python script editor opened
  from the Browse tree, its sample table opened, the script run for that table. Each step is claimed
  as ticked and as done — the editor, the table, the run dialog and its end.
  Translated from playwright-tests/e2e/tutorials/scripting.test.ts.

  The script runs in Jupyter on the server, and the tutorial will not start without it (its
  prerequisite): where the stand does not run it, the walk is skipped, not failed.
  Not translated: the rest of the tutorial (the console, a second output typed into the script, the
  second run). Each of those steps waits for a Python result, and what a server-side script returns is
  not a feature's claim (the library's hard rule); the script belongs to a package test. Walked by
  hand on dev: the outputs show under the editor ("Results" count 510, then a "Clone" grid), not in
  the console the tutorial sends the learner to.
  Fixed in the tutorial for this translation: the counter was 11 for 10 steps.

  Serial, because a finished tutorial writes its completion record into the account's settings,
  which every page syncs whole.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the "tutorials" user settings are put back at feature end
    And the "achievement-badges" user settings are put back at feature end
    And the "Scripting" tutorial is not completed yet
    And the Tutorials app is open

  Scenario: A learner runs the Scripting tutorial's script for the first time
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
    And no errors should have been logged
