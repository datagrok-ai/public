@journey @full-stand @serial @realizes:views.scripts
Feature: The Octave template counts the cells of its sample table
  The Octave template of the New menu takes the sample table and answers the number of its cells,
  like every other language's (files/TestTrack/Scripts/create.md 13-20). With cars — 30 rows by 17
  columns, 510 cells, which R, Python, NodeJS, Julia, Grok and Pyodide answer — Octave answers 527:
  the table reaches Octave with 31 rows (probed on dev 22 Sep through the JS API as well, so it is
  the Octave handler's conversion, not the editor). That is GROK-17456 "Octave: Incorrect rows
  count", closed Won't fix, so 527 is what the product promises today and what the scenario claims;
  should the handler ever be fixed, this scenario is the one that says so.

  Octave runs in a container, so the feature is @full-stand. The run and its answer are two
  scenarios: the first waits for the run to end, so a failure of the second is the value, not a
  timeout.

  Not translated, and why: nothing of the Octave row of the md is left out.

  Serial: every scenario here works in the Scripts view, whose search text and view mode are the
  account's own settings — two features searching it at the same time would see each other's text.

  Background:
    Given user is logged in
    And user opens the Scripts view

  Scenario: The Octave template runs with cars
    When user clicks on New button
    And user picks "Octave Script..." from the open menu
    Then the "Template" view should be current
    And code editor should contain the text "#language: octave"
    When user clicks on "Open script sample table" icon
    Then table "cars" should be open
    When user clicks on "Run script (F5)" icon
    Then "Template" dialog should be visible
    When user selects "cars" in Table input in "Template" dialog
    And user clicks on OK button in "Template" dialog
    Then the "Template" dialog should close
    And the script results should list "count"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Octave counts the cells of cars its own way
    Then the script results should show "count" as "527"
