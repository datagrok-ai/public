@journey @realizes:views.functions
Feature: Run history records every run of the function view and brings any back
  Every RUN in the function view is recorded on the server with its own number, shown in the header
  above the results. The ribbon's history icon docks the Run history panel in the right column: the
  function's runs, newest first, each with what changed against the run before it, a search and the
  Mine / Starred / Failed filters. The pencil names a run, the star keeps it, a click on a row brings the
  run's inputs back at its own address, and the author can delete a run. The panel leaves with its view
  and comes back with it. A function that reruns on every input change records one run once the inputs
  rest, Save to space links a run into a space (the Space row is PowerPack's), and a package function's
  history lists its runs (Compute2's ObjectCooling2; only the runs this feature makes are deleted).

  Two JavaScript scripts and a space are made through the API; the scripts, every run they record and
  the space are deleted at feature start and end.

  Background:
    Given user is logged in
    And simple mode is off
    And the runs of the "BDDRunHistory{time}" function are deleted at the end
    And the runs of the "BDDRunLive{time}" function are deleted at the end
    And a script "BDDRunHistory{time}" is on the server:
      """
      //language: javascript
      //input: int rows = 5 {caption: Rows}
      //output: int total {caption: Total}
      total = rows * 10;
      """
    And a script "BDDRunLive{time}" is on the server:
      """
      //language: javascript
      //input: int rows = 5 {caption: Rows}
      //output: int total {caption: Total}
      //meta.runOnInput: true
      total = rows * 10;
      """
    And a space named "BDDRunSpace{time}" is on the server

  Scenario: Every RUN is recorded with its own number
    When user opens the address "/func/Admin.BDDRunHistory{time}?rows=5"
    And user clicks on RUN button
    Then run header should contain text "Run #1"
    When user enters "7" into "Rows" input
    And user clicks on RUN button
    Then run header should contain text "Run #2"
    And "Total" input should have value "70"

  Scenario: The history icon docks the runs, each with what changed against the run before it
    When user clicks on "Run history" icon
    Then run history should be visible
    And "Run #2" run should be visible
    And "Run #1" run should be visible
    And changes of "Run #2" run should contain text "rows 5 → 7"
    And changes of "Run #2" run should contain text "total = 70"

  Scenario: The pencil names a run, and the name shows in the header and in the history
    When user clicks on "Edit title, description and tags" icon in run header
    Then "Run #2" dialog should be visible
    When user enters "Baseline" into "Title" input in "Run #2" dialog
    And user clicks on OK button in "Run #2" dialog
    Then run header should contain text "Baseline"
    And "Baseline" run should be visible

  Scenario: The star keeps a run, and the Starred filter shows only the starred ones
    When user clicks on run header star
    And user clicks on "Starred" run filter
    Then there should be 1 visible run in run history
    And "Baseline" run should be visible
    When user clicks on "Starred" run filter
    Then there should be 2 visible run in run history

  Scenario: The search finds a run by its name
    When user types "Baseline" into run history search
    Then there should be 1 visible run in run history
    When user clears run history search
    Then there should be 2 visible run in run history

  Scenario: A click on a run brings its inputs back at the run's own address
    When user clicks on "Run #1" run
    Then "Rows" input should have value "5"
    And run header should contain text "Run #1"
    And the page address should contain "/runs/1"

  Scenario: The author deletes a run of their own
    When user hovers over "Run #1" run
    And user clicks on delete of "Run #1" run
    Then "Are you sure?" dialog should be visible
    When user clicks on DELETE button in "Are you sure?" dialog
    Then "Run #1" run should not be visible
    And "Baseline" run should be visible

  Scenario: The panel leaves with its view and comes back with it
    Given user opens demog-1000 dataset
    Then run history should not be visible
    When user clicks on "BDDRunHistory{time}" view tab
    Then run history should be visible

  Scenario: Save to space links the run into the space picked
    Given the "PowerPack" package is installed
    When user clicks on "BDDRunHistory{time}" view tab
    And user clicks on "Baseline" run
    Then run header should contain text "Baseline"
    When user clicks on "Save this run to a space" icon in run header
    Then "Save Baseline to a space" dialog should be visible
    And space input in "Save Baseline to a space" dialog should contain text "Your own projects"
    When user clicks on space input in "Save Baseline to a space" dialog
    And user clicks on "BDDRunSpace{time}" tree node
    And user clicks on OK button in "Save Baseline to a space" dialog
    Then the "Save Baseline to a space" dialog should close
    And the space "BDDRunSpace{time}" on the server should link a run of the "BDDRunHistory{time}" function

  Scenario: A function that reruns on input records one run once the inputs rest
    When user opens the address "/func/Admin.BDDRunLive{time}?rows=5"
    And user enters "6" into "Rows" input
    And user enters "8" into "Rows" input
    Then "Total" input should have value "80"
    And run header should contain text "Run #"
    When user clicks on "Run history" icon
    Then "Run #1" run should be visible
    And changes of "Run #1" run should contain text "rows 8"

  Scenario: A package function's history lists its runs
    Given the "Compute2" package is installed
    And the runs of the "ObjectCooling2" function made from now on are deleted at the end
    When user opens the address "/func/Compute2.ObjectCooling2"
    And user clicks on RUN button
    Then run header should contain text "Run #"
    When user clicks on "Run history" icon
    Then run history should be visible
    And first run in run history should be visible
