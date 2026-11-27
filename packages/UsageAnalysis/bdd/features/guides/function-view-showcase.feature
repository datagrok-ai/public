@guide @demo
Feature: Function view · What's new
  A demo clip set for the core team (filmed with `grok-bdd guide`: the @demo tag leaves the set-up and the checks off
  camera), one scenario per clip, cut into one reel with `guide-render.py --reel`. Compute2's ObjectCooling2 model
  (Newton's law of cooling) opens in the core function view at /func/Compute2.ObjectCooling2: the form above the
  results with the viewers the script declares, the "Inputs changed" note and a re-run in place, the link with the
  current inputs, the Excel / PowerPoint exports, fitting with its fitted coefficient applied back to the form, the
  Run history panel (name, star, bring an older run back) and Save to space. ObjectCooling2 is a shared model, so
  only the runs this feature makes are deleted at its end, by id; the space is a fixture deleted with it.

  Background:
    Given user is logged in
    And simple mode is off
    And the "Compute2" package is installed
    And the runs of the "ObjectCooling2" function made from now on are deleted at the end
    And the browse panel is closed
    And the context panel is closed

  Scenario: A model opens as a form with its results
    The form on top, the results below in the viewers the model declares.

    When user opens the address "/func/Compute2.ObjectCooling2"
    Then "Heat transfer coefficient" input should have value "8.3"
    When the video says "Any function is a form, and its results"
    And user clicks on RUN button
    Then "Line chart" viewer should be visible
    And "Grid" viewer should be visible
    And "Time to cool" input should be visible

  Scenario: Change an input, run again in place
    A changed input marks the results stale; RUN updates them in the same viewers.

    Given user opens the address "/func/Compute2.ObjectCooling2"
    And user clicks on RUN button
    And "Line chart" viewer should be visible
    When user enters "20" into "Heat transfer coefficient" input
    # caption: The results say they are stale
    Then function results should contain text "Inputs changed"
    When user clicks on RUN button
    Then function results should not contain text "Inputs changed"
    And "Line chart" viewer should be visible

  Scenario: Share a link with the inputs
    The copy icon puts the inputs into the address.

    Given user opens the address "/func/Compute2.ObjectCooling2"
    When user enters "0.12" into "Surface area" input
    And user clicks on copy icon
    # caption: The link carries the inputs
    Then the clipboard should contain text "area=0.12"

  Scenario: Export a run
    A run goes to Excel or PowerPoint, its viewers included.

    Given user opens the address "/func/Compute2.ObjectCooling2"
    And user clicks on RUN button
    And "Line chart" viewer should be visible
    When user clicks on "Export the results" icon
    Then "Excel" menu item should be visible
    And "PowerPoint" menu item should be visible
    When the video says "Excel and PowerPoint, viewers included"
    And user hovers over "PowerPoint" menu item
    And user presses Escape

  Scenario: Fit a coefficient to a target
    Fitting finds the heat transfer coefficient that cools the object in an hour, and puts it into the form.

    Given user opens the address "/func/Compute2.ObjectCooling2"
    And user clicks on RUN button
    And "Line chart" viewer should be visible
    When user clicks on "Fitting" icon
    Then "Apply fitted parameters to the model" icon should be visible
    When user switches on "Heat transfer coefficient" in the analysis form
    And user enters "10" into "Heat transfer coefficient (min)" input
    And user enters "50" into "Heat transfer coefficient (max)" input
    And user switches off "Temp. vs time" in the analysis form
    And user switches on "Time to cool" in the analysis form
    And user enters "3600" into "Time to cool" input
    When the video says "Cool it down in one hour"
    And user clicks on "Run" icon
    Then status bar should contain text "Rows: 1"
    When user clicks on "Apply fitted parameters to the model" icon
    Then "Select fitted parameters" dialog should be visible
    When user clicks on OK button in "Select fitted parameters" dialog
    # caption: The fitted coefficient is in the form
    Then "Heat transfer coefficient" input should have a value between 44 and 44.6

  Scenario: Every run is kept in Run history
    Name a run, star it, and bring an older one back.

    Given user opens the address "/func/Compute2.ObjectCooling2"
    And user clicks on "Run history" icon
    And run history should be visible
    When user clicks on RUN button
    And user enters "20" into "Heat transfer coefficient" input
    And user clicks on RUN button
    # caption: What changed against the run before
    Then changes of first run should contain text "heatTransferCoeff 8.300 → 20"
    When user clicks on "Edit title, description and tags" icon in run header
    And user enters "Fast cooling" into "Title" input
    And user clicks on OK button
    Then run header should contain text "Fast cooling"
    When user clicks on run header star
    And user clicks on second run
    # caption: The older run is back
    Then "Heat transfer coefficient" input should have value "8.3"

  Scenario: Save a run to a space
    A run goes into a space next to the team's projects.

    Given the "PowerPack" package is installed
    And a space named "BDDDemoSpace{time}" is on the server
    And user opens the address "/func/Compute2.ObjectCooling2"
    And user clicks on RUN button
    And "Line chart" viewer should be visible
    And user clicks on "Edit title, description and tags" icon in run header
    And user enters "Fast cooling" into "Title" input
    And user clicks on OK button
    And run header should contain text "Fast cooling"
    When user clicks on "Save this run to a space" icon in run header
    Then "Save Fast cooling to a space" dialog should be visible
    When user clicks on space input in "Save Fast cooling to a space" dialog
    And user clicks on "BDDDemoSpace{time}" tree node
    And user clicks on OK button in "Save Fast cooling to a space" dialog
    Then the "Save Fast cooling to a space" dialog should close
    And the space "BDDDemoSpace{time}" on the server should link a run of the "ObjectCooling2" function
