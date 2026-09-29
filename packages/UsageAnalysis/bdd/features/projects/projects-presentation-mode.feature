@journey @serial @realizes:views.projects
Feature: The Presentation mode switch of a saved project's Save dialog
  demog.csv from Browse > Files > Demo with a scatter plot is saved as a project; reopened from its
  Dashboards card it opens in design mode (the top menu, the left sidebar, the toolbox and the
  status bar shown), and its Save dialog offers the Presentation mode switch with a tooltip that
  tells how to get back. Translated from the TestTrack case Projects/project-presentation-mode,
  steps 1 and 2 up to the switch.

  Parked (see the request document): everything that puts the shell into presentation mode — the
  switch turned on and saved, the reopen in presentation mode, "back to design mode" and F7, the
  new project saved from the Dashboards panel with a description and the switch on — until the
  library can put presentation mode back at the feature end (a failure in presentation mode would
  leave every later feature of the worker without menus); the project's link and
  "?mode=presentation" in a new tab (no step opens an address); and the second project the md
  saves for the "?mode=presentation" case.

  The project is named with the run's time (letters and digits only: the Dashboards search misses
  "-" and "_") and removed with its table and view when the feature starts and ends. It is serial:
  the Dashboards search is shared with every feature that saves a project.

  Background:
    Given user is logged in
    And the browse panel is open
    And no project named "BDDPresentProj{time}" is on the server

  Scenario: demog with a scatter plot is saved as a project
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---demog.csv tree node inside browse tree
    Then the "demog" view should be current
    Given the toolbox pane is shown
    When user clicks on "scatter plot" icon in toolbox
    Then scatter plot viewer should be visible
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDPresentProj{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDPresentProj{time}" uploaded' should have been shown
    And "Share BDDPresentProj{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDPresentProj{time}" dialog
    Then the "Share BDDPresentProj{time}" dialog should close

  Scenario: The project opens in design mode
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDPresent" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDPresentProj{time} gallery card
    Then the "demog" view should be current
    And scatter plot viewer should be visible
    And "Data" menu item should be visible
    And browse tab should be visible
    And status bar should be visible
    And toolbox tab should be visible

  Scenario: The Save dialog offers Presentation mode, and its tooltip tells how to get back
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And "Presentation mode" input in "Save project" dialog should be visible
    When user hovers over "Presentation mode" input in "Save project" dialog
    Then tooltip should contain text "visualization"
    And tooltip should contain text "Design mode"
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close

  Scenario: The project is deleted from its card
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDPresentProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user picks "Delete Project" from the context menu of BDDPresentProj{time} gallery card
    And user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 projects named "BDDPresentProj{time}" should be on the server
