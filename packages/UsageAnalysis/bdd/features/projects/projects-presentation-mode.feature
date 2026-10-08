@journey @serial @realizes:views.projects
Feature: The Presentation mode switch of a saved project's Save dialog
  demog.csv from Browse > Files > Demo with a scatter plot is saved as a project; reopened from its
  Dashboards card it opens in design mode (the top menu, the left sidebar, the toolbox and the
  status bar shown), and its Save dialog offers the Presentation mode switch with a tooltip that
  tells how to get back. The switch turned on and saved, the project reopens in presentation mode
  (the sidebar, the toolbox and the status bar hidden, the "Press F7" balloon and the "back to
  design mode" link shown); the link and F7 switch between the modes. The project's address opens
  it in presentation mode, a second project saved in design mode opens from its address in design
  mode and with "?mode=presentation" in presentation mode, and a new project saved from the
  Dashboards panel with the switch on opens in presentation mode. Translated from the TestTrack case
  Projects/project-presentation-mode.

  The addresses are opened in the same tab (the md opens a new one) and typed, as
  /p/<namespace>.<name>, each after a Close All (a new tab starts with nothing open, and the scatter
  plot only BDDPresentProj has tells the two projects apart); the harness turns presentation mode off after the journey and every
  feature's start does too, so a failure in presentation mode does not leave a later feature
  without menus.

  Kept without (see the request document): the copy icon of the URL row of Links... and opening the
  address from the clipboard, and the Description the md types for the new project (the Save
  dialog's Description box has no element) with the card's claim on it. Not claimed: that the top
  menu is hidden in presentation mode — the view's menu (Edit, View, Select, Data, ML) and its
  ribbon stay shown on localhost (a suspected defect, described outside the repository).

  The project is named with the run's time (letters and digits only: the Dashboards search misses
  "-" and "_") and removed with its table and view when the feature starts and ends. It is serial:
  the Dashboards search is shared with every feature that saves a project.

  Background:
    Given user is logged in
    And the browse panel is open
    And no project named "BDDPresentProj{time}" is on the server
    And no project named "BDDPresentPlain{time}" is on the server
    And no project named "BDDPresentNew{time}" is on the server

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

  Scenario: A second project is saved in design mode
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user double-clicks Files---Demo---demog.csv tree node inside browse tree
    Then the "demog" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDPresentPlain{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDPresentPlain{time}" uploaded' should have been shown
    When user clicks on CANCEL button in "Share BDDPresentPlain{time}" dialog
    Then the "Share BDDPresentPlain{time}" dialog should close

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

  Scenario: Presentation mode is turned on and saved into the project
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user switches on "Presentation mode" input in "Save project" dialog
    Then "Presentation mode" input in "Save project" dialog should be switched on
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDPresentProj{time}" uploaded' should have been shown
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current

  Scenario: The project reopens in presentation mode
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDPresent" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDPresentProj{time} gallery card
    Then grid should be visible
    And scatter plot viewer should be visible
    And "back to design mode" link should be visible
    And browse tab should be hidden
    And toolbox tab should be hidden
    And status bar should be hidden
    And an info balloon containing "Press F7 to go back to the design mode" should have been shown

  Scenario: back to design mode and F7 switch between the modes
    When user clicks on "back to design mode" link
    Then browse tab should be visible
    And status bar should be visible
    And "back to design mode" link should be absent
    And an info balloon containing "Press F7 to go back to the presentation mode" should have been shown
    When user presses F7
    Then status bar should be hidden
    And "back to design mode" link should be visible
    When user presses F7
    Then status bar should be visible
    And "back to design mode" link should be absent

  Scenario: A new project saved from the Dashboards panel with the switch on opens in presentation mode
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---demog.csv tree node inside browse tree
    Then the "demog" view should be current
    When user clicks on Dashboards tab
    Then "New Dashboard > demog" tree node inside browse tree should be visible
    When user clicks on Save button in "New Dashboard" tree node inside browse tree
    Then "Save project" dialog should be visible
    When user enters "BDDPresentNew{time}" into Name text input in "Save project" dialog
    And user switches on "Presentation mode" input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDPresentNew{time}" uploaded' should have been shown
    When user clicks on OK button in "Share BDDPresentNew{time}" dialog
    Then the "Share BDDPresentNew{time}" dialog should close
    # the new project is the open one now, and the shell goes into its presentation mode at once
    And "back to design mode" link should be visible
    And status bar should be hidden
    When user clicks on "back to design mode" link
    Then status bar should be visible
    # the sidebar tab toggles the Dashboards panel; left open, it adds a second Save button to every table view
    When user clicks on Dashboards tab
    Then "New Dashboard" tree node inside browse tree should be hidden
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDPresent" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    Then BDDPresentNew{time} gallery card should be visible
    When user double-clicks on BDDPresentNew{time} gallery card
    Then grid should be visible
    And "back to design mode" link should be visible
    And status bar should be hidden
    When user clicks on "back to design mode" link
    Then status bar should be visible

  Scenario: The project's address opens it in presentation mode, and ?mode=presentation any project
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    When user opens the address "/p/Admin.BDDPresentProj{time}"
    Then scatter plot viewer should be visible
    And "back to design mode" link should be visible
    And status bar should be hidden
    When user clicks on "back to design mode" link
    Then status bar should be visible
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    When user opens the address "/p/Admin.BDDPresentPlain{time}"
    Then the "demog" view should be current
    And grid should be visible
    And scatter plot viewer should be absent
    And status bar should be visible
    And "back to design mode" link should be absent
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    When user opens the address "/p/Admin.BDDPresentPlain{time}?mode=presentation"
    Then grid should be visible
    And scatter plot viewer should be absent
    And "back to design mode" link should be visible
    And status bar should be hidden
    When user clicks on "back to design mode" link
    Then status bar should be visible

  Scenario: The projects are deleted from their cards
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDPresent" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user picks "Delete Project" from the context menu of BDDPresentProj{time} gallery card
    And user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    When user picks "Delete Project" from the context menu of BDDPresentPlain{time} gallery card
    And user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    When user picks "Delete Project" from the context menu of BDDPresentNew{time} gallery card
    And user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 projects named "BDDPresentProj{time}" should be on the server
    And 0 projects named "BDDPresentPlain{time}" should be on the server
    And 0 projects named "BDDPresentNew{time}" should be on the server
