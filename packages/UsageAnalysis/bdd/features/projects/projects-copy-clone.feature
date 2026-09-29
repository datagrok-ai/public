@journey @serial @realizes:views.projects @realizes:sharing.share-dialog
Feature: A project saved as the original, as a copy with link and as a copy with clone
  demog.csv with a bar chart is saved as a project, shared with the second account from its
  Dashboards card and saved again with a scatter plot ("Save original project"). A line chart is
  added and the project saved as a copy whose table is linked to the original's; the original still
  has no line chart (GROK-19750). A histogram is added and the project saved as a copy with its own
  (cloned) table. The copies are shared too. The link copy's Content pane marks demog as linked;
  deleting the original leaves the demo file and the clone copy working, while the link copy loses
  the table. Translated from the TestTrack case Projects/projects-copy-clone.

  Parked (see the request document): the card's thumbnail (no reading of a card's picture); the
  personal view customizations (sort, hidden column, filter, "Save personal view customizations",
  the Custom views pane after the save, Reset, and the Save dialog starting in that mode), which
  also need their user-data entry swept at the start and the end; opening each project from the URL
  of Links... in a new tab; the second account's opens; and what opening the link copy says once
  the original is gone ("a dialog or a balloon", which no single claim reads — on localhost 1.28.0
  the double click left the Projects view current and showed nothing, a suspected defect).

  Every save goes through the ribbon's Save dialog and every reopen through the Dashboards gallery.
  Names are letters and digits only (the Dashboards search misses "-" and "_") and carry the run's
  time. The three projects are deleted with their tables, views and grants when the feature starts
  and ends. It is serial: the Dashboards search and the uploads are shared with every feature that
  saves a project.

  Background:
    Given user is logged in
    And simple mode is off
    And the browse panel is open
    And no project named "BDDCopyClone{time}" is on the server
    And no project named "BDDCopyCloneLink{time}" is on the server
    And no project named "BDDCopyCloneClone{time}" is on the server

  Scenario: The original is saved from demog with a bar chart
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---demog.csv tree node inside browse tree
    Then the "demog" view should be current
    Given the toolbox pane is shown
    When user clicks on "bar chart" icon in toolbox
    Then bar chart viewer should be visible
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And "Creation script" button in "demog" project table in "Save project" dialog should be visible
    And Data sync switch in "demog" project table in "Save project" dialog should be checked
    When user enters "BDDCopyClone{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDCopyClone{time}" uploaded' should have been shown
    And "Share BDDCopyClone{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDCopyClone{time}" dialog
    Then the "Share BDDCopyClone{time}" dialog should close

  Scenario: The card's context panel names the project, its author and its table, with no custom views
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDCopyClone{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user clicks on BDDCopyClone{time} gallery card
    Then the context panel should show "BDDCopyClone{time}"
    And Details pane in context panel should contain text "Created by"
    And Details pane in context panel should contain text "Admin"
    Given Content pane in context panel is expanded
    Then Content pane in context panel should contain text "demog"
    Given "Custom views" pane in context panel is expanded
    Then "Custom views" pane in context panel should contain text "No personal view customizations"

  Scenario: The original is shared with the second account
    When user picks "Share..." from the context menu of BDDCopyClone{time} gallery card
    Then "Share BDDCopyClone{time}" dialog should be visible
    And share access selector should contain text "View and use"
    When user picks the sharing user in "User, group, or email" input in "Share BDDCopyClone{time}" dialog
    And user clicks on OK button in "Share BDDCopyClone{time}" dialog
    Then the "Share BDDCopyClone{time}" dialog should close
    When user clicks on BDDCopyClone{time} gallery card
    Then the context panel should show "BDDCopyClone{time}"
    And the sharing pane should list the sharing user

  Scenario: The original is saved again with a scatter plot
    When user double-clicks on BDDCopyClone{time} gallery card
    Then the "demog" view should be current
    And bar chart viewer should be visible
    Given the toolbox pane is shown
    When user clicks on "scatter plot" icon in toolbox
    Then scatter plot viewer should be visible
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDCopyClone{time}" uploaded' should have been shown

  Scenario: A copy with link is saved with a line chart
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDCopyClone{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDCopyClone{time} gallery card
    Then the "demog" view should be current
    Given the toolbox pane is shown
    When user clicks on "line chart" icon in toolbox
    Then line chart viewer should be visible
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user selects "Save a copy" in radio input in "Save project" dialog
    Then Name text input in "Save project" dialog should have value "Copy of BDDCopyClone{time}"
    When user enters "BDDCopyCloneLink{time}" into Name text input in "Save project" dialog
    And user selects "Link" in choice input in "demog" project table in "Save project" dialog
    Then choice input in "demog" project table in "Save project" dialog should have value "Link"
    And "Local data changes will not be saved" text in "demog" project table in "Save project" dialog should be visible
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDCopyCloneLink{time}" uploaded' should have been shown

  Scenario: The copy with link opens with all four viewers
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDCopyCloneLink{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDCopyCloneLink{time} gallery card
    Then the "demog" view should be current
    And grid should be visible
    And line chart viewer should be visible
    And bar chart viewer should be visible
    And scatter plot viewer should be visible
    And the "rows" reading of grid should be at least 1

  Scenario: The original still has no line chart (GROK-19750)
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDCopyClone{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDCopyClone{time} gallery card
    Then the "demog" view should be current
    And grid should be visible
    And bar chart viewer should be visible
    And scatter plot viewer should be visible
    And line chart viewer should be absent

  Scenario: A copy with clone is saved with a histogram and opens with its own table
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDCopyClone{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDCopyClone{time} gallery card
    Then the "demog" view should be current
    Given the toolbox pane is shown
    When user clicks on "histogram" icon in toolbox
    Then histogram viewer should be visible
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user selects "Save a copy" in radio input in "Save project" dialog
    And user enters "BDDCopyCloneClone{time}" into Name text input in "Save project" dialog
    Then choice input in "demog" project table in "Save project" dialog should have value "Clone"
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDCopyCloneClone{time}" uploaded' should have been shown
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDCopyCloneClone{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDCopyCloneClone{time} gallery card
    Then the "demog" view should be current
    And histogram viewer should be visible
    And bar chart viewer should be visible
    And scatter plot viewer should be visible
    And the table should have 5850 rows

  Scenario: The copies are shared with the second account
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDCopyCloneLink{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user picks "Share..." from the context menu of BDDCopyCloneLink{time} gallery card
    Then "Share BDDCopyCloneLink{time}" dialog should be visible
    When user picks the sharing user in "User, group, or email" input in "Share BDDCopyCloneLink{time}" dialog
    And user clicks on OK button in "Share BDDCopyCloneLink{time}" dialog
    Then the "Share BDDCopyCloneLink{time}" dialog should close
    When user enters "BDDCopyCloneClone{time}" into gallery search
    And user picks "Share..." from the context menu of BDDCopyCloneClone{time} gallery card
    Then "Share BDDCopyCloneClone{time}" dialog should be visible
    When user picks the sharing user in "User, group, or email" input in "Share BDDCopyCloneClone{time}" dialog
    And user clicks on OK button in "Share BDDCopyCloneClone{time}" dialog
    Then the "Share BDDCopyCloneClone{time}" dialog should close
    When user clicks on BDDCopyCloneClone{time} gallery card
    Then the context panel should show "BDDCopyCloneClone{time}"
    And the sharing pane should list the sharing user
    When user enters "BDDCopyCloneLink{time}" into gallery search
    And user clicks on BDDCopyCloneLink{time} gallery card
    Then the context panel should show "BDDCopyCloneLink{time}"
    And the sharing pane should list the sharing user

  Scenario: The copy with link refers to the original's table
    When user enters "BDDCopyClone" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user clicks on BDDCopyCloneLink{time} gallery card
    Then the context panel should show "BDDCopyCloneLink{time}"
    Given Content pane in context panel is expanded
    Then Content pane in context panel should contain text "demog"
    When user hovers over link icon in demog tree node in Content pane in context panel
    Then tooltip should contain text "This entity is not included to this project, but linked."

  Scenario: The original is deleted; the copies stay, and so does the demo file
    When user picks "Delete Project" from the context menu of BDDCopyClone{time} gallery card
    Then "Are you sure?" dialog should contain text 'Delete project "BDDCopyClone{time}"?'
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 projects named "BDDCopyClone{time}" should be on the server
    When user clicks on "Refresh" icon inside gallery toolbar
    Then BDDCopyClone{time} gallery card should be absent
    And BDDCopyCloneLink{time} gallery card should be visible
    And BDDCopyCloneClone{time} gallery card should be visible
    Given Files tree node inside browse tree is expanded
    When user clicks on "Files > Demo" tree node inside browse tree
    Then the "Demo" view should be current
    And demog.csv link in gallery should be visible

  Scenario: The copy with clone still opens with its data
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDCopyClone" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDCopyCloneClone{time} gallery card
    Then the "demog" view should be current
    And the table should have 5850 rows
    And "Data loading error" dialog should be absent
    And no error or warning balloon should have been shown

  Scenario: The copy with link lost the original's table
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDCopyClone" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user clicks on BDDCopyCloneLink{time} gallery card
    Then the context panel should show "BDDCopyCloneLink{time}"
    Given Content pane in context panel is expanded
    # the pane lists the Demo files connection the table was read from; once it shows, demog is not there
    Then Content pane in context panel should contain text "Demo"
    And Content pane in context panel should not contain text "demog"
