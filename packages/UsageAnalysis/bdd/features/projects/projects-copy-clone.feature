@journey @serial @realizes:views.projects @realizes:sharing.share-dialog
Feature: A project saved as itself, as a copy with link, as a copy with clone and as personal views
  demog.csv opened from Browse > Files > Demo gets a bar chart and is saved through the ribbon's Save
  dialog; its card in the Dashboards gallery is previewed and shared with the second account. The
  original is then saved again with a scatter plot, saved as a copy with link (with a line chart the
  original must not get, GROK-19750), saved as a copy with clone (with a histogram), and its view is
  customized and saved as personal view customizations. Every result is reopened from its gallery
  card; the copies are shared too, and the second account opens all three. Last, each project's
  link is copied from Context Panel > Details > Links... and opened in a new tab of the same
  browser, where the project it names must open. Translated from the TestTrack cases
  Projects/projects-copy-clone, Projects/project-url and the automatable half of
  Projects/projects-copy-clone-ui.

  From projects-copy-clone-ui: the card's thumbnail is claimed as the project's own picture (the
  server's entity picture, decoded by the browser, not the gallery placeholder); the personal view
  customizations are claimed to bring back a range filter, the sort, the hidden column and a
  viewer's dock position, each read from the view's state. Not translated from it: tile sizes,
  truncation, overlap and jitter while scrolling (a judgement of the picture, which the suite does
  not make), and "the rest of the workspace has no drift" (no criterion to read). The case also
  names the personal customizations as a copy ("name the copy test_copy_clone_pvc"); the product
  does not make one: that mode disables the name and saves the views into the account's settings
  under the original project, which is what is claimed here.

  projects-copy-clone step 8 says the personal-customizations dialog "has no name field". It has
  one, disabled (probed 2026-09-25): the claim is that the Name field is disabled. The Order or Hide
  Columns dialog closes with CLOSE, not OK. Reopening a project with personal customizations shows
  the warning balloon "Project has personal view customizations", which is claimed; the Context
  Panel's Custom views pane says so too. The dialog keeps the customizations in the page's copy of
  the account's settings, which the page sends to the server on a 10 s timer; the save is claimed on
  the server, since the second account's sign-in reloads the page and the new tab is another session. For the second account the same original opens unsorted:
  the customizations are the author's alone.

  The second account is the library's sharing user; it is switched to by a session of its own, not
  by the platform's Logout (which would end the session every worker shares). The link of a project
  reads <server>/p/<owner namespace>.<name> — on the stand the namespace is "Admin", with a capital,
  where the md says "your login". The preview's "author" is not claimed (the account differs per
  stand); the name, and the Content pane listing demog once it is opened, are. The Save dialog's
  preview logs "Unable to find element in cloned iframe" (GROK-18606, known noise), so a save is
  claimed to log nothing else; every reopen is claimed to log nothing at all. The new tab's balloons
  are recorded from its first script; its console is not, so "no errors should have been logged"
  after a tab closes reads the main page only.

  md's "a copy with clone has its own copy of the data" is claimed on the server: the clone copy
  holds a TableInfo of its own, while the link copy holds the original's. Every Close All is claimed
  to leave no table open, so a reopen cannot be answered by a table left in the workspace.

  Everything this feature saves goes when it ends, and anything an earlier run left goes when it
  starts: the three projects with their tables and views (the grants go with them), and the
  personal view customizations the account keeps for projects that are gone. It is @serial: the
  saves upload a table and the second account reads it back.

  Background:
    Given user is logged in
    And the browse panel is open
    And no project named "BDDCopy{time}" is on the server
    And no project named "BDDCopyLink{time}" is on the server
    And no project named "BDDCopyClone{time}" is on the server
    And no personal view customizations of deleted projects are kept

  Scenario: The original is saved from a file with a bar chart
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---demog.csv tree node inside browse tree
    Then the table should have 5850 rows
    When user clicks on "bar chart" icon in toolbox
    Then the open tableview should have 1 bar chart viewer
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "demog" project table in "Save project" dialog should be switched on
    When user enters "BDDCopy{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDCopy{time}" uploaded' should have been shown
    And 1 project named "BDDCopy{time}" should be on the server
    And "Share BDDCopy{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDCopy{time}" dialog
    Then the "Share BDDCopy{time}" dialog should close
    When user picks "Close All" from the context menu of left sidebar
    Then the "Home" view should be current
    And no table should be left in the workspace
    And no errors but the project preview's should have been logged

  Scenario: The card previews the project
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    Then the "Projects" view should be current
    When user enters "BDDCopy{time}" into gallery search
    Then BDDCopy{time} gallery card should be visible
    And BDDCopy{time} gallery card should show the project's own picture
    When user clicks on BDDCopy{time} gallery card
    Then the context panel should show "BDDCopy{time}"
    And Details section in context panel should be visible
    And Content section in context panel is collapsed
    And context panel should not contain text "demog"
    When user clicks on Content section in context panel
    Then context panel should contain text "demog"
    And no errors should have been logged

  Scenario: The original is shared with the second account from its card
    When user picks "Share..." from the context menu of BDDCopy{time} gallery card
    Then "Share BDDCopy{time}" dialog should be visible
    And share access selector should contain text "View and use"
    When user picks the sharing user in "User, group, or email" input in "Share BDDCopy{time}" dialog
    And user clicks on OK button in "Share BDDCopy{time}" dialog
    Then the "Share BDDCopy{time}" dialog should close
    When user clicks on BDDCopy{time} gallery card
    Then the context panel should show "BDDCopy{time}"
    And the sharing pane should list the sharing user
    And no errors should have been logged

  Scenario: The original is saved again with a scatter plot
    When user double-clicks on BDDCopy{time} gallery card
    Then the table should have 5850 rows
    And the open tableview should have 1 bar chart viewer
    When user clicks on "scatter plot" icon in toolbox
    Then the open tableview should have 1 scatter plot viewer
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    When user picks "Close All" from the context menu of left sidebar
    Then the "Home" view should be current
    And no table should be left in the workspace
    And no errors but the project preview's should have been logged

  Scenario: A copy with link is saved with a line chart
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDCopy{time}" into gallery search
    And user double-clicks on BDDCopy{time} gallery card
    Then the table should have 5850 rows
    And the open tableview should have 1 scatter plot viewer
    When user clicks on "line chart" icon in toolbox
    Then the open tableview should have 1 line chart viewer
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user selects "Save a copy" in radio input in "Save project" dialog
    Then Name text input in "Save project" dialog should have value "Copy of BDDCopy{time}"
    And choice input in "demog" project table in "Save project" dialog should have value "Clone"
    And "Local data changes will not be saved" text in "demog" project table in "Save project" dialog should be hidden
    When user enters "BDDCopyLink{time}" into Name text input in "Save project" dialog
    And user selects "Link" in choice input in "demog" project table in "Save project" dialog
    Then "Local data changes will not be saved" text in "demog" project table in "Save project" dialog should be visible
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And 1 project named "BDDCopyLink{time}" should be on the server
    And the "demog" table of the "BDDCopyLink{time}" project should be the one of the "BDDCopy{time}" project
    When user picks "Close All" from the context menu of left sidebar
    Then the "Home" view should be current
    And no table should be left in the workspace
    And no errors but the project preview's should have been logged

  Scenario: The copy with link opens with all four viewers
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDCopyLink{time}" into gallery search
    And user double-clicks on BDDCopyLink{time} gallery card
    Then the table should have 5850 rows
    And the open tableview should have 1 line chart viewer
    And the open tableview should have 1 bar chart viewer
    And the open tableview should have 1 scatter plot viewer
    When user picks "Close All" from the context menu of left sidebar
    Then the "Home" view should be current
    And no table should be left in the workspace
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The original has no line chart after the copy with link (GROK-19750)
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDCopy{time}" into gallery search
    And user double-clicks on BDDCopy{time} gallery card
    Then the table should have 5850 rows
    And the open tableview should have 1 bar chart viewer
    And the open tableview should have 1 scatter plot viewer
    And the open tableview should have 0 line chart viewers
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A copy with clone is saved with a histogram and a table of its own
    When user clicks on "histogram" icon in toolbox
    Then the open tableview should have 1 histogram viewer
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user selects "Save a copy" in radio input in "Save project" dialog
    And user enters "BDDCopyClone{time}" into Name text input in "Save project" dialog
    Then choice input in "demog" project table in "Save project" dialog should have value "Clone"
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And 1 project named "BDDCopyClone{time}" should be on the server
    And the "demog" table of the "BDDCopyClone{time}" project should be its own, not the one of the "BDDCopy{time}" project
    And no errors but the project preview's should have been logged
    When user picks "Close All" from the context menu of left sidebar
    Then no table should be left in the workspace
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDCopyClone{time}" into gallery search
    And user double-clicks on BDDCopyClone{time} gallery card
    Then the table should have 5850 rows
    And the open tableview should have 1 histogram viewer
    And the open tableview should have 1 bar chart viewer
    And the open tableview should have 1 scatter plot viewer
    And the open tableview should have 0 line chart viewers
    When user picks "Close All" from the context menu of left sidebar
    Then the "Home" view should be current
    And no table should be left in the workspace
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Personal view customizations bring back the filter, the sort, the hidden column and the dock
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDCopy{time}" into gallery search
    And user double-clicks on BDDCopy{time} gallery card
    Then the table should have 5850 rows
    When user adds a range filter on "AGE" from 50 to 100
    Then 2345 rows should pass the filter
    And the "sort column" reading of grid should be ""
    When user picks "Sort > Ascending" from the context menu of the "header AGE" area of grid
    Then the "sort column" reading of grid should be "AGE"
    And the "sort direction" reading of grid should be "ascending"
    When user picks "Order or Hide Columns..." from the context menu of the "header AGE" area of grid
    Then Order or Hide Columns dialog should be visible
    And the grid should show the "DIS_POP" column
    When user toggles the "DIS_POP" column in the column list of Order or Hide Columns dialog
    And user clicks on CLOSE button in Order or Hide Columns dialog
    Then the grid should hide the "DIS_POP" column
    And bar chart viewer should not be docked along the left edge of the view
    When user docks bar chart viewer to the left edge of the view
    Then bar chart viewer should be docked along the left edge of the view
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user selects "Save personal view customizations" in radio input in "Save project" dialog
    Then Name text input in "Save project" dialog should be disabled
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And the personal view customizations of the "BDDCopy{time}" project should be on the server
    And no errors but the project preview's should have been logged
    When user picks "Close All" from the context menu of left sidebar
    Then no table should be left in the workspace
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDCopy{time}" into gallery search
    And user clicks on BDDCopy{time} gallery card
    Then the context panel should show "BDDCopy{time}"
    When user clicks on "Custom views" section in context panel
    Then context panel should contain text "Project has personal view customizations"
    When user double-clicks on BDDCopy{time} gallery card
    Then the table should have 5850 rows
    And a warning balloon containing "Project has personal view customizations" should have been shown
    And the filter panel should have 1 filter
    And 2345 rows should pass the filter
    And the "sort column" reading of grid should be "AGE"
    And the "sort direction" reading of grid should be "ascending"
    And the grid should hide the "DIS_POP" column
    And bar chart viewer should be docked along the left edge of the view
    And the open tableview should have 0 line chart viewers
    And the open tableview should have 0 histogram viewers
    When user picks "Close All" from the context menu of left sidebar
    Then the "Home" view should be current
    And no table should be left in the workspace
    And no errors should have been logged

  Scenario: The copies are shared with the second account from their cards
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDCopyLink{time}" into gallery search
    And user picks "Share..." from the context menu of BDDCopyLink{time} gallery card
    Then "Share BDDCopyLink{time}" dialog should be visible
    When user picks the sharing user in "User, group, or email" input in "Share BDDCopyLink{time}" dialog
    And user clicks on OK button in "Share BDDCopyLink{time}" dialog
    Then the "Share BDDCopyLink{time}" dialog should close
    When user clicks on BDDCopyLink{time} gallery card
    Then the context panel should show "BDDCopyLink{time}"
    And the sharing pane should list the sharing user
    When user enters "BDDCopyClone{time}" into gallery search
    And user picks "Share..." from the context menu of BDDCopyClone{time} gallery card
    Then "Share BDDCopyClone{time}" dialog should be visible
    When user picks the sharing user in "User, group, or email" input in "Share BDDCopyClone{time}" dialog
    And user clicks on OK button in "Share BDDCopyClone{time}" dialog
    Then the "Share BDDCopyClone{time}" dialog should close
    When user clicks on BDDCopyClone{time} gallery card
    Then the context panel should show "BDDCopyClone{time}"
    And the sharing pane should list the sharing user
    And no errors should have been logged

  Scenario: The second account opens the original and both copies
    Given user signs in as the sharing user
    And the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDCopy{time}" into gallery search
    And user double-clicks on BDDCopy{time} gallery card
    Then the table should have 5850 rows
    And the open tableview should have 1 bar chart viewer
    And the open tableview should have 1 scatter plot viewer
    And the open tableview should have 0 line chart viewers
    And the "sort column" reading of grid should be ""
    And the grid should show the "DIS_POP" column
    When user picks "Close All" from the context menu of left sidebar
    Then no table should be left in the workspace
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDCopyLink{time}" into gallery search
    And user double-clicks on BDDCopyLink{time} gallery card
    Then the table should have 5850 rows
    And the open tableview should have 1 line chart viewer
    And the open tableview should have 1 bar chart viewer
    And the open tableview should have 1 scatter plot viewer
    When user picks "Close All" from the context menu of left sidebar
    Then no table should be left in the workspace
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDCopyClone{time}" into gallery search
    And user double-clicks on BDDCopyClone{time} gallery card
    Then the table should have 5850 rows
    And the open tableview should have 1 histogram viewer
    And the open tableview should have 1 bar chart viewer
    And the open tableview should have 1 scatter plot viewer
    When user picks "Close All" from the context menu of left sidebar
    Then the "Home" view should be current
    And no table should be left in the workspace
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The link of the original opens it in a new tab, with its personal views
    Given user signs in as themselves again
    And the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDCopy{time}" into gallery search
    And user clicks on BDDCopy{time} gallery card
    Then the context panel should show "BDDCopy{time}"
    When user clicks on Links... link in context panel
    Then "Links to BDDCopy{time}" dialog should be visible
    When user clicks on copy icon in URL input in "Links to BDDCopy{time}" dialog
    Then the clipboard should hold the URL of the "BDDCopy{time}" project
    When user presses Escape
    And user opens the copied link in a new tab
    Then the new tab should show the "BDDCopy{time}" project with the viewers "Grid, Bar chart, Scatter plot"
    And the new tab should have shown a warning balloon containing "Project has personal view customizations"
    And the table in the new tab should have 2345 rows passing the filter
    And the grid in the new tab should be sorted by "AGE ascending"
    And the grid in the new tab should hide the "DIS_POP" column
    And the "Bar chart" viewer in the new tab should be docked along the left edge of the view
    And the new tab should have shown no error balloon
    When user closes the new tab
    Then no errors should have been logged

  Scenario: The links of the copies open the copies, not the original
    When user enters "BDDCopyLink{time}" into gallery search
    And user clicks on BDDCopyLink{time} gallery card
    Then the context panel should show "BDDCopyLink{time}"
    When user clicks on Links... link in context panel
    Then "Links to BDDCopyLink{time}" dialog should be visible
    When user clicks on copy icon in URL input in "Links to BDDCopyLink{time}" dialog
    Then the clipboard should hold the URL of the "BDDCopyLink{time}" project
    When user presses Escape
    And user opens the copied link in a new tab
    Then the new tab should show the "BDDCopyLink{time}" project with the viewers "Grid, Bar chart, Scatter plot, Line chart"
    And the new tab should have shown no error balloon
    When user closes the new tab
    And user enters "BDDCopyClone{time}" into gallery search
    And user clicks on BDDCopyClone{time} gallery card
    Then the context panel should show "BDDCopyClone{time}"
    When user clicks on Links... link in context panel
    Then "Links to BDDCopyClone{time}" dialog should be visible
    When user clicks on copy icon in URL input in "Links to BDDCopyClone{time}" dialog
    Then the clipboard should hold the URL of the "BDDCopyClone{time}" project
    When user presses Escape
    And user opens the copied link in a new tab
    Then the new tab should show the "BDDCopyClone{time}" project with the viewers "Grid, Bar chart, Scatter plot, Histogram"
    And the new tab should have shown no error balloon
    When user closes the new tab
    Then no errors should have been logged
