@journey @serial @realizes:views.projects
Feature: A project built from the user's own script: shared, renamed, broken
  A JavaScript script is written in a new editor (Scripts > New > JavaScript Script...) and saved;
  Run... from the Scripts gallery opens its demog output, which is saved as a project with Data
  sync on and the script as its creation script. The project is shared with the second account
  through the Dashboards tile's Share..., which gives it the script too (the script's own Share dialog
  is read, as the md says of today's product), and the second account gets the data (GROK-19403). The
  script is renamed in its editor and the project still opens; the project is renamed and still opens.
  Then the script is broken (it throws): the second account, with View and use access, gets the Data
  loading error dialog asking it to turn to the owner and no EDIT SCRIPT... (GROK-19728); the owner
  then gets the same dialog with OPEN ANYWAY, EDIT SCRIPT... and CLOSE PROJECT. Last, Delete Project
  and the script's Delete remove both. Translated from the TestTrack case
  Projects/projects-lifecycle-script.

  Order: the md has the owner open the broken project first and the second account after. Here the
  second account goes first, so the owner's open follows a switch of accounts, which reloads the page.
  A session that ran the script under its new name keeps that body for a while after the script is
  saved again (the script save reaches the page with a delay), and the owner's session ran it at the
  opens after the renames; after the reload the page takes the script from the server.

  The script runs in the page (grok.data.getDemoTable reads the demo file from the server), so no
  script container is involved. The script and the project are named with the run's time and removed,
  under both their names, at the start and at the end; the grant the Share dialog adds lives on the
  project and the script and goes with them. @serial: the Dashboards search and the Scripts view's
  search (an account setting, cleared at the end) are shared with the features that use them.

  Not translated, and why: Logout and signing in with the second user's credentials, and the reload
  of the browser tab — the platform's Logout ends every session of the account, which all workers of
  a run share, so the second account is entered through its own session, which reloads the page;
  Close All from the left sidebar's context menu is done through the shell (closing views is not the
  claim). The md changes the script's first line and inserts the throw before the df line; the editor
  has no line-replacing gesture, so the whole text is typed again with the new name, or with the throw
  in that place.

  Background:
    Given user is logged in
    And the browse panel is open
    And no project named "BDDLifeScriptProj{time}" is on the server
    And no project named "BDDLifeScriptProjRenamed{time}" is on the server
    And no script named "BDDLifeScript{time}" is on the server
    And no script named "BDDLifeScriptRenamed{time}" is on the server

  Scenario: A JavaScript script is written in a new editor and saved
    Given Platform tree node inside browse tree is expanded
    And Platform---Functions tree node inside browse tree is expanded
    When user clicks on Platform---Functions---Scripts tree node inside browse tree
    Then the "Scripts" view should be current
    When user clicks on New button
    And user picks "JavaScript Script..." from the open menu
    Then the "Template" view should be current
    And code editor should contain the text "Hello World"
    When user replaces the code of code editor with "//name: BDDLifeScript{time}"
    And user appends "//language: javascript" to code editor
    And user appends "//output: dataframe df" to code editor
    And user appends "df = await grok.data.getDemoTable('demog.csv');" to code editor
    And user saves the script
    Then an info balloon containing "Script saved." should have been shown
    And 1 script named "BDDLifeScript{time}" should be on the server
    And no errors should have been logged

  Scenario: Run... from the Scripts gallery opens the script's table
    When user closes the current view
    Then the "Scripts" view should be current
    When user clears gallery search
    And user types "BDDLifeScript{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user picks "Run..." from the context menu of "BDDLifeScript{time}" link in gallery
    Then the current view should be a TableView view
    And the "demog" view should be current
    And the table should have 5850 rows
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Saved with Data sync, the table's creation script calls the script
    When user clicks on Save button
    Then "Save project" dialog should be visible
    And there should be 1 visible saved table row
    And Data sync switch in "demog" project table in "Save project" dialog should be checked
    And "Save project" dialog should contain text "Some tables require this script for data sync."
    When user clicks on "Creation script" button in "demog" project table in "Save project" dialog
    Then creation script text in "demog" project table in "Save project" dialog should be visible
    And creation script text in "demog" project table in "Save project" dialog should contain text ":BDDLifeScript{time}()"
    When user enters "BDDLifeScriptProj{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And 1 project named "BDDLifeScriptProj{time}" should be on the server
    And the "demog" table of the "BDDLifeScriptProj{time}" project should be saved with data sync
    And "Share BDDLifeScriptProj{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDLifeScriptProj{time}" dialog
    Then the "Share BDDLifeScriptProj{time}" dialog should close
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The project's share reaches the script too
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeScriptProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user picks "Share..." from the context menu of BDDLifeScriptProj{time} gallery card
    Then "Share BDDLifeScriptProj{time}" dialog should be visible
    # the dialog fetches the project's grants after it opens; OK before that throws "Not initialized"
    And "Share BDDLifeScriptProj{time}" dialog should contain text "Full access"
    When user picks the sharing user in "User, group, or email" input in "Share BDDLifeScriptProj{time}" dialog
    Then share access selector should contain text "View and use"
    When user clicks on OK button in "Share BDDLifeScriptProj{time}" dialog
    Then the "Share BDDLifeScriptProj{time}" dialog should close
    Given the context panel is open
    When user clicks on BDDLifeScriptProj{time} gallery card
    Then the context panel should show "BDDLifeScriptProj{time}"
    And the sharing pane should list the sharing user
    Given user opens the Scripts view
    When user clears gallery search
    And user types "BDDLifeScript{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user picks "Share..." from the context menu of "BDDLifeScript{time}" link in gallery
    Then "Share BDDLifeScript{time}" dialog should be visible
    And the access level of the sharing user in "Share BDDLifeScript{time}" dialog should be "View and use"
    When user clicks on CANCEL button in "Share BDDLifeScript{time}" dialog
    Then the "Share BDDLifeScript{time}" dialog should close
    When user clears gallery search
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The second account gets the data from the shared project (GROK-19403)
    Given user signs in as the sharing user
    And the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeScriptProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDLifeScriptProj{time} gallery card
    Then the current view should be a TableView view
    And the "demog" view should be current
    And the table should have 5850 rows
    And the table should have been reloaded by data sync
    And "Data loading error" dialog should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Renamed in its editor, the script still feeds the project
    Given user signs in as themselves again
    And the browse panel is open
    Given user opens the Scripts view
    When user clears gallery search
    And user types "BDDLifeScript{time}" into gallery search
    And user picks "Edit..." from the context menu of "BDDLifeScript{time}" link in gallery
    Then the "BDDLifeScript{time}" view should be current
    When user replaces the code of code editor with "//name: BDDLifeScriptRenamed{time}"
    And user appends "//language: javascript" to code editor
    And user appends "//output: dataframe df" to code editor
    And user appends "df = await grok.data.getDemoTable('demog.csv');" to code editor
    And user saves the script
    Then 1 script named "BDDLifeScriptRenamed{time}" should be on the server
    And 0 scripts named "BDDLifeScript{time}" should be on the server
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeScriptProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDLifeScriptProj{time} gallery card
    Then the current view should be a TableView view
    And the "demog" view should be current
    And the table should have 5850 rows
    And the table should have been reloaded by data sync
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Renamed, the project opens under its new name
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeScriptProj{time}" into gallery search
    And user picks "Rename..." from the context menu of BDDLifeScriptProj{time} gallery card
    Then "Rename project" dialog should be visible
    When user enters "BDDLifeScriptProjRenamed{time}" into Name input in "Rename project" dialog
    And user clicks on OK button in "Rename project" dialog
    Then the "Rename project" dialog should close
    And 1 project named "BDDLifeScriptProjRenamed{time}" should be on the server
    And 0 projects named "BDDLifeScriptProj{time}" should be on the server
    When user enters "BDDLifeScriptProjRenamed{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDLifeScriptProjRenamed{time} gallery card
    Then the current view should be a TableView view
    And the "demog" view should be current
    And the table should have 5850 rows
    And the table should have been reloaded by data sync
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The script is broken to throw
    When user closes all views
    Given user opens the Scripts view
    When user clears gallery search
    And user types "BDDLifeScriptRenamed{time}" into gallery search
    And user picks "Edit..." from the context menu of "BDDLifeScriptRenamed{time}" link in gallery
    Then the "BDDLifeScriptRenamed{time}" view should be current
    When user replaces the code of code editor with "//name: BDDLifeScriptRenamed{time}"
    And user appends "//language: javascript" to code editor
    And user appends "//output: dataframe df" to code editor
    And user appends "throw new Error('intentional break');" to code editor
    And user appends "df = await grok.data.getDemoTable('demog.csv');" to code editor
    And user saves the script
    Then the script "BDDLifeScriptRenamed{time}" on the server should contain "throw new Error('intentional break');"
    And no errors should have been logged

  Scenario: The second account gets no way to edit the broken script (GROK-19728)
    Given user signs in as the sharing user
    And the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeScriptProjRenamed{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDLifeScriptProjRenamed{time} gallery card
    Then "Data loading error" dialog should be visible
    And "Data loading error" dialog should contain text "Error: intentional break"
    And "Data loading error" dialog should contain text "Ask the project owner to fix the script."
    And the following elements should be visible:
      | OPEN ANYWAY button in "Data loading error" dialog   |
      | CLOSE PROJECT button in "Data loading error" dialog |
    And "EDIT SCRIPT..." button in "Data loading error" dialog should be absent
    When user clicks on CLOSE PROJECT button in "Data loading error" dialog
    Then the "Data loading error" dialog should close

  Scenario: The owner is offered OPEN ANYWAY, EDIT SCRIPT... and CLOSE PROJECT
    Given user signs in as themselves again
    And the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeScriptProjRenamed{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDLifeScriptProjRenamed{time} gallery card
    Then "Data loading error" dialog should be visible
    And "Data loading error" dialog should contain text "Project \"BDDLifeScriptProjRenamed{time}\" could not load some of its data."
    And "Data loading error" dialog should contain text "Error: intentional break"
    And "Data loading error" dialog should contain text "edit the script to fix it"
    And the following elements should be visible:
      | OPEN ANYWAY button in "Data loading error" dialog      |
      | "EDIT SCRIPT..." button in "Data loading error" dialog |
      | CLOSE PROJECT button in "Data loading error" dialog    |
    When user clicks on CLOSE PROJECT button in "Data loading error" dialog
    Then the "Data loading error" dialog should close

  Scenario: Delete Project and the script's Delete remove both
    When user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeScriptProjRenamed{time}" into gallery search
    And user picks "Delete Project" from the context menu of BDDLifeScriptProjRenamed{time} gallery card
    Then "Are you sure?" dialog should be visible
    And "Are you sure?" dialog should contain text "Delete project \"BDDLifeScriptProjRenamed{time}\"?"
    When user clicks on DELETE button in "Are you sure?" dialog
    # the dialog stays open, its button disabled, until the server has deleted the project
    Then the "Are you sure?" dialog should close
    And 0 projects named "BDDLifeScriptProjRenamed{time}" should be on the server
    Given user opens the Scripts view
    When user clears gallery search
    And user types "BDDLifeScriptRenamed{time}" into gallery search
    And user picks "Delete" from the context menu of "BDDLifeScriptRenamed{time}" link in gallery
    Then "Are you sure?" dialog should be visible
    And "Are you sure?" dialog should contain text "Delete script \"BDDLifeScriptRenamed{time}\"?"
    When user clicks on YES button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 scripts named "BDDLifeScriptRenamed{time}" should be on the server
    And "BDDLifeScriptRenamed{time}" link in gallery should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown
    # the gallery keeps its search text for the next visit: leave it empty
    When user clears gallery search
