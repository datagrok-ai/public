@journey @serial @realizes:views.projects @realizes:sharing.share-dialog
Feature: A project built from the user's own script, shared, then the script renamed and broken
  A JavaScript script is written in a new editor (Scripts > New > JavaScript Script...) and saved;
  Run... from the Scripts gallery opens its demog output, which is saved as a project with Data
  sync on and the script as its creation script. The project alone is shared with the second
  account from its Dashboards card. The script is then renamed in its editor with a new body
  (demog(100)) and broken with a throw; Delete Project and the script's Delete remove both.
  Translated from the TestTrack case Projects/projects-lifecycle-script.

  The second account signs in on the feature's page (the library's sharing user) for the
  recipient's opens: GROK-19403 (it gets the data) and GROK-19728 (for the broken script it is told
  to ask the owner and gets no EDIT SCRIPT...). After the rename the page is reloaded before the
  owner's reopen, as the md does: the session that saved the script keeps running its old body for
  a while; the owner's open of the broken script follows signing back in, which loads a fresh page.

  The script runs in the page (grok.data reads demog), so no script container is involved. The
  editor has no line-replacing gesture, so each edit types the whole text again (with the new name
  and body, or with the throw before the df line). The script and the project are named with the
  run's time (letters and digits only: the Dashboards search misses "-" and "_") and removed, the
  script under both its names, when the feature starts and ends. It is serial: the Dashboards
  search and the Scripts view's search are shared with the features that use them.

  Background:
    Given user is logged in
    And the browse panel is open
    And no project named "BDDLifeScriptProj{time}" is on the server
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
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And "Creation script" button in "demog" project table in "Save project" dialog should be visible
    And Data sync switch in "demog" project table in "Save project" dialog should be checked
    And "Save project" dialog should contain text "Some tables require this script for data sync."
    When user clicks on "Creation script" button in "demog" project table in "Save project" dialog
    Then "demog" project table in "Save project" dialog should contain text ":BDDLifeScript{time}()"
    When user enters "BDDLifeScriptProj{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDLifeScriptProj{time}" uploaded' should have been shown
    And 1 project named "BDDLifeScriptProj{time}" should be on the server
    And the "demog" table of the "BDDLifeScriptProj{time}" project should be saved with data sync
    And "Share BDDLifeScriptProj{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDLifeScriptProj{time}" dialog
    Then the "Share BDDLifeScriptProj{time}" dialog should close
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Only the project is shared with the second account
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeScript" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user picks "Share..." from the context menu of BDDLifeScriptProj{time} gallery card
    Then "Share BDDLifeScriptProj{time}" dialog should be visible
    And share access selector should contain text "View and use"
    When user picks the sharing user in "User, group, or email" input in "Share BDDLifeScriptProj{time}" dialog
    And user clicks on OK button in "Share BDDLifeScriptProj{time}" dialog
    Then the "Share BDDLifeScriptProj{time}" dialog should close
    And an info balloon containing "Shared" should have been shown
    When user clicks on BDDLifeScriptProj{time} gallery card
    Then the context panel should show "BDDLifeScriptProj{time}"
    And the sharing pane should list the sharing user

  Scenario: The second account gets the data of the shared project (GROK-19403)
    When user picks "Close All" from the context menu of browse tab
    And user signs in as the sharing user
    Then the sharing user should be signed in
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeScriptProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDLifeScriptProj{time} gallery card
    Then the "demog" view should be current
    And the table should have 5850 rows
    And the table should have been reloaded by data sync
    And "Data loading error" dialog should be absent
    And no error or warning balloon should have been shown
    When user picks "Close All" from the context menu of browse tab
    And user signs in as themselves again
    Then the running account should be signed in

  Scenario: The script is renamed in its editor with a new body
    Then the running account should be signed in
    When user picks "Close All" from the context menu of browse tab
    Given the browse panel is open
    And user opens the Scripts view
    When user clears gallery search
    And user types "BDDLifeScript{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user picks "Edit..." from the context menu of "BDDLifeScript{time}" link in gallery
    Then the "BDDLifeScript{time}" view should be current
    When user replaces the code of code editor with "//name: BDDLifeScriptRenamed{time}"
    And user appends "//language: javascript" to code editor
    And user appends "//output: dataframe df" to code editor
    And user appends "df = grok.data.demo.demog(100);" to code editor
    And user saves the script
    Then 1 script named "BDDLifeScriptRenamed{time}" should be on the server
    And 0 scripts named "BDDLifeScript{time}" should be on the server
    And the script "BDDLifeScriptRenamed{time}" on the server should contain "grok.data.demo.demog(100)"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: After a reload, the project runs the renamed script's new body
    When user reloads the page
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeScriptProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDLifeScriptProj{time} gallery card
    Then the "demog" view should be current
    And the table should have 100 rows
    And the table should have been reloaded by data sync
    And "Data loading error" dialog should be absent
    And no error or warning balloon should have been shown

  Scenario: The script is broken with a throw before its df line
    When user picks "Close All" from the context menu of browse tab
    Given the browse panel is open
    And user opens the Scripts view
    When user clears gallery search
    And user types "BDDLifeScriptRenamed{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user picks "Edit..." from the context menu of "BDDLifeScriptRenamed{time}" link in gallery
    Then the "BDDLifeScriptRenamed{time}" view should be current
    When user replaces the code of code editor with "//name: BDDLifeScriptRenamed{time}"
    And user appends "//language: javascript" to code editor
    And user appends "//output: dataframe df" to code editor
    And user appends "throw new Error('intentional break');" to code editor
    And user appends "df = grok.data.demo.demog(100);" to code editor
    And user saves the script
    Then the script "BDDLifeScriptRenamed{time}" on the server should contain "throw new Error('intentional break');"
    And no errors should have been logged

  Scenario: The second account is told to ask the owner and gets no EDIT SCRIPT... (GROK-19728)
    When user picks "Close All" from the context menu of browse tab
    And user signs in as the sharing user
    Then the sharing user should be signed in
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeScriptProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDLifeScriptProj{time} gallery card
    Then "Data loading error" dialog should be visible
    And "Data loading error" dialog should contain text "Ask the project owner to fix the script"
    And "OPEN ANYWAY" button in "Data loading error" dialog should be visible
    And "CLOSE PROJECT" button in "Data loading error" dialog should be visible
    And "EDIT SCRIPT..." button in "Data loading error" dialog should be absent
    When user clicks on "CLOSE PROJECT" button in "Data loading error" dialog
    Then the "Data loading error" dialog should close
    When user signs in as themselves again
    Then the running account should be signed in

  Scenario: The owner is offered OPEN ANYWAY, EDIT SCRIPT... and CLOSE PROJECT
    Then the running account should be signed in
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeScriptProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDLifeScriptProj{time} gallery card
    Then "Data loading error" dialog should be visible
    And "Data loading error" dialog should contain text "could not load some of its data"
    And "Data loading error" dialog should contain text "Error: intentional break"
    And "OPEN ANYWAY" button in "Data loading error" dialog should be visible
    And "EDIT SCRIPT..." button in "Data loading error" dialog should be visible
    And "CLOSE PROJECT" button in "Data loading error" dialog should be visible
    When user clicks on "CLOSE PROJECT" button in "Data loading error" dialog
    Then the "Data loading error" dialog should close

  Scenario: Delete Project and the script's Delete remove both
    When user picks "Close All" from the context menu of browse tab
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLifeScriptProj{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user picks "Delete Project" from the context menu of BDDLifeScriptProj{time} gallery card
    Then "Are you sure?" dialog should be visible
    And "Are you sure?" dialog should contain text 'Delete project "BDDLifeScriptProj{time}"?'
    When user clicks on DELETE button in "Are you sure?" dialog
    # the dialog stays open, its button disabled, until the server has deleted the project
    Then the "Are you sure?" dialog should close
    And 0 projects named "BDDLifeScriptProj{time}" should be on the server
    Given user opens the Scripts view
    When user clears gallery search
    And user types "BDDLifeScriptRenamed{time}" into gallery search
    And user picks "Delete" from the context menu of "BDDLifeScriptRenamed{time}" link in gallery
    Then "Are you sure?" dialog should be visible
    And "Are you sure?" dialog should contain text 'Delete script "BDDLifeScriptRenamed{time}"?'
    When user clicks on YES button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 scripts named "BDDLifeScriptRenamed{time}" should be on the server
    And "BDDLifeScriptRenamed{time}" link in gallery should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown
    # the gallery keeps its search text for the next visit: leave it empty
    When user clears gallery search
