@journey @serial @realizes:views.projects @realizes:sharing.share-dialog
Feature: A project through its card in the Dashboards gallery
  demog.csv opened from Browse > Files > Demo is saved through the ribbon's Save dialog, and every
  later step goes through the project's card in Browse > Dashboards: the Share dialog, Rename, the
  Copy items, Add To Favorites (a submenu: Only for me, or a group), Save as Zip, reopening by a double click and Delete Project.
  Translated from the TestTrack case Projects/projects-ui-smoke.

  Parked (see the request document): the tag added in Context Panel > Details and the "#tag"
  search (the Add tag box is a bare input no element reaches); everything about the description:
  typing it into the Save dialog (its Description box is a bare text area no element reaches), the
  second project saved without one, the card and the Details pane showing it, and changing it by
  saving again; the tag and the description surviving a page reload (no reload step); Copy > ID
  (a UUID the feature cannot know without reading the project on the server); the filled star in
  the context panel's header after Add To Favorites. The words of the recipient's line in the
  Sharing pane ("has special permissions") are claimed on the pane as a whole.

  The recipient of the share is the library's sharing user (the md's "qa_playwright" is any user
  other than yourself). The Grok name, markup and URL are claimed for the admin account the suite
  runs as (namespace "Admin"). Names are letters and digits only (the Dashboards search misses "-"
  and "_") and carry the run's time. The project, under both names, is deleted
  when the feature starts and ends, with its table, view and grant; the favorite the scenario
  adds is removed by the scenario. It is serial: the save uploads a table the reopen reads back.

  Background:
    Given user is logged in
    And the browse panel is open
    And no project named "BDDSmoke{time}" is on the server
    And no project named "BDDSmokeRenamed{time}" is on the server

  Scenario: A file opened from its context menu is saved as a project
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user picks "Open" from the context menu of Files---Demo---demog.csv tree node inside browse tree
    Then the "demog" view should be current
    And the table should have 5850 rows
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And "Creation script" button in "demog" project table in "Save project" dialog should be visible
    And Data sync switch in "demog" project table in "Save project" dialog should be checked
    When user enters "BDDSmoke{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDSmoke{time}" uploaded' should have been shown
    And 1 project named "BDDSmoke{time}" should be on the server
    And "Share BDDSmoke{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDSmoke{time}" dialog
    Then the "Share BDDSmoke{time}" dialog should close
    And no errors should have been logged

  Scenario: The card is found by the gallery search
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    Then the "Projects" view should be current
    When user remembers the gallery counter
    And user enters "BDDSmoke{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    Then the gallery counter should be lower than remembered
    And BDDSmoke{time} gallery card should be visible

  Scenario: The project is shared from its card
    When user clicks on BDDSmoke{time} gallery card
    Then the context panel should show "BDDSmoke{time}"
    When user picks "Share..." from the context menu of BDDSmoke{time} gallery card
    Then "Share BDDSmoke{time}" dialog should be visible
    And share access selector should contain text "View and use"
    When user picks the sharing user in "User, group, or email" input in "Share BDDSmoke{time}" dialog
    And user unchecks "Send notifications" input in "Share BDDSmoke{time}" dialog
    When user clicks on OK button in "Share BDDSmoke{time}" dialog
    Then the "Share BDDSmoke{time}" dialog should close
    When user clicks on BDDSmoke{time} gallery card
    Then the context panel should show "BDDSmoke{time}"
    And the sharing pane should list the sharing user
    And Sharing pane in context panel should contain text "has special permissions"

  Scenario: The project is renamed from its card
    When user picks "Rename..." from the context menu of BDDSmoke{time} gallery card
    Then Rename project dialog should be visible
    And Name input in Rename project dialog should have value "BDDSmoke{time}"
    When user enters "BDDSmokeRenamed{time}" into Name input in Rename project dialog
    And user clicks on OK button in Rename project dialog
    Then the "Rename project" dialog should close
    And 1 project named "BDDSmokeRenamed{time}" should be on the server
    And 0 projects named "BDDSmoke{time}" should be on the server
    When user enters "BDDSmokeRenamed{time}" into gallery search
    And user clicks on Refresh icon in gallery toolbar
    Then BDDSmokeRenamed{time} gallery card should be visible

  Scenario: The Copy items put the Grok name, the markup and the URL on the clipboard
    When user picks "Copy > Grok name" from the context menu of BDDSmokeRenamed{time} gallery card
    Then the clipboard should have text "Admin:BDDSmokeRenamed{time}"
    When user picks "Copy > Markup" from the context menu of BDDSmokeRenamed{time} gallery card
    Then the clipboard should have text '#{x.Admin:BDDSmokeRenamed{time}."BDDSmokeRenamed{time}"}'
    When user picks "Copy > URL" from the context menu of BDDSmokeRenamed{time} gallery card
    Then the clipboard should contain text "/p/Admin.BDDSmokeRenamed{time}"

  Scenario: The project is added to favorites and taken out again
    Given "BDDSmokeRenamed{time}" is not in favorites
    And "My stuff" tree node inside browse tree is expanded
    And "My stuff > Favorites" tree node inside browse tree is expanded
    Then "My stuff > Favorites > BDDSmokeRenamed{time}" tree node inside browse tree should be absent
    When user picks "Add To Favorites > Only for me" from the context menu of BDDSmokeRenamed{time} gallery card
    Then "BDDSmokeRenamed{time}" should be in favorites on the server
    When user collapses "My stuff > Favorites" tree node inside browse tree
    And user expands "My stuff > Favorites" tree node inside browse tree
    Then "My stuff > Favorites > BDDSmokeRenamed{time}" tree node inside browse tree should be visible
    When user picks "Add To Favorites > Only for me" from the context menu of BDDSmokeRenamed{time} gallery card
    Then "BDDSmokeRenamed{time}" should not be in favorites on the server
    When user collapses "My stuff > Favorites" tree node inside browse tree
    And user expands "My stuff > Favorites" tree node inside browse tree
    Then "My stuff > Favorites > BDDSmokeRenamed{time}" tree node inside browse tree should be absent

  Scenario: The project is saved as a zip file
    Given user watches downloads
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDSmokeRenamed{time}" into gallery search
    And user picks "Save as Zip" from the context menu of BDDSmokeRenamed{time} gallery card
    Then a file "BDDSmokeRenamed{time}.zip" should have been downloaded
    And the downloaded file "BDDSmokeRenamed{time}.zip" should contain text "demog"

  Scenario: The card reopens the project with its data
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDSmokeRenamed{time}" into gallery search
    And user double-clicks on BDDSmokeRenamed{time} gallery card
    Then the "demog" view should be current
    And the table should have 5850 rows
    And status bar should contain text "5,850"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The project is deleted from its card
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDSmokeRenamed{time}" into gallery search
    Then BDDSmokeRenamed{time} gallery card should be visible
    When user remembers the gallery counter
    And user picks "Delete Project" from the context menu of BDDSmokeRenamed{time} gallery card
    Then "Are you sure?" dialog should contain text 'Delete project "BDDSmokeRenamed{time}"?'
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 projects named "BDDSmokeRenamed{time}" should be on the server
    When user clicks on Refresh icon in gallery toolbar
    Then the gallery counter should be lower than remembered
    And BDDSmokeRenamed{time} gallery card should be absent
