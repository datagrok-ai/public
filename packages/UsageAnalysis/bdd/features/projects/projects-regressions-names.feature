@serial @realizes:views.projects
Feature: Projects regressions: names of projects and their tables
  Regression guards for fixed Jira bugs about names: GROK-19726 (a table renamed to digits only),
  GROK-19788 and GROK-18472 (a project name with "|" and "_", renamed from the Dashboards gallery),
  GROK-15135 and GROK-20197 (a name reused after the project was deleted) and GROK-17700 (saving a
  data-sync project renamed in the gallery). GROK-19792 (copies saved without renaming) is in
  projects-regressions-copies.feature. Saves go through the ribbon's Save dialog; renames, deletes
  and reopens through the Dashboards gallery. demog and cars stand for the tables of the reports; a
  table whose Data sync matters is opened from Browse > Files > Demo, the others through the JS API.

  What each scenario fails on: the digits-named table not opening or erroring (GROK-19726); the
  Rename dialog erroring or refusing "_" (GROK-19788, GROK-18472); the reused name bringing back
  the deleted project's table or a "_1" in its address (GROK-15135, GROK-20197); a Save dialog with
  no table, an error on OK, or the save landing under the old name (GROK-17700).

  The console is claimed clean after every save as well (the Save dialog's preview noise, GROK-18606,
  did not show for these single-table views). Every project is named with the run's time and removed
  with its tables and views when each scenario starts and when the feature ends. The gallery's Delete
  Project removes the project alone: the table and view of the project GROK-15135 deletes stay on the
  server until the library sweeps them (requested in the request document). That the reused name brings back
  only the new table is read in the reopened workspace (no demog view); the server-side claim on the
  project's tables is parked in the request document.

  Background:
    Given user is logged in
    And the browse panel is open

  @realizes:GROK-19726
  Scenario Outline: A project whose table was renamed to digits opens again, with Data sync <mode>
    Given no project named "BDDRegDigits<mode>{time}" is on the server
    And simple mode is off
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---demog.csv tree node inside browse tree
    Then the "demog" view should be current
    When user clicks on browse tab
    Then the browse tree should be visible
    When user picks "Table > Rename..." from the context menu of demog view
    Then "Rename table" dialog should be visible
    When user enters "12345" into "New name:" text input in "Rename table" dialog
    And user clicks on OK button in "Rename table" dialog
    Then the "Rename table" dialog should close
    And table "12345" should be open
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDRegDigits<mode>{time}" into Name text input in "Save project" dialog
    And user <switch> Data sync switch in "12345" project table in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegDigits<mode>{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegDigits<mode>{time}" dialog
    Then the "Share BDDRegDigits<mode>{time}" dialog should close
    And the "12345" table of the "BDDRegDigits<mode>{time}" project should be saved <kept>
    And no errors should have been logged
    When user closes all views
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegDigits<mode>{time}" into gallery search
    And user double-clicks on BDDRegDigits<mode>{time} gallery card
    Then the "12345" view should be current
    And the table should have 5850 rows
    And the table should have been <loaded>
    And no error or warning balloon should have been shown
    And no errors should have been logged

    Examples:
      | mode | switch   | kept           | loaded               |
      | On   | checks   | with data sync | reloaded by data sync |
      | Off  | unchecks | as a snapshot  | loaded as a snapshot |

  @realizes:GROK-19788 @realizes:GROK-18472
  Scenario: A project named with "|" and "_" is renamed from the gallery to a name with "_"
    Given no project named "BDDRegPipe{time}|A_B" is on the server
    And no project named "BDDRegPiped{time}_x" is on the server
    And user opens demog dataset
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDRegPipe{time}|A_B" into Name text input in "Save project" dialog
    Then OK button in "Save project" dialog should be enabled
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegPipe{time}|A_B" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegPipe{time}|A_B" dialog
    Then the "Share BDDRegPipe{time}|A_B" dialog should close
    And 1 project named "BDDRegPipe{time}|A_B" should be on the server
    And no errors should have been logged
    When user closes all views
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegPipe{time}" into gallery search
    And user picks "Rename..." from the context menu of BDDRegPipe{time}AB gallery card
    Then "Rename project" dialog should be visible
    And Name text input in "Rename project" dialog should have value "BDDRegPipe{time}|A_B"
    And no errors should have been logged
    When user enters "BDDRegPiped{time}_x" into Name text input in "Rename project" dialog
    Then OK button in "Rename project" dialog should be enabled
    When user clicks on OK button in "Rename project" dialog
    Then the "Rename project" dialog should close
    And 1 project named "BDDRegPiped{time}_x" should be on the server
    And 0 projects named "BDDRegPipe{time}|A_B" should be on the server
    And no error or warning balloon should have been shown
    And no errors should have been logged

  @realizes:GROK-15135 @realizes:GROK-20197
  Scenario: A name reused after its project was deleted brings back only the new table, at the same address
    Given no project named "BDDRegReuse{time}" is on the server
    And user opens demog dataset
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDRegReuse{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegReuse{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegReuse{time}" dialog
    Then the "Share BDDRegReuse{time}" dialog should close
    And no errors should have been logged
    And 1 project named "BDDRegReuse{time}" should be on the server
    When user closes all views
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegReuse{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user picks "Delete Project" from the context menu of BDDRegReuse{time} gallery card
    Then "Are you sure?" dialog should be visible
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 projects named "BDDRegReuse{time}" should be on the server
    Given user opens cars dataset
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDRegReuse{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegReuse{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegReuse{time}" dialog
    Then the "Share BDDRegReuse{time}" dialog should close
    And no errors should have been logged
    When user closes all views
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegReuse{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDRegReuse{time} gallery card
    Then the "cars" view should be current
    And table "cars" should be open
    And demog view should be absent
    And the page address should contain ".bddregreuse{time}/"
    And the page address should not contain "bddregreuse{time}_"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  @realizes:GROK-17700
  Scenario: A data-sync project renamed in the gallery is saved again with its table
    Given no project named "BDDRegRename{time}" is on the server
    And no project named "BDDRegRenamed{time}" is on the server
    And simple mode is off
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---demog.csv tree node inside browse tree
    Then the "demog" view should be current
    When user clicks on browse tab
    Then the browse tree should be visible
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And "Creation script" button in "demog" project table in "Save project" dialog should be visible
    And Data sync switch in "demog" project table in "Save project" dialog should be checked
    When user enters "BDDRegRename{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegRename{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegRename{time}" dialog
    Then the "Share BDDRegRename{time}" dialog should close
    And no errors should have been logged
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegRename{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user picks "Rename..." from the context menu of BDDRegRename{time} gallery card
    Then "Rename project" dialog should be visible
    When user enters "BDDRegRenamed{time}" into Name text input in "Rename project" dialog
    And user clicks on OK button in "Rename project" dialog
    Then the "Rename project" dialog should close
    And 1 project named "BDDRegRenamed{time}" should be on the server
    When user clicks on demog view
    Then the "demog" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And "demog" project table in "Save project" dialog should be visible
    And Data sync switch in "demog" project table in "Save project" dialog should be checked
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDRegRenamed{time}" uploaded' should have been shown
    And no error or warning balloon should have been shown
    And no errors should have been logged
    And 1 project named "BDDRegRenamed{time}" should be on the server
    And 0 projects named "BDDRegRename{time}" should be on the server
    And the "demog" table of the "BDDRegRenamed{time}" project should be saved with data sync
    When user closes all views
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegRenamed{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDRegRenamed{time} gallery card
    Then the "demog" view should be current
    And the table should have 5850 rows
    And the table should have been reloaded by data sync
    And no error or warning balloon should have been shown
    And no errors should have been logged
