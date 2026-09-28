@serial @realizes:views.projects @realizes:file.menu.save.tables-as-project
Feature: Projects uploaded from an SDF file, a Get Top 100 result and the scratchpad
  The sources the automated uploading matrix leaves out: an SDF file of molecules saved with Data
  sync, a database table's Get Top 100 result saved with Data sync on and off, and the upload dialog
  the Dashboards panel opens from the New Dashboard row (the scratchpad), with a description and
  the Presentation mode switch. Each is saved, reopened from its Dashboards card and claimed by how
  its table came back. Translated from the TestTrack case Projects/uploading-ui.

  Not translated: "Project from two local files" — it goes through the operating system's file
  picker with files of the local machine, which a CI agent does not have.

  Fixtures substituted: Get Top 100 runs on System:Datagrok's public.entity_types instead of the
  Northwind orders the md names (NorthwindTest exists only on dev); the table has fewer than 100
  rows on every stand, so its count is remembered after the first run and compared on reopen. The
  md asks that "molecules render": the claim is that the reopened table's molecule column is typed
  Molecule, that its grid column resolved Chem's Molecule cell renderer, and that nothing is logged
  — the drawing itself is a picture, which the suite does not judge.

  A reopen is claimed to be a reopen: the tables open before the save are marked in memory, Close
  All is claimed to leave no table, and every table that comes back must be a frame without the
  mark.

  The scratchpad upload also puts the shell itself into presentation mode as soon as it is saved;
  the "back to design mode" link the shell shows is what the claims read, and what leaves the mode.
  Presentation mode is switched off again when the feature ends.

  Every project (with its tables and views) is named with the run's time and removed before its
  scenario starts and when the feature ends. The Save dialog's preview logs "Unable to find element
  in cloned iframe" (GROK-18606, known noise), so a save is claimed to log nothing else; the reopen
  is claimed to log nothing at all. @serial: the Dashboards search and the uploads are shared with
  every feature that saves a project.

  Background:
    Given user is logged in
    And the browse panel is open

  Scenario: mol1K.sdf from App Data is saved with Data sync and reopens with its molecules
    Given no project named "BDDUpXSdf{time}" is on the server
    And the package autostarts have completed
    And Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---Chem tree node inside browse tree is expanded
    When user double-clicks Files---App-Data---Chem---mol1K.sdf tree node inside browse tree
    Then the "mol1K" view should be current
    And table "mol1K" should have 1000 rows
    And "molecule" column should have semantic type "Molecule"
    When user marks the open tables as the frames in memory
    And user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "mol1K" project table in "Save project" dialog should be switched on
    When user enters "BDDUpXSdf{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "BDDUpXSdf{time}" uploaded' should have been shown
    And 1 project named "BDDUpXSdf{time}" should be on the server
    And the "mol1K" table of the "BDDUpXSdf{time}" project should be saved with data sync
    And "Share BDDUpXSdf{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDUpXSdf{time}" dialog
    Then the "Share BDDUpXSdf{time}" dialog should close
    And no errors but the project preview's should have been logged
    When user picks "Close All" from the context menu of left sidebar
    Then no table should be left in the workspace
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDUpXSdf{time}" into gallery search
    And user double-clicks on BDDUpXSdf{time} gallery card
    Then table "mol1K" should have been reloaded by data sync with 1000 rows
    And "molecule" column should have semantic type "Molecule"
    And the grid should draw the "molecule" column with the "Molecule" renderer
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario Outline: A Get Top 100 result of a database table is saved with Data sync <Sync>
    Given no project named "<Project>" is on the server
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok---Schemas---public tree node inside browse tree is expanded
    When user picks "Get Top 100" from the context menu of Databases---Postgres---Datagrok---Schemas---public---entity_types tree node inside browse tree
    Then the "entity_types" view should be current
    When user remembers the row count of table "entity_types" as "top 100"
    And user marks the open tables as the frames in memory
    And user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user clicks on "Creation script" button in "entity_types" project table in "Save project" dialog
    Then "entity_types" project table in "Save project" dialog should contain text "limit = 100"
    When user enters "<Project>" into Name text input in "Save project" dialog
    And user <switch> Data sync switch in "entity_types" project table in "Save project" dialog
    Then "Creation script" button in "entity_types" project table in "Save project" dialog should be <script>
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "<Project>" uploaded' should have been shown
    And 1 project named "<Project>" should be on the server
    And the "entity_types" table of the "<Project>" project should be saved <saved>
    And "Share <Project>" dialog should be visible
    When user clicks on CANCEL button in "Share <Project>" dialog
    Then the "Share <Project>" dialog should close
    And no errors but the project preview's should have been logged
    When user picks "Close All" from the context menu of left sidebar
    Then no table should be left in the workspace
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "<Project>" into gallery search
    And user double-clicks on <Project> gallery card
    Then table "entity_types" should have been <how> with the "top 100" row count
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | Sync | Project                | switch   | script  | saved          | how               |
      | ON   | BDDUpXTopSync{time}    | checks   | visible | with data sync | reloaded by data sync |
      | OFF  | BDDUpXTopNoSync{time}  | unchecks | hidden  | as a snapshot  | loaded as a snapshot |

  Scenario: The scratchpad's upload dialog saves a description and presentation mode
    Given no project named "BDDUpXScratch{time}" is on the server
    And presentation mode is off, now and when the feature ends
    And user opens demog dataset
    When user marks the open tables as the frames in memory
    Given the dashboards panel of the left sidebar is open
    Then New-Dashboard tree node inside browse tree should be visible
    When user clicks on new dashboard save button
    Then "Save project" dialog should be visible
    When user enters "BDDUpXScratch{time}" into Name text input in "Save project" dialog
    And user enters "Uploaded from the scratchpad" into project description field
    And user switches on Presentation mode switch in "Save project" dialog
    Then Presentation mode switch in "Save project" dialog should be switched on
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And 1 project named "BDDUpXScratch{time}" should be on the server
    And presentation back link should be visible
    And "Share BDDUpXScratch{time}" dialog should be visible
    When user clicks on OK button in "Share BDDUpXScratch{time}" dialog
    Then the "Share BDDUpXScratch{time}" dialog should close
    And no errors but the project preview's should have been logged
    When user clicks on presentation back link
    Then presentation back link should be absent
    When user picks "Close All" from the context menu of left sidebar
    Then no table should be left in the workspace
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDUpXScratch{time}" into gallery search
    Then BDDUpXScratch{time} gallery card should contain text "BDDUpXScratch{time}"
    And BDDUpXScratch{time} gallery card should contain text "Uploaded from the scratchpad"
    When user double-clicks on BDDUpXScratch{time} gallery card
    Then table "demog" should have been loaded as a snapshot with 5850 rows
    And presentation back link should be visible
    When user clicks on presentation back link
    Then presentation back link should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown
