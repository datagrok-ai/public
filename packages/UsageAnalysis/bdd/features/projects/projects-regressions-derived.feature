@serial @realizes:views.projects
Feature: Projects regressions: derived and linked tables in projects
  Regression guards for fixed Jira bugs about tables made from other tables: GROK-19664 (a pivot
  published from a second table of the same name goes into the other open project), GROK-19580
  (Save as Zip of a data-sync project with a published pivot) and GROK-19135 (linking works only in
  the first of two open projects with the same linked tables). GROK-20013 (a copy saved after an
  in-place join) is in projects-regressions-join.feature. Saves go through the ribbon's Save dialog,
  links through the Data menu's dialog, reopens, Save as Zip and deletes through the Dashboards
  gallery. demog stands for SPGI in every report. Where the ticket opens its tables from files
  (GROK-19664, GROK-19580), they are opened by a double-click in Browse > Files > Demo; the second
  demog opened that way keeps the name "demog", the same name as in the ticket.

  What each scenario fails on: the pivot listed under the saved project in the Dashboards panel
  (GROK-19664); no zip downloaded, a zip without the pivot table, or an error (GROK-19580); the
  second project's link doing nothing, or moving the first project's tables (GROK-19135).

  In GROK-19135 the rows are selected through the JS API: both open projects hold tables named
  "demog" and "demog (2)" in views of the same names, so no grid gesture can name the one it means;
  the claim is what the link does with the selection. The ticket's two projects use two link types,
  and so do these (selection to selection, selection to filter).

  The Save dialog's preview logs "Unable to find element in cloned iframe" for some views
  (GROK-18606, won't fix, known noise by the operator's ruling): the check right after a save lets
  that one message through; every check after a reopen is strict. Every project is named with the
  run's time and removed with its tables and views when each scenario starts and when the feature
  ends.

  Background:
    Given user is logged in
    And the browse panel is open

  @realizes:GROK-19664
  Scenario: A pivot of a second, same-named table goes into the new project, not the saved one
    Given no project named "BDDRegFirst{time}" is on the server
    And user clears the saved pivot table parameters
    Given the browse panel is open
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---demog.csv tree node inside browse tree
    Then the "demog" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDRegFirst{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegFirst{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegFirst{time}" dialog
    Then the "Share BDDRegFirst{time}" dialog should close
    And no errors but the project preview's should have been logged
    Given the browse panel is open
    And user double-clicks Files---Demo---demog.csv tree node inside browse tree
    Then the open tables should be exactly "demog, demog"
    And the "demog" view should be current
    Given the toolbox pane is shown
    When user clicks on "pivot table" icon in toolbox
    Then pivot table viewer should be visible
    When user clicks on the "add to workspace" area of pivot table viewer
    Then the "demog aggregation" view should be current
    Given the dashboards panel of the left sidebar is open
    Then there should be 2 visible dashboards project node
    And New-Dashboard---demog-aggregation tree node inside browse tree should be visible
    Given BDDRegFirst{time} tree node inside browse tree is expanded
    Then BDDRegFirst{time}---demog tree node inside browse tree should be visible
    And BDDRegFirst{time}---demog-aggregation tree node inside browse tree should be absent
    And no errors should have been logged

  @realizes:GROK-19580
  Scenario: A data-sync project with a published pivot is saved as a zip from the gallery
    Given no project named "BDDRegZip{time}" is on the server
    And user clears the saved pivot table parameters
    Given the browse panel is open
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---demog.csv tree node inside browse tree
    Then the "demog" view should be current
    Given the toolbox pane is shown
    When user clicks on "pivot table" icon in toolbox
    Then pivot table viewer should be visible
    When user clicks on the "add to workspace" area of pivot table viewer
    Then the "demog aggregation" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Data sync switch in "demog" project table in "Save project" dialog should be checked
    And Data sync switch in "demog aggregation" project table in "Save project" dialog should be checked
    When user enters "BDDRegZip{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegZip{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegZip{time}" dialog
    Then the "Share BDDRegZip{time}" dialog should close
    And the "demog aggregation" table of the "BDDRegZip{time}" project should be saved with data sync
    And no errors but the project preview's should have been logged
    When user closes all views
    And the browse panel is open
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegZip{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    Given user watches downloads
    When user picks "Save as Zip" from the context menu of BDDRegZip{time} gallery card
    Then a file "BDDRegZip{time}.zip" should have been downloaded
    And the downloaded file "BDDRegZip{time}.zip" should contain the text "DemogAggregation = Aggregate("
    And no error or warning balloon should have been shown
    And no errors should have been logged

  @realizes:GROK-19135
  Scenario: Two open projects with the same linked tables each keep their own link
    Given no project named "BDDRegLinkA{time}" is on the server
    And no project named "BDDRegLinkB{time}" is on the server
    And user opens demog dataset
    And user opens demog dataset
    When user picks "Data > Link Tables..." from the top menu
    Then "Link Tables" dialog should be visible
    When user selects "selection to selection" in "Link Type" input in "Link Tables" dialog
    And user clicks on LINK button in "Link Tables" dialog
    And user clicks on CLOSE button in "Link Tables" dialog
    Then the "Link Tables" dialog should close
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDRegLinkA{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegLinkA{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegLinkA{time}" dialog
    Then the "Share BDDRegLinkA{time}" dialog should close
    And no errors but the project preview's should have been logged
    When user closes all views
    Given user opens demog dataset
    And user opens demog dataset
    When user picks "Data > Link Tables..." from the top menu
    Then "Link Tables" dialog should be visible
    When user selects "selection to filter" in "Link Type" input in "Link Tables" dialog
    And user clicks on LINK button in "Link Tables" dialog
    And user clicks on CLOSE button in "Link Tables" dialog
    Then the "Link Tables" dialog should close
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDRegLinkB{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDRegLinkB{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDRegLinkB{time}" dialog
    Then the "Share BDDRegLinkB{time}" dialog should close
    And no errors but the project preview's should have been logged
    When user closes all views
    And the browse panel is open
    And user clicks on Dashboards tree node inside browse tree
    And user enters "BDDRegLink" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    Given user watches the task bar
    When user double-clicks on BDDRegLinkA{time} gallery card
    Then the task bar should have finished "Opening project"
    Then the "demog" table of the open "BDDRegLinkA{time}" project should have 0 selected rows
    Given the browse panel is open
    And user clicks on Dashboards tree node inside browse tree
    Given user watches the task bar
    When user double-clicks on BDDRegLinkB{time} gallery card
    Then the task bar should have finished "Opening project"
    Then the "demog (2)" table of the open "BDDRegLinkB{time}" project should have 5850 rows passing the filter
    When user selects the first 10 rows of the "demog" table of the open "BDDRegLinkB{time}" project
    Then the "demog (2)" table of the open "BDDRegLinkB{time}" project should have 10 rows passing the filter
    And the "demog (2)" table of the open "BDDRegLinkA{time}" project should have 0 selected rows
    When user selects the first 10 rows of the "demog" table of the open "BDDRegLinkA{time}" project
    Then the "demog (2)" table of the open "BDDRegLinkA{time}" project should have 10 selected rows
    And the "demog (2)" table of the open "BDDRegLinkB{time}" project should have 10 rows passing the filter
    And no error or warning balloon should have been shown
    And no errors should have been logged
