@browse @realizes:views.browse
Feature: The Files section of the Browse tree
  What the Files section holds, what a folder opens as, and what a tabular file previews as.
  Translated from the manual cases Browse-Files-01, -02 and -03 (playwright-public/browse/
  files.test.ts, browse_manual_tests2.md section 7).

  A file preview does not become the shell's current table — `grok.shell.t` stays null — so the
  size of the preview is claimed through the grid's own reading rather than through the table.

  Browse-Files-04 (a file written on the server appears after Refresh, GROK-19844) is claimed with a
  file the feature writes into its own package's App Data folder and deletes again: Refresh reports
  when it is done, so the tree is read after it. Browse-Files-06 (download from a file's menu) is
  claimed on the file that arrives. Browse-Files-05 (Shared with me grouped by sharer) needs the
  second account to own and share something first, which no step can make it do.

  The shares claimed are the two every stand is provisioned with, App Data and Demo. The user's
  own home share is named by the stand ("My files" on dev) and a stand whose users have no home
  storage lists none, so it is not claimed.

  A node below the top level is named by its full tree path ("Files---Demo"), which is what the
  platform writes into its own `name` attribute. Several sections carry a node called Demo, Files
  or App Data, and a bare name matches whichever of them another feature happened to leave
  open: the tree remembers its expanded set per user, across features and across runs.

  Background:
    Given user is logged in
    And the browse panel is open
    And Files tree node inside browse tree is expanded

  Scenario: The Files section lists its file shares
    Then the following elements should be visible:
      | Files---App-Data tree node inside browse tree |
      | Files---Demo tree node inside browse tree     |
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A folder opens as a folder view of its own
    When user clicks on Files---Demo tree node inside browse tree
    Then the "Demo" view should be current
    # the title flips before anything is drawn, so the folder's contents are claimed too
    And gallery should contain text "chem"
    And the gallery counter should show as many items as the "System:DemoFiles/" folder holds on the server
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A tabular file opens as a preview with its rows
    Given Files---Demo tree node inside browse tree is expanded
    When user clicks on Files---Demo---demog.csv tree node inside browse tree
    Then the "demog" view should be current
    And grid should show 5850 rows
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Download from a file's menu hands over the file
    Given user watches downloads
    And Files---Demo tree node inside browse tree is expanded
    When user opens the context menu of Files---Demo---demog.csv tree node inside browse tree
    And user downloads a file through "Download" menu item in context menu
    Then a file "demog.csv" should have been downloaded
    And the downloaded file "demog.csv" should contain text "USUBJID"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A file written on the server shows in the tree after Refresh
    Given a file "System:AppData/UsageAnalysis/bdd-browse-refresh.txt" with text "written by a feature" is on the server
    And Files---App-Data tree node inside browse tree is expanded
    When user refreshes the browse tree
    # Refresh brings a folder open inside a section back collapsed (a candidate finding, see
    # browse-navigation.feature), so the path is opened again
    And Files---App-Data tree node inside browse tree is expanded
    And user expands Files---App-Data---UsageAnalysis tree node inside browse tree
    Then Files---App-Data---UsageAnalysis---bdd-browse-refresh.txt tree node inside browse tree should be visible
    When user collapses Files---App-Data tree node inside browse tree
    Then no errors should have been logged
    And no error or warning balloon should have been shown
