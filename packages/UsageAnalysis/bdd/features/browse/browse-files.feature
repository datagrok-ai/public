@browse @realizes:views.browse
Feature: The Files section of the Browse tree
  What the Files section holds, what a folder opens as, and what a tabular file previews as.
  Translated from the manual cases Browse-Files-01, -02 and -03 (playwright-public/browse/
  files.test.ts, browse_manual_tests2.md section 7).

  A file preview does not become the shell's current table — `grok.shell.t` stays null — so the
  size of the preview is claimed through the grid's own reading rather than through the table.

  Browse-Files-04 (a file added on the server appears after Refresh, GROK-19844) is not translated:
  the old spec waived it in its own comment and asserted instead that the API sees what the API just
  wrote. A real claim needs Refresh to say when it has finished, and it says nothing.
  Browse-Files-05 (Shared with me grouped by sharer) depends on what other people happen to have
  shared with this account. Browse-Files-06 (download) is writable — the old spec used a plain
  download event and needed no configured path — and is simply not written yet.

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
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A tabular file opens as a preview with its rows
    Given Files---Demo tree node inside browse tree is expanded
    When user clicks on Files---Demo---demog.csv tree node inside browse tree
    Then the "demog" view should be current
    And grid should show 5850 rows
    And no errors should have been logged
    And no error or warning balloon should have been shown
