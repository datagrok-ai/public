@browse @realizes:views.browse
Feature: The Browse panel and the icons of its toolbar
  What the Browse panel opens with, and what each icon of its header does. Translated from the
  manual cases Browse-Nav-01..08 (files/TestTrack sources: playwright-public/browse/nav.test.ts and
  browse_manual_tests2.md section 1).

  The icons are named by the tooltip the platform gives them (browse_panel.dart initRibbon), which
  is not what the manual cases call them: "Import file" is "Open local file", "Import text" is
  "Open text", "Collapse all" is "Collapse tree" and "Locate current object" is "Find path".

  Collapse is claimed on a child that was visible and stops being visible, not on a count of
  expanded arrows: the old spec allowed one arrow to remain and would have passed with a node
  left open.

  A node below the top level is named by its full tree path ("Files---Demo"), which is what the
  platform writes into its own `name` attribute: several sections carry a node called Demo.

  Refresh (Browse-Nav-06, Nav-09) is a gesture with an end: the tree reports itself rebuilt
  (`onBrowseTreeRefreshed`) and every group it reopened has its children again, so both halves of
  Nav-06 are claimed — the script saved on the server meanwhile appears, and an open section stays
  open — as is the preview that was current before the refresh (browse.md 6). A folder open inside
  that section comes back collapsed (the GROK-16261 shape); that half waits for a filed ticket.

  Browse-Nav-03 (the search box of the header) is not claimed: the PowerPack box carries no name —
  its only handle is a placeholder with quotes in it ("Search everywhere. Try "aspirin" or "7JZK"")
  that an element phrase cannot quote — and a search reports no end and no empty state.

  Background:
    Given user is logged in
    And the browse panel is open

  Scenario: The tree opens with every top-level section
    Then the following elements should be visible:
      | My stuff tree node inside browse tree  |
      | Spaces tree node inside browse tree    |
      | Apps tree node inside browse tree      |
      | Files tree node inside browse tree     |
      | Dashboards tree node inside browse tree |
      | Databases tree node inside browse tree |
      | Platform tree node inside browse tree  |
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The Browse tab toggles the panel
    When user clicks on browse tab
    Then the browse tree should be hidden
    When user clicks on browse tab
    Then the browse tree should be visible
    And Apps tree node inside browse tree should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The Home icon returns to the Home view
    # Home is already current after login, so the claim is made from another view: otherwise it
    # holds whether or not the icon does anything.
    Given user clicks on "Open text" icon inside browse toolbar
    And the "Import text" view should be current
    When user clicks on "Home" icon inside browse toolbar
    Then the "Home" view should be current
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Collapse tree closes an open section
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree should be visible
    When user clicks on "Collapse tree" icon inside browse toolbar
    Then Files---Demo tree node inside browse tree should be hidden
    And Files tree node inside browse tree should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  # the node was selected by the click that opened the file, so its selection is not Find path's to claim
  Scenario: Find path reveals the node of the object that is open
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user clicks on Files---Demo---demog.csv tree node inside browse tree
    Then the "demog" view should be current
    When user clicks on "Collapse tree" icon inside browse toolbar
    Then Files---Demo---demog.csv tree node inside browse tree should be hidden
    When user clicks on "Find path" icon inside browse toolbar
    Then Files---Demo---demog.csv tree node inside browse tree should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Open local file imports a CSV into a table view
    When user uploads "fixtures/browse-import.csv" through "Open local file" icon inside browse toolbar
    Then the "browse-import" view should be current
    And the table should have 5 rows
    And the table should have 3 columns
    And the table should have a column "population"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Open text opens the text import view
    When user clicks on "Open text" icon inside browse toolbar
    Then the "Import text" view should be current
    And code editor should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Refresh brings in a script saved on the server meanwhile
    Given no script named "BDD-Browse-Script-{run}" is on the server
    And My stuff tree node inside browse tree is expanded
    # a bucket loads its children when it opens: open before the script exists, it can learn of it only by Refresh
    And user expands the "Scripts" bucket of My stuff
    And a script "BDD-Browse-Script-{run}" is on the server:
      """
      //language: javascript
      let x = 1;
      """
    Then "BDD-Browse-Script-{run}" should not be listed in the "Scripts" bucket of My stuff
    When user refreshes the browse tree
    And user expands the "Scripts" bucket of My stuff
    Then "BDD-Browse-Script-{run}" should be listed in the "Scripts" bucket of My stuff
    And no errors should have been logged
    And no error or warning balloon should have been shown
    # the tree remembers what was open per user: put My stuff back as the other features expect it
    When user collapses My stuff tree node inside browse tree

  # The folder inside it (Files > Demo) comes back collapsed after Refresh, 2 of 2 on the local stand:
  # a candidate finding (GROK-16261 shape), not claimed here until it is walked by hand and filed
  Scenario: Refresh keeps an open section open
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree should be visible
    When user refreshes the browse tree
    Then Files tree node inside browse tree should be expanded
    And Files---Demo tree node inside browse tree should be visible
    And Files---App-Data tree node inside browse tree should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Refresh keeps the preview that is open
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user clicks on Files---Demo---demog.csv tree node inside browse tree
    Then the "demog" view should be current
    When user refreshes the browse tree
    Then the "demog" view should be current
    And demog view should be visible
    And grid should show 5850 rows
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The panel's close icon hides it and the Browse tab brings it back
    When user clicks on browse panel close icon
    Then the browse tree should be hidden
    When user clicks on browse tab
    Then the browse tree should be visible
    And Apps tree node inside browse tree should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown
