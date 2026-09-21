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

  Not translated here, and each for its own reason. Browse-Nav-06 has two halves: that a newly
  created object appears after Refresh, and that the expanded set survives it (GROK-16261). Only
  the second needs a signal Refresh does not give — the icon's handler awaits the reload and then
  chains the path parse without awaiting it, and nothing observable marks either boundary — while
  the first is writable today and is simply not written yet. Browse-Nav-09 walks every icon of the header; the
  scenarios below cover Home, Open local file, Open text, Collapse tree and Find path, and leave
  Refresh and the panel's own close icon uncovered.

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

  Scenario: Find path reveals the node of the object that is open
    Given Apps tree node inside browse tree is expanded
    When user clicks on Tutorials tree node inside browse tree
    And the browse panel is open
    And user clicks on "Collapse tree" icon inside browse toolbar
    Then Tutorials tree node inside browse tree should be hidden
    When user clicks on "Find path" icon inside browse toolbar
    Then Tutorials tree node inside browse tree should be visible
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
    And no errors should have been logged
    And no error or warning balloon should have been shown
