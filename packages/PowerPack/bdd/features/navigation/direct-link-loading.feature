@journey @serial
Feature: A project opened by its direct link
  A saved project opens when the page is loaded on the project's direct link — the address the
  platform gives it, `/p/<namespace>.<name>` — as pasting the link into the address bar does
  (GROK-18721). The project's table view
  becomes the current view with its grid, columns and rows, the Home page yields to it, and nothing
  is logged and no balloon shown. The control case opens the same project from inside the platform,
  through the Browse tree, with the same claims. Translated from TestTrack
  PowerPack/direct-link-loading.md.

  Not translated: how the loading window looks while the page loads (its size, cropping and the
  layout under it) is a claim about pixels only, and stays a manual check. The page is loaded on
  the link in the same browser the feature runs in, not in a fresh profile: what the browser has
  cached from the earlier load stays.

  It runs @serial: a page loading while home-widgets.feature hides or shows a widget can write the
  widget settings it read at its start back to the server (see that feature).

  Background:
    Given user is logged in

  Scenario: A project with the demog table is saved from the toolbar
    Given user opens demog dataset
    And no project named "bdd-direct-link-{time}" is on the server
    And the browse panel is open
    When user clicks on Save button in toolbar
    Then "Save project" dialog should be visible
    When user enters "bdd-direct-link-{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then "Save project" dialog should be hidden
    And 1 project named "bdd-direct-link-{time}" should be on the server
    And no error or warning balloon should have been shown

  Scenario: The direct link opens the project's table view
    When user closes all views
    Then the project "bdd-direct-link-{time}" should not be open
    When user loads the direct link of project "bdd-direct-link-{time}"
    Then the project "bdd-direct-link-{time}" should be open
    And the "demog" view should be current
    And the table should have 5850 rows
    And the table should have a column "AGE"
    And the table should have a column "RACE"
    And grid should show 5850 rows
    And the "Home" view should not be current
    And no loading indicator should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The same project opened from the Browse tree shows the same view
    When user closes all views
    Then the "Home" view should be current
    And the project "bdd-direct-link-{time}" should not be open
    Given the browse panel is open
    When user refreshes the browse tree
    And user expands "My stuff" tree node inside browse tree
    And user double-clicks on "My stuff > bdd-direct-link-{time}" tree node inside browse tree
    Then the project "bdd-direct-link-{time}" should be open
    And the "demog" view should be current
    And the table should have 5850 rows
    And the table should have a column "AGE"
    And grid should show 5850 rows
    And no errors should have been logged
    And no error or warning balloon should have been shown
