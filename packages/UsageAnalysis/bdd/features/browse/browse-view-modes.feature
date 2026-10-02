@browse @realizes:views.browse
Feature: Browsing mode and persistent views
  A single click on a tree node opens a view that the next single click replaces; a double click
  makes it persistent and later clicks leave it alone. Translated from the manual cases
  Browse-View-01 and Browse-View-04 (playwright-public/browse/view.test.ts, browse_manual_tests2.md
  section 3).

  Both sides of the rule are claimed every time — the view that must go and the view that must
  stay. A one-sided claim passes on a workspace that simply closes everything.

  The old spec resolved a nested Chem app at run time and skipped itself when none was deployed.
  Dashboards is used instead: it is a top-level node of the tree, present on every stand, so the
  case has no reason to skip.

  Browse-View-03 (the pin on a preview tab keeps the view) is claimed with simple mode off, since
  the pin lives on the view's tab handle, which simple mode hides. Browse-View-02 (an edit pins the
  preview) needs an edit gesture on a previewed view and is not written; Browse-View-05 (the badge
  counting the open views) has no text and no name: the count is a CSS `content: attr(data-count)`
  on a sidebar header simple mode keeps hidden.

  Background:
    Given user is logged in
    And the browse panel is open
    And Apps tree node inside browse tree is expanded

  Scenario: A single click replaces the view the previous single click opened
    When user clicks on Tutorials tree node inside browse tree
    Then the "Tutorials" view should be current
    And Tutorials view should be visible
    When user clicks on Dashboards tree node inside browse tree
    Then Projects view should be visible
    And Tutorials view should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A double click keeps the view through the next single click
    When user double-clicks on Tutorials tree node inside browse tree
    Then the "Tutorials" view should be current
    When user clicks on Dashboards tree node inside browse tree
    Then Projects view should be visible
    And Tutorials view should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A new browsing session does not unpin what is already persistent
    # the third node is a folder rather than another top-level group, so a nested node is covered too
    Given Files tree node inside browse tree is expanded
    When user double-clicks on Tutorials tree node inside browse tree
    Then Tutorials view should be visible
    When user clicks on Dashboards tree node inside browse tree
    Then Projects view should be visible
    When user clicks on Files---Demo tree node inside browse tree
    Then Demo view should be visible
    And Tutorials view should be visible
    And Projects view should be absent
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The pin on a preview tab keeps the view through the next single click
    Given simple mode is off
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user clicks on Files---Demo---demog.csv tree node inside browse tree
    Then the "demog" view should be current
    And "This is Browse preview. Click to keep it open" icon should be visible
    When user clicks on "This is Browse preview. Click to keep it open" icon
    Then "This is Browse preview. Click to keep it open" icon should be absent
    # the pinned view is a table view now, and its Toolbox takes the Browse tab's place
    Given the toolbox pane is hidden
    And the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    Then Projects view should be visible
    And demog view should be present
    And no errors should have been logged
    And no error or warning balloon should have been shown
