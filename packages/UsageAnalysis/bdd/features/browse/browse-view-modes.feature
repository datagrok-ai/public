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

  Browse-View-02, -03 and -05 are not translated, and there is no other file holding them.
  -03 (the pin control) is writable today — the preview tab's pin carries an aria-label — and is
  the first thing to add to this file. -02 needs an edit gesture on the previewed view, and -05
  needs the sidebar's count badge, which carries no name at all.

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
    # The third node is a folder rather than another top-level group: a group that has been
    # expanded far enough to paginate grows a "Show more" item that carries the group's own
    # name (tree_view.dart addMoreLink calls addItem('')), so "Databases tree node" matches
    # two elements once another feature has opened it.
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
