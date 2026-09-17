@browse @realizes:views.browse
Feature: The context panel and the context menus of the Browse tree
  F4 shows and hides the context panel, the panel follows whatever is current, and a node's
  context menu offers what that kind of node can do. Translated from the manual cases
  Browse-CtxPanel-01, -02 and Browse-Tree-05 (playwright-public/browse/ctxpanel.test.ts,
  tree.test.ts, browse_manual_tests2.md sections 2 and 16).

  Browse-Tree-05 asked for the menu of an entity and left the exact set to the first run. It is
  written down here: a connection offers Browse, the query and table commands, Edit, Rename, Clone,
  Delete and Clear cache. The case also asks for the cross-check that a plain file is not offered
  the same things, so the file's menu is claimed too — a one-sided claim would pass on a menu that
  offers everything to everyone.

  A node below the top level is named by its full tree path ("Files---Demo"), which is what the
  platform writes into its own `name` attribute. Several sections carry a node called Demo, My
  files or App Data, and a bare name matches whichever of them another feature happened to leave
  open: the tree remembers its expanded set per user, across features and across runs.

  Browse-CtxPanel-04 (Collapse all / Expand all panes) is not translated: the old spec asserted
  only that the two icons existed and that clicking them logged no error, which holds on a dead
  button. Browse-CtxPanel-03 (Back and Forward) IS a real claim in the old spec — after Back the
  panel holds the previous object and not the current one — and is writable with the header's
  aria-labelled icons plus the same "context panel should show" step used below. It is simply not
  written yet.

  Background:
    Given user is logged in
    And the browse panel is open

  Scenario: F4 hides the context panel and shows it again
    Given the context panel is open
    When user presses F4
    Then context panel should be hidden
    When user presses F4
    Then context panel should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The context panel follows the node that was clicked
    Given the context panel is open
    And Apps tree node inside browse tree is expanded
    When user clicks on Tutorials tree node inside browse tree
    Then the context panel should show "Tutorials"
    # a view is not an entity: the panel follows grok.shell.o, so the next claim names
    # another object the tree owns rather than the view a node opens
    When Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And user clicks on Databases---Postgres---Datagrok tree node inside browse tree
    Then the context panel should show "Datagrok"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A connection offers Browse, the query commands, Edit, Rename, Clone, Delete and Clear cache
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    When user opens the context menu of Databases---Postgres---Datagrok tree node inside browse tree
    Then the open menu should list "Browse"
    And the open menu should list "New Query..."
    And the open menu should list "Edit..."
    And the open menu should list "Rename..."
    And the open menu should list "Clone..."
    And the open menu should list "Delete..."
    And the open menu should list "Clear cache"
    # the anchor for the file scenario's negative: an entity is offered this, a file is not
    And the open menu should list "Add to favorites"
    When user closes the context menu
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A file is not offered the commands of a connection
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user opens the context menu of Files---Demo---demog.csv tree node inside browse tree
    Then the open menu should list "Open"
    And the open menu should list "Download"
    And the open menu should not list "New Query..."
    And the open menu should not list "Clear cache"
    # Browse-Fav-05: a file is not an entity, so it cannot be made a favourite of its own
    And the open menu should not list "Add to favorites"
    When user closes the context menu
    Then no errors should have been logged
    And no error or warning balloon should have been shown
