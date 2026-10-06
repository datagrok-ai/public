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
  platform writes into its own `name` attribute. Several sections carry a node called Demo, Files
  or App Data, and a bare name matches whichever of them another feature happened to leave
  open: the tree remembers its expanded set per user, across features and across runs.

  Browse-CtxPanel-03 (Back and Forward) is claimed on what the panel shows: Back renders the previous
  object with the current object cleared (property_panel.dart sets it to null without notice), so
  "the context panel should show" cannot back it, and the claim is the panel's own text — the object
  that must be back and the one that must be gone. Browse-CtxPanel-04 (Collapse all / Expand all) is
  claimed on a pane's own expanded state. The title bar's icons are reached by their labels: the help
  panel carries a Back and a Forward too, hidden while it is closed, and the icons show only while
  the pointer is on the bar. Browse-Fav-02 (the star beside the object's name toggles it in and out
  of favorites) is claimed on the account's favorites on the server — the star publishes its state
  only as a font class. Since GROK-21108 a file is a favorite too: Browse-Fav-05 and -05b, which
  denied a file the menu item and the star, are claimed the other way round. For an account that
  administers a group, as the running one does, "Add To Favorites" is a submenu of "Only for me"
  and the groups. The title bar's own "Favorites" icon is not a toggle: it lists the favorites.

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
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user clicks on Files---Demo---demog.csv tree node inside browse tree
    Then the context panel should show "demog.csv"
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
    And the open menu should list "Add To Favorites > Only for me"
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
    # Browse-Fav-05: a file is a favorite of its own since GROK-21108
    And the open menu should list "Add To Favorites > Only for me"
    When user closes the context menu
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Back and Forward walk the panel's history
    Given the context panel is open
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user clicks on Files---Demo---demog.csv tree node inside browse tree
    Then the context panel should show "demog.csv"
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    When user clicks on Databases---Postgres---Datagrok tree node inside browse tree
    Then the context panel should show "Datagrok"
    # the title bar shows its history icons only while the pointer is on it
    When user hovers over "Expand all" icon
    And user clicks on "Back" icon
    Then context panel should contain text "demog.csv"
    And context panel should not contain text "Datagrok"
    When user hovers over "Expand all" icon
    And user clicks on "Forward" icon
    Then context panel should contain text "Datagrok"
    And context panel should not contain text "demog.csv"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Collapse all and Expand all fold every pane of the panel
    Given the context panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    When user clicks on Databases---Postgres---Datagrok tree node inside browse tree
    Then the context panel should show "Datagrok"
    When user clicks on "Expand all" icon
    Then "Details" accordion header in context panel should be expanded
    When user clicks on "Collapse all" icon
    Then "Details" accordion header in context panel should be collapsed
    When user clicks on "Expand all" icon
    Then "Details" accordion header in context panel should be expanded
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The star beside a connection's name adds it to favorites and takes it out
    Given a "Postgres" connection named "BDD-Browse-Star-{run}" is on the server
    And "BDD-Browse-Star-{run}" is not in favorites
    And the context panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    When user clicks on Databases---Postgres---BDD-Browse-Star-{run} tree node inside browse tree
    Then the context panel should show "BDD-Browse-Star-{run}"
    When user clicks on favorite star in context panel
    Then "BDD-Browse-Star-{run}" should be in favorites on the server
    When user clicks on favorite star in context panel
    Then "BDD-Browse-Star-{run}" should not be in favorites on the server
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A file has a favorite star, as a connection has
    Given the context panel is open
    And Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    When user clicks on Databases---Postgres---Datagrok tree node inside browse tree
    Then the context panel should show "Datagrok"
    And favorite star in context panel should be visible
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user clicks on Files---Demo---demog.csv tree node inside browse tree
    Then the context panel should show "demog.csv"
    And favorite star in context panel should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown
