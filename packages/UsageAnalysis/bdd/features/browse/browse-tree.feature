@browse @realizes:views.browse
Feature: Working with the nodes of the Browse tree
  Expanding and collapsing a node, driving the tree from the keyboard, and the expanded set
  surviving a trip to another view. Translated from the manual cases Browse-Tree-01, -02 and -03
  (playwright-public/browse/tree.test.ts, browse_manual_tests2.md section 2).

  Browse-Tree-04 (the sidebar must not leave nested nodes stuck open, GROK-19802) is not here. The
  old spec allowed the expanded count to grow by one — which is exactly the defect the case is
  named for, so it passed with the bug present. A node-level claim is writable (collapse a child,
  toggle the sidebar, claim that child still collapsed); what it needs is a positive anchor that
  the rebuild publishes, and the tree publishes nothing when a group has finished loading.

  Browse-Tree-05 (the context menu of an entity) is in browse-context-panel-and-menus.feature,
  on a connection and on a file. The dashboard the manual case names is not used, so the entity
  menu is claimed on a different kind of entity than the case asks for. Browse-Tree-06 (drag and
  drop) and -07 (a node the user may not read) are not translated: the first asserted nothing
  about the drop, the second needs the second account.

  A node below the top level is named by its full tree path ("Files---Demo"), which is what the
  platform writes into its own `name` attribute. Several sections carry a node called Demo, Files
  or App Data, and a bare name matches whichever of them another feature happened to leave
  open: the tree remembers its expanded set per user, across features and across runs.

  Background:
    Given user is logged in
    And the browse panel is open

  Scenario: A node opens and closes on its own twistie
    # The tree remembers what was open, so the section is closed first: otherwise "is expanded"
    # returns without touching anything and only the closing half is ever driven by a gesture.
    Given user collapses Files tree node inside browse tree
    And Files---Demo tree node inside browse tree should be hidden
    When user expands Files tree node inside browse tree
    Then Files---Demo tree node inside browse tree should be visible
    And Files---App-Data tree node inside browse tree should be visible
    When user collapses Files tree node inside browse tree
    Then Files---Demo tree node inside browse tree should be hidden
    And Files---App-Data tree node inside browse tree should be hidden
    And Files tree node inside browse tree should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The arrow keys walk the tree and open a node
    # ArrowDown descends into an expanded group, so the pre-state has to be stated: with Files
    # left open by an earlier feature the next node is its first child, not Dashboards.
    Given user collapses Files tree node inside browse tree
    When user clicks on Files tree node inside browse tree
    Then Files tree node inside browse tree should be selected
    When user presses ArrowDown
    Then Dashboards tree node inside browse tree should be selected
    And Files tree node inside browse tree should not be selected
    When user presses ArrowUp
    Then Files tree node inside browse tree should be selected
    # Right opens the group and descends into its first child (whichever share the stand lists
    # first), so the first Left climbs back to the group and only the second one closes it
    # (tree_view.dart, the root key handler).
    When user presses ArrowRight
    Then Files---Demo tree node inside browse tree should be visible
    And Files tree node inside browse tree should not be selected
    When user presses ArrowLeft
    Then Files tree node inside browse tree should be selected
    And Files---Demo tree node inside browse tree should be visible
    When user presses ArrowLeft
    Then Files---Demo tree node inside browse tree should be hidden
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The tree keeps what it had open across a visit to another view
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree should be visible
    When user clicks on "Open text" icon inside browse toolbar
    Then the "Import text" view should be current
    # the panel has to actually go away and come back: opening another view leaves it standing,
    # so re-asserting here would claim a state nothing had disturbed
    When user clicks on browse tab
    Then the browse tree should be hidden
    When user clicks on browse tab
    Then the browse tree should be visible
    And Files---Demo tree node inside browse tree should be visible
    And Files---App-Data tree node inside browse tree should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown
