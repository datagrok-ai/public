@browse @realizes:views.browse
Feature: Working with the nodes of the Browse tree
  Expanding and collapsing a node, driving the tree from the keyboard, and the expanded set
  surviving a trip to another view. Translated from the manual cases Browse-Tree-01, -02 and -03
  (playwright-public/browse/tree.test.ts, browse_manual_tests2.md section 2).

  Browse-Tree-04 (hiding and showing the panel must not leave a nested node stuck open, GROK-19802)
  is claimed on the node itself: a child closed by hand stays closed once the panel is back and the
  tree reports that no group is still fetching its children — the old spec allowed the expanded count
  to grow by one, which is the very defect.

  Browse-Tree-02 step 5 (Enter or Space opens the selected item) is not claimed: in the tree an arrow
  key already opens the file it moves to, and Enter on a selected item whose view was closed opens
  nothing — a candidate finding, walked by hand before anything is claimed or filed. Browse-Tree-05 (the context menu of an
  entity) is in browse-context-panel-and-menus.feature, on a connection and on a file — the dashboard
  the manual case names is not used. Browse-Tree-06 (drag and drop) is claimed in
  spaces-drag-and-drop.feature; Browse-Tree-07 (a node the user may not read) in
  browse-platform-and-databases.feature, through the second account.

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

  Scenario: Hiding and showing the panel leaves a closed child closed
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    And Files---Demo---chem tree node inside browse tree should be visible
    When user collapses Files---Demo tree node inside browse tree
    Then Files---Demo---chem tree node inside browse tree should be hidden
    When user clicks on browse tab
    Then the browse tree should be hidden
    When user clicks on browse tab
    Then the browse tree should be visible
    And the browse tree should have finished loading
    And Files tree node inside browse tree should be expanded
    And Files---Demo tree node inside browse tree should be collapsed
    And Files---Demo---chem tree node inside browse tree should be hidden
    And Files---App-Data tree node inside browse tree should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: An arrow key opens the file it moves to
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user clicks on Files---Demo---demog-1000.csv tree node inside browse tree
    Then the "demog-1000" view should be current
    When user presses ArrowDown
    Then Files---Demo---demog.csv tree node inside browse tree should be selected
    And the "demog" view should be current
    And grid should show 5850 rows
    And no errors should have been logged
