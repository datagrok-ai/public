@browse @realizes:views.browse
Feature: Every section of the Browse tree opens without an error
  The matrix of manual case Browse-Node-* (browse_manual_tests2.md section 18): click each node of
  the tree and require that nothing was logged. Translated from playwright-public/browse/
  section18.test.ts.

  The old spec walked a hard-coded list of providers guarded by `if (!(await node.isVisible()))
  continue;`. If the whole subtree failed to load, every node was invisible, the loop body never
  ran and the test was green. Here each row names its node, so a node that is not there fails the
  row instead of skipping it.

  The selection claim says the click reached the node, and no more: the tree sets that class
  inside its own click handler, before it opens anything or asks the server for anything. An error
  raised later in the same gesture lands in the next scenario's floor rather than this one's,
  because the floor is a snapshot rather than a poll — so these rows catch a node that throws on
  click, not one that fails to load afterwards.

  Background:
    Given user is logged in
    And the browse panel is open
    # open so that a pane which throws while rendering the clicked object lands in the floor
    And the context panel is open

  # "first": a section that has been expanded far enough to paginate grows a "Show more" item
  # carrying the section's own name (tree_view.dart addMoreLink calls addItem('')), so the bare
  # name matches two elements once another feature has opened it.
  Scenario Outline: Clicking <node> logs no error
    When user clicks on first <node> tree node inside browse tree
    Then first <node> tree node inside browse tree should be selected
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | node      |
      | My stuff  |
      | Spaces    |
      | Apps      |
      | Files     |
      | Dashboards|
      | Databases |
      | Platform  |

  # The tree remembers what was open, so each row closes its section before opening it: an
  # "expand" on an already-open node returns without touching the tree. The twistie alone would
  # not be evidence either — it flips on the click, before the children are asked for — so each
  # row names a child that has to arrive. Spaces and Dashboards are not here: the first has
  # nothing under it on a stand without spaces, the second is a leaf that opens a view.
  Scenario Outline: Opening the <section> section logs no error
    Given user collapses first <section> tree node inside browse tree
    And first <section> tree node inside browse tree should be collapsed
    When user expands first <section> tree node inside browse tree
    Then first <section> tree node inside browse tree should be expanded
    And <child> tree node inside browse tree should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | section   | child                     |
      | My stuff  | My-stuff---Recent         |
      | Apps      | Apps---Compute            |
      | Files     | Files---Demo              |
      | Databases | Databases---Postgres      |
      | Platform  | Platform---Users          |
