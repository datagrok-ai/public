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
    And the context panel is open

  Scenario Outline: Clicking <node> logs no error
    When user clicks on <node> tree node inside browse tree
    Then <node> tree node inside browse tree should be selected
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
  # "expand" on an already-open node returns without touching the tree.
  Scenario Outline: Opening the <section> section logs no error
    Given user collapses <section> tree node inside browse tree
    And <section> tree node inside browse tree should be collapsed
    When user expands <section> tree node inside browse tree
    Then <section> tree node inside browse tree should be expanded
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | section   |
      | My stuff  |
      | Apps      |
      | Files     |
      | Databases |
      | Platform  |
