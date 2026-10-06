@demo @realizes:u2.accordion @realizes:u2.breadcrumbs @realizes:u2.toolbar @realizes:u2.splitter
Feature: Layout containers
  Accordion panes, breadcrumbs, a panel-local toolbar and a splitter, each with the readout that
  mirrors its signal.

  Background:
    Given user opens the "Layout" demo page

  Scenario: Accordion panes expand and collapse
    Then General pane should be expanded
    And Advanced pane should be collapsed
    When user expands Advanced pane
    Then Advanced pane should be expanded
    And Advanced pane should contain text "Lazily built on first expand"
    When user collapses Advanced pane
    Then Advanced pane should be collapsed

  Scenario: Breadcrumbs truncate the path
    Then value of path readout should have text "Home > Projects > Demo > Tables > demog"
    When user clicks on Projects breadcrumb
    Then value of path readout should have text "Home > Projects"
    And breadcrumbs should not contain text "demog"

  Scenario: A panel-local toolbar with buttons, a toggle and a menu
    When user clicks on Run button in toolbar
    Then value of toolbar readout should contain text "Run clicked"
    When user clicks on Save button in toolbar
    Then value of toolbar readout should contain text "Saved"
    When user clicks on Compact button in toolbar
    Then Compact button in toolbar should be selected
    And value of toolbar readout should contain text "compact true"
    When user clicks on Help button in toolbar
    And user clicks on "About u2" menu item
    Then value of toolbar readout should contain text "u2: next-gen Datagrok UI library"

  Scenario: The splitter reports its sizes
    Then value of sizes readout should have text "0.3 / 0.7"
    When user drags sash in splitter to Content heading
    Then value of sizes readout should not contain text "0.3 / 0.7"
