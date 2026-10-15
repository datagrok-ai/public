@journey @viewers @realizes:viewers.tooltips
Feature: The Tooltip group of a viewer's context menu
  Every viewer of a table view, the view's own grid included, carries a Tooltip group in its context
  menu that includes the same four actions: Hide, Edit..., Use as Group Tooltip and Remove Group Tooltip.
  Translated from the TestTrack case Tooltips/actions-in-the-context-menu (on demog-1000) with a histogram, a
  line chart, a bar chart and a trellis plot, and the grid.

  What the actions do is the subject of the other tooltip features (Hide and Show Custom in
  default-tooltip-visibility and tooltip-properties, Edit... in edit-tooltip); the group tooltip
  actions are walked by the Embedded Viewers tutorial (Tutorials/bdd).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a histogram viewer
    And user adds a line chart viewer
    And user adds a bar chart viewer
    And user adds a trellis plot viewer

  Scenario Outline: The context menu of the <viewer> offers the four Tooltip actions
    When user opens the context menu of <viewer>
    Then the open menu should list "Tooltip > Hide"
    And the open menu should list "Tooltip > Edit..."
    And the open menu should list "Tooltip > Use as Group Tooltip"
    And the open menu should list "Tooltip > Remove Group Tooltip"
    When user closes the context menu
    Then no errors should have been logged

    Examples:
      | viewer              |
      | grid                |
      | histogram viewer    |
      | line chart viewer   |
      | bar chart viewer    |
      | trellis plot viewer |
