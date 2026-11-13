@journey @viewers @realizes:viewers.matrix-plot
Feature: Matrix plot — a filter-panel filter in every cell, and the Show X Axes checkbox a person clicks
  Two anchors of the old specs the other matrix plot features leave open. "Filter Panel narrows the
  shared filter" (matrixplot-row-source-filter.md): a SEX = M card in the filter panel narrows every
  cell, each to its own rows, and closing the card puts every cell back as it was — claimed on the
  cells' own `cell rows shown` and on the cell's picture, where the old spec checked only the table's
  filter count. And "Show X Axes checkbox under Auto Layout" (GROK-19106): the real checkbox in the
  context panel, not the property set through the API that `matrix-plot-layout.feature` uses, takes
  the axis strip away and brings it back.
  Not translated: hovering, clicking or dragging over the points inside a cell
  (`matrix-plot-ui.md`, "Interact with Plot Elements") — see MISSING.md.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a matrix plot viewer
    Then 1000 rows should pass the filter
    And the "rows shown" reading of matrix plot viewer should be 1000
    And the "cells" reading of matrix plot viewer should be 16

  Scenario: A SEX = M card narrows every cell to its own rows, and closing it restores them
    # The panel is opened first: docking it narrows the plot, and a signature taken before that
    # would differ after the card is closed for the size alone.
    When user opens an empty filter panel
    Then the "distinct cell signatures" reading of matrix plot viewer should be 16
    And the "cell rows shown of HEIGHT x AGE" reading of matrix plot viewer should be 872
    When user remembers the "cell signature HEIGHT x AGE" reading of matrix plot viewer
    And user adds a categorical filter on "SEX" keeping "M"
    Then 447 rows should pass the filter
    And the "rows shown" reading of matrix plot viewer should be 447
    And the "cell rows shown of AGE x AGE" reading of matrix plot viewer should be 447
    And the "cell rows shown of HEIGHT x AGE" reading of matrix plot viewer should be 373
    And the "cell signature HEIGHT x AGE" reading of matrix plot viewer should not be as remembered
    When user hovers over "SEX" filter card
    And user clicks on close of "SEX" filter card
    Then 1000 rows should pass the filter
    And the "rows shown" reading of matrix plot viewer should be 1000
    And the "cell rows shown of HEIGHT x AGE" reading of matrix plot viewer should be 872
    And the "cell signature HEIGHT x AGE" reading of matrix plot viewer should be as remembered
    And no errors should have been logged

  Scenario: The Show X Axes checkbox in the context panel takes the axis strip away and brings it back (GROK-19106)
    Then the "auto layout" reading of matrix plot viewer should be "true"
    And the "x axis shown" reading of matrix plot viewer should be "true"
    And matrix plot viewer should have an "x axis" area
    When user clicks on settings icon of matrix plot viewer
    Given "Axes" category in context panel is expanded
    When user unchecks "Show X Axes" property in context panel
    Then "showXAxes" property of matrix plot viewer should be "false"
    And the "x axis shown" reading of matrix plot viewer should be "false"
    And matrix plot viewer should not have an "x axis" area
    And the "y axis shown" reading of matrix plot viewer should be "true"
    When user checks "Show X Axes" property in context panel
    Then the "x axis shown" reading of matrix plot viewer should be "true"
    And matrix plot viewer should have an "x axis" area
    And no errors should have been logged
