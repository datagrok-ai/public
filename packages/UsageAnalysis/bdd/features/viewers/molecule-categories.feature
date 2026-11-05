@journey @viewers
Feature: Molecules as the categories of five viewers
  A molecule column can be the category axis of the bar chart, the box plot, the pie chart, the
  trellis plot and the scatter plot. On spgi-100 Core holds 7 cores and an empty value: the bar
  chart splits into 7 bars, the box plot lays out a column per core, the pie chart cuts 7 slices,
  the trellis plot labels 7 columns and the scatter plot spreads its x axis over the 7 categories,
  each one named by its core's SMILES. Structure holds 100 different molecules, more than a category
  axis shows at once: dragging the bar chart's category slider shows fewer bars, and dragging the
  scatter plot's x slider moves the start of its axis, with the categories still laid out and no
  error (the regression of the old manual case).

  That the labels are drawn as structures, not as SMILES text, is Chem's renderer at work: the
  viewers report no reading of the renderer that draws a category label, so the drawing is the
  renderer's own package tests and the manual checklist. The pie chart has no axis slider. Needs
  Chem on the stand: without it the columns are plain strings.

  Background:
    Given user is logged in
    And the "Chem" package is installed
    And user opens spgi dataset
    Then "Core" column should have semantic type "Molecule"
    And "Structure" column should have semantic type "Molecule"

  Scenario: The bar chart splits into a bar per core
    When user adds a bar chart viewer with:
      | splitColumnName | Core |
    Then the "bars" reading of bar chart viewer should be 7
    And bar chart viewer should have a "bar n1nc([*:3])n([*:2])c1[*:1]" area
    And bar chart viewer should have a "y axis" area
    And no errors should have been logged

  Scenario: Dragging the bar chart's category slider over the molecules shows fewer bars
    When user sets "splitColumnName" property of bar chart viewer to "Structure"
    Then the "bars" reading of bar chart viewer should be 19
    When user drags the max handle of the "y" range slider of bar chart viewer by -120 pixels
    Then the "bars" reading of bar chart viewer should be lower than before
    And the "bars" reading of bar chart viewer should be at least 1
    And bar chart viewer should have a "y axis" area
    And no errors should have been logged
    When user clicks on close icon of bar chart viewer
    Then bar chart viewer should be absent

  Scenario: The box plot lays out a column per core
    When user adds a box plot viewer with:
      | categoryColumnNames | Core             |
      | valueColumnName     | Chemical Space X |
    Then box plot viewer should have a "category n1nc([*:3])n([*:2])c1[*:1]" area
    And box plot viewer should have a "category CC1(N[*:1])CCN([*:2])C1" area
    And no errors should have been logged
    When user clicks on close icon of box plot viewer
    Then box plot viewer should be absent

  Scenario: The pie chart cuts a slice per core
    When user adds a pie chart viewer with:
      | categoryColumnName | Core |
    Then the "slices" reading of pie chart viewer should be 7
    And pie chart viewer should have a 'slice "n1nc([*:3])n([*:2])c1[*:1]"' area
    And no errors should have been logged
    When user clicks on close icon of pie chart viewer
    Then pie chart viewer should be absent

  Scenario: The trellis plot labels a column per core
    When user adds a trellis plot viewer with:
      | xColumnNames | Core |
    Then the "x categories" reading of trellis plot viewer should be 7
    And the "x labels shown" reading of trellis plot viewer should be 7
    And no errors should have been logged
    When user clicks on close icon of trellis plot viewer
    Then trellis plot viewer should be absent

  Scenario: The scatter plot spreads its x axis over the cores, and its x slider moves over the molecules
    When user adds a scatter plot viewer with:
      | xColumnName | Core             |
      | yColumnName | Chemical Space Y |
    Then the "x axis span" reading of scatter plot viewer should be 6.6
    When user sets "xColumnName" property of scatter plot viewer to "Structure"
    Then the "x axis span" reading of scatter plot viewer should be at least 90
    When user drags the min handle of the "x" range slider of scatter plot viewer by 150 pixels
    Then the "x axis min" reading of scatter plot viewer should be higher than before
    And scatter plot viewer should have a "x axis" area
    And no errors should have been logged
