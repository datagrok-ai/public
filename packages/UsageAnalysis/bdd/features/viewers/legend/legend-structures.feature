@journey @viewers @realizes:viewers.legend
Feature: Molecules in the legend
  A legend whose column holds molecules draws each category through the molecule renderer instead
  of printing its SMILES, on every viewer that takes the column as its legend source; a scatter
  plot whose Markers column holds molecules gets structure items for the markers next to plain
  items for its colors; and the structures come back from a saved layout and a saved project.
  One journey on spgi-100, the one dataset of the legend features that is not demog-1000: its Core
  column is detected as Molecule and holds 7 cores, Series holds 5 series and Id 100 identifiers.
  Each viewer takes Core through the property that feeds its legend, as in the legend-across-viewers
  feature (scatter plot Color, histogram and line chart Split, bar chart Stack, pie chart Category,
  trellis plot X with an inner scatter plot colored by Core, box plot Category and Marker Color),
  each legend docked on the right with Visibility Always. "Drawn as a structure" is read from every
  item, each list scrolled through (the legend renders only the rows in view): the renderer's
  canvas in place of the label, its paint spanning a figure rather than one line of text. The
  layout scenario turns every legend to Series first, so each one has to come back to structures.
  Not translated: that the scatter plot's markers on the canvas are drawn as molecules — the scatter
  plot reports no reading of the marker it draws, and a pixel count would not tell a molecule from
  a large marker; needs a marker reading in the scatter plot's status — for the same reason the
  project scenario claims the legends, not the markers. The structure filter of the filtering case
  belongs to the Chem package.

  Background:
    Given user is logged in
    And user opens spgi dataset
    Then "Core" column should have semantic type "Molecule"
    When user adds a scatter plot viewer with:
      | colorColumnName   | Core   |
      | Legend Visibility | Always |
      | Legend Position   | Right  |
    And user adds a histogram viewer with:
      | splitColumnName   | Core   |
      | Legend Visibility | Always |
      | Legend Position   | Right  |
    And user adds a line chart viewer with:
      | splitColumnNames  | Core   |
      | Legend Visibility | Always |
      | Legend Position   | Right  |
    And user adds a bar chart viewer with:
      | splitColumnName   | Series |
      | stackColumnName   | Core   |
      | Legend Visibility | Always |
      | Legend Position   | Right  |
    And user adds a pie chart viewer with:
      | categoryColumnName | Core   |
      | Legend Visibility  | Always |
      | Legend Position    | Right  |
    And user adds a trellis plot viewer with:
      | xColumnNames      | Core         |
      | Viewer Type       | Scatter plot |
      | Legend Visibility | Always       |
      | Legend Position   | Right        |
    And user sets "colorColumnName" inner property of trellis plot viewer to "Core"
    And user adds a box plot viewer with:
      | categoryColumnNames   | Core   |
      | markerColorColumnName | Core   |
      | Legend Visibility     | Always |
      | Legend Position       | Right  |

  Scenario Outline: The <viewer> legend draws the seven cores as structures
    Then the legend of <viewer> viewer should list 7 items
    And every item in the legend of <viewer> viewer should be drawn as a structure
    And no errors should have been logged

    Examples:
      | viewer       |
      | scatter plot |
      | histogram    |
      | line chart   |
      | bar chart    |
      | pie chart    |
      | trellis plot |
      | box plot     |

  Scenario: Markers on Core add structure items next to the color items of Series
    When user sets "Color" property of scatter plot viewer to "Series"
    Then the legend of scatter plot viewer should list 5 items
    And every item in the legend of scatter plot viewer should be drawn as text
    When user sets "Markers" property of scatter plot viewer to "Core"
    Then the legend of scatter plot viewer should list 12 items
    And thumbnail of last legend item in legend of scatter plot viewer should be visible
    And thumbnail of first legend item in legend of scatter plot viewer should be absent
    When user sets "Color" property of scatter plot viewer to "Id"
    Then the legend of scatter plot viewer should list 107 items
    And thumbnail of last legend item in legend of scatter plot viewer should be visible
    And thumbnail of first legend item in legend of scatter plot viewer should be absent
    When user sets properties of scatter plot viewer:
      | Color   | Core |
      | Markers |      |
    Then the legend of scatter plot viewer should list 7 items
    And every item in the legend of scatter plot viewer should be drawn as a structure
    And no errors should have been logged

  Scenario: The structures come back from a saved layout
    When user saves the layout of the current table view to the server
    And user sets "Color" property of scatter plot viewer to "Series"
    And user sets "splitColumnName" property of histogram viewer to "Series"
    And user sets "splitColumnNames" property of line chart viewer to "Series"
    And user sets "stackColumnName" property of bar chart viewer to "Series"
    And user sets "Category" property of pie chart viewer to "Series"
    And user sets "colorColumnName" inner property of trellis plot viewer to "Series"
    And user sets "markerColorColumnName" property of box plot viewer to "Series"
    Then every item in the legend of scatter plot viewer should be drawn as text
    And every item in the legend of histogram viewer should be drawn as text
    And every item in the legend of line chart viewer should be drawn as text
    And every item in the legend of bar chart viewer should be drawn as text
    And every item in the legend of pie chart viewer should be drawn as text
    And every item in the legend of trellis plot viewer should be drawn as text
    And every item in the legend of box plot viewer should be drawn as text
    When user loads the saved layout
    Then "Category" property of pie chart viewer should be "Core"
    And the legend of pie chart viewer should list 7 items
    And every item in the legend of pie chart viewer should be drawn as a structure
    And every item in the legend of scatter plot viewer should be drawn as a structure
    And every item in the legend of histogram viewer should be drawn as a structure
    And every item in the legend of line chart viewer should be drawn as a structure
    And every item in the legend of bar chart viewer should be drawn as a structure
    And every item in the legend of trellis plot viewer should be drawn as a structure
    And every item in the legend of box plot viewer should be drawn as a structure
    And no errors should have been logged

  Scenario: The structures come back from a saved project
    When user saves the current view as project "bdd-legend-structures"
    And user closes all views
    And user opens the "bdd-legend-structures" project
    Then "Core" column should have semantic type "Molecule"
    And the open tableview should have 1 box plot viewer
    And the legend of scatter plot viewer should list 7 items
    And every item in the legend of scatter plot viewer should be drawn as a structure
    And every item in the legend of histogram viewer should be drawn as a structure
    And every item in the legend of line chart viewer should be drawn as a structure
    And every item in the legend of bar chart viewer should be drawn as a structure
    And every item in the legend of pie chart viewer should be drawn as a structure
    And every item in the legend of trellis plot viewer should be drawn as a structure
    And every item in the legend of box plot viewer should be drawn as a structure
    And no errors should have been logged
