@journey @viewers @realizes:viewers.bar-chart
Feature: Bar chart property surface
  The bar chart's property surface: coloring by a value column and its aggregation, the
  missing-values bar (Include Nulls), bar style, labels, the on-chart selectors and the axes
  (Controls), the value aggregation and the value column, the stack legend and its position, title
  and description, values in place of category names, the context menu as a path to properties,
  and the Data panel — row source, table, filter and color column through a layout round-trip.
  Selection overlays, sorting and orientation, stacking and the value axis have features of their
  own. One journey on demog-1000 with a bar chart of average AGE by RACE; every scenario puts back
  what it changed.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a bar chart viewer with:
      | Split           | RACE |
      | Value           | AGE  |
      | Value Aggr Type | avg  |
    Then the "bars" reading of bar chart viewer should be 4
    And bar chart viewer should show 1000 rows

  Scenario: Color coding by a value column
    Then "Color" property of bar chart viewer should be ""
    And the "bar Caucasian" area of bar chart viewer should contain the color "#96D794"
    When user sets "Color" property of bar chart viewer to "HEIGHT"
    Then "Color" property of bar chart viewer should be "HEIGHT"
    And bar chart viewer should have repainted by at least 2000 pixels
    When user sets "Color Aggr Type" property of bar chart viewer to "min"
    Then bar chart viewer should have repainted
    When user sets "Color Aggr Type" property of bar chart viewer to "max"
    Then bar chart viewer should have repainted
    When user sets "Color Aggr Type" property of bar chart viewer to "med"
    Then "Color Aggr Type" property of bar chart viewer should be "med"
    When user sets "Invert Color Scheme" property of bar chart viewer to "true"
    Then bar chart viewer should have repainted by at least 2000 pixels
    When user sets properties of bar chart viewer:
      | Color               |       |
      | Color Aggr Type     | avg   |
      | Invert Color Scheme | false |
    Then "Color" property of bar chart viewer should be ""
    And the "bar Caucasian" area of bar chart viewer should contain the color "#96D794"
    And no errors should have been logged

  Scenario: Include Nulls gates the missing-values bar
    Given table "demog-1000" should have missing values in "HEIGHT" column
    When user sets "Split" property of bar chart viewer to "HEIGHT"
    Then "Include Nulls" property of bar chart viewer should be "true"
    And the "bars" reading of bar chart viewer should be higher than before
    And bar chart viewer should show 1000 rows
    When user sets "Include Nulls" property of bar chart viewer to "false"
    Then bar chart viewer should show 872 rows
    And bar chart viewer should have repainted
    When user sets "Include Nulls" property of bar chart viewer to "true"
    Then bar chart viewer should show 1000 rows
    And bar chart viewer should have repainted
    When user sets "Split" property of bar chart viewer to "RACE"
    Then the "bars" reading of bar chart viewer should be 4
    And no errors should have been logged

  Scenario: Bar style
    When user sets "Bar Border Line Width" property of bar chart viewer to "2"
    Then bar chart viewer should have repainted by at least 300 pixels
    When user sets "Max Bar Height" property of bar chart viewer to "20"
    Then the "bar Caucasian" area of bar chart viewer should have less ink than before
    When user sets "Bar Corner Radius" property of bar chart viewer to "10"
    Then bar chart viewer should have repainted
    When user sets "Vertical Align" property of bar chart viewer to "Top"
    Then "Vertical Align" property of bar chart viewer should be "Top"
    When user sets "Vertical Align" property of bar chart viewer to "Bottom"
    Then "Vertical Align" property of bar chart viewer should be "Bottom"
    When user sets "Vertical Align" property of bar chart viewer to "Center"
    Then "Vertical Align" property of bar chart viewer should be "Center"
    When user sets "Show Category Zero Baseline" property of bar chart viewer to "false"
    Then "Show Category Zero Baseline" property of bar chart viewer should be "false"
    When user sets properties of bar chart viewer:
      | Bar Border Line Width       | 0    |
      | Bar Corner Radius           | 0    |
      | Max Bar Height              | 50   |
      | Show Category Zero Baseline | true |
    Then bar chart viewer should have repainted
    And no errors should have been logged

  Scenario: Labels
    When user sets "Show Labels" property of bar chart viewer to "inside"
    Then "Show Labels" property of bar chart viewer should be "inside"
    When user sets "Show Labels" property of bar chart viewer to "never"
    Then bar chart viewer should have repainted by at least 300 pixels
    When user sets "Show Labels" property of bar chart viewer to "outside"
    Then bar chart viewer should have repainted by at least 300 pixels
    When user sets "Show Labels" property of bar chart viewer to "auto"
    Then "Show Labels" property of bar chart viewer should be "auto"
    And no errors should have been logged

  Scenario: Controls visibility
    When user hovers over bar chart viewer
    Then Value column input in bar chart viewer should be visible
    And Split column input in bar chart viewer should be visible
    And Stack column input in bar chart viewer should be visible
    When user sets properties of bar chart viewer:
      | Show Value Selector    | false |
      | Show Category Selector | false |
      | Show Stack Selector    | false |
    And user hovers over bar chart viewer
    Then Value column input in bar chart viewer should be hidden
    And Split column input in bar chart viewer should be hidden
    And Stack column input in bar chart viewer should be hidden
    When user sets "Show Value Axis" property of bar chart viewer to "false"
    Then bar chart viewer should not have an "x axis" area
    And bar chart viewer should have repainted by at least 500 pixels
    When user sets "Show Category Values" property of bar chart viewer to "false"
    Then bar chart viewer should not have a "y axis" area
    And bar chart viewer should have repainted by at least 500 pixels
    When user sets properties of bar chart viewer:
      | Show Value Selector    | true |
      | Show Category Selector | true |
      | Show Stack Selector    | true |
      | Show Value Axis        | true |
      | Show Category Values   | true |
    Then bar chart viewer should have an "x axis" area
    And bar chart viewer should have a "y axis" area
    When user hovers over bar chart viewer
    Then Value column input in bar chart viewer should be visible
    And Split column input in bar chart viewer should be visible
    And Stack column input in bar chart viewer should be visible
    And no errors should have been logged

  Scenario: Aggregation and the value column
    When user sets "Value Aggr Type" property of bar chart viewer to "max"
    Then bar chart viewer should have repainted by at least 500 pixels
    And the bars of bar chart viewer should differ in length
    When user sets "Value" property of bar chart viewer to "WEIGHT"
    Then bar chart viewer should have repainted by at least 500 pixels
    When user sets properties of bar chart viewer:
      | Value           | AGE |
      | Value Aggr Type | avg |
    Then no errors should have been logged

  Scenario: The stack legend and its position
    Then legend of bar chart viewer should be hidden
    When user sets properties of bar chart viewer:
      | Value Aggr Type   | count  |
      | Stack             | SEX    |
      | Legend Visibility | Always |
    Then legend of bar chart viewer should be visible
    And legend of bar chart viewer should have 2 items
    And the "stack segments" reading of bar chart viewer should be higher than before
    When user sets "Legend Position" property of bar chart viewer to "Left"
    Then the legend of bar chart viewer should be on the left
    When user sets "Legend Position" property of bar chart viewer to "Right"
    Then the legend of bar chart viewer should be on the right
    When user sets "Legend Position" property of bar chart viewer to "Top"
    Then the legend of bar chart viewer should be on the top
    When user sets "Legend Position" property of bar chart viewer to "Bottom"
    Then "Legend Position" property of bar chart viewer should be "Bottom"
    And legend of bar chart viewer should be visible
    When user sets properties of bar chart viewer:
      | Stack             |      |
      | Legend Position   | Auto |
      | Legend Visibility | Auto |
      | Value Aggr Type   | avg  |
    Then legend of bar chart viewer should be hidden
    And the "stack segments" reading of bar chart viewer should be 0
    And no errors should have been logged

  Scenario: Title and description
    When user sets properties of bar chart viewer:
      | Show Title | true         |
      | Title      | Demographics |
    Then title of bar chart viewer should have text "Demographics"
    When user sets properties of bar chart viewer:
      | Description                 | By race |
      | Description Visibility Mode | Always  |
    Then description of bar chart viewer should have text "By race"
    When user sets "Description Position" property of bar chart viewer to "Bottom"
    Then description of bar chart viewer should be visible
    When user sets "Description Position" property of bar chart viewer to "Left"
    Then description of bar chart viewer should be visible
    When user sets "Description Position" property of bar chart viewer to "Right"
    Then description of bar chart viewer should be visible
    When user sets "Description Visibility Mode" property of bar chart viewer to "Never"
    Then description of bar chart viewer should be absent
    And no errors should have been logged
    When user sets properties of bar chart viewer:
      | Show Title                  | false |
      | Title                       |       |
      | Description                 |       |
      | Description Visibility Mode | Auto  |
      | Description Position        | Top   |

  Scenario: Values in place of category names
    When user sets "Show Values Instead Of Categories" property of bar chart viewer to "true"
    Then "Show Values Instead Of Categories" property of bar chart viewer should be "true"
    And bar chart viewer should have repainted by at least 300 pixels
    When user sets "Show Values Instead Of Categories" property of bar chart viewer to "false"
    Then bar chart viewer should have repainted by at least 300 pixels
    And no errors should have been logged

  Scenario: The context menu as a path to properties
    When user opens the context menu of bar chart viewer
    Then "Reset View" menu item in context menu should be visible
    And "Orientation" menu item in context menu should be visible
    And "Data" menu item in context menu should be visible
    And "Order" menu item in context menu should be visible
    And "Controls" menu item in context menu should be visible
    And "Selection" menu item in context menu should be visible
    When user hovers over "Controls" menu item in context menu
    Then "Show Value Axis" menu item in context menu should be visible
    And "Show Category Values" menu item in context menu should be visible
    When user hovers over "Data" menu item in context menu
    Then "Relative Values" menu item in context menu should be visible
    And "Include Nulls" menu item in context menu should be visible
    When user hovers over "Selection" menu item in context menu
    Then "Show Selected Rows" menu item in context menu should be visible
    When user closes the context menu
    And user picks "Controls > Show Value Axis" from the context menu of bar chart viewer
    Then "Show Value Axis" property of bar chart viewer should be "false"
    And bar chart viewer should not have an "x axis" area
    And bar chart viewer should have repainted by at least 500 pixels
    When user picks "Controls > Show Value Axis" from the context menu of bar chart viewer
    Then "Show Value Axis" property of bar chart viewer should be "true"
    And bar chart viewer should have an "x axis" area
    And bar chart viewer should have repainted by at least 500 pixels
    When user picks "Data > Include Nulls" from the context menu of bar chart viewer
    Then "Include Nulls" property of bar chart viewer should be "false"
    When user picks "Data > Include Nulls" from the context menu of bar chart viewer
    Then "Include Nulls" property of bar chart viewer should be "true"
    And no errors should have been logged

  Scenario: The Data panel through a layout round-trip on the server
    Given user opens spgi dataset
    And user switches to the "demog-1000" table view
    When user sets "Row Source" property of bar chart viewer to "Filtered"
    Then "Row Source" property of bar chart viewer should be "Filtered"
    And bar chart viewer should show 1000 rows
    When user sets "Row Source" property of bar chart viewer to "Selected"
    Then bar chart viewer should show 0 rows
    When user sets "Row Source" property of bar chart viewer to "All"
    Then bar chart viewer should show 1000 rows
    When user sets "Table" property of bar chart viewer to "spgi-100"
    Then bar chart viewer should be bound to table "spgi-100"
    When user sets properties of bar chart viewer:
      | Split           | Primary Series Name |
      | Value           | CAST Idea ID        |
      | Value Aggr Type | count               |
    Then the "bars" reading of bar chart viewer should be 5
    And bar chart viewer should show 100 rows
    When user sets "Filter" property of bar chart viewer to "${CAST Idea ID} < 634835"
    Then bar chart viewer should show 50 rows
    When user sets "Color" property of bar chart viewer to "Chemical Space Y"
    Then bar chart viewer should have repainted by at least 2000 pixels
    When user saves the layout of the current table view to the server
    And user clicks on close icon of bar chart viewer
    Then bar chart viewer should be absent
    When user loads the saved layout
    Then bar chart viewer should be visible
    And properties of bar chart viewer should be:
      | Table  | spgi-100                 |
      | Filter | ${CAST Idea ID} < 634835 |
      | Color  | Chemical Space Y         |
    And bar chart viewer should show 50 rows
    And no errors should have been logged
    When user sets properties of bar chart viewer:
      | Filter |            |
      | Color  |            |
      | Table  | demog-1000 |
    And user sets properties of bar chart viewer:
      | Split           | RACE |
      | Value           | AGE  |
      | Value Aggr Type | avg  |
    Then bar chart viewer should be bound to table "demog-1000"
    And bar chart viewer should show 1000 rows
