@journey @viewers @realizes:viewers.box-plot
Feature: Box plot property surface
  The box plot's property surface: context menus as paths to properties, the statistics and
  group-comparison hit regions, auto layout and resize stability, the Show Markers gate, whisker
  and control-band style, controls visibility, title and description, axis font, date category
  mapping, the custom row tooltip, table switching, coloring, and the double-click view reset.
  Selection, persistence, filtering, statistics and group comparison have features of their own.
  One journey: demog-1000 (a stratified 1000-row demog: one marker per row per paint) with a box
  plot of AGE by SEX is opened once, and every scenario puts back what it changed.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a box plot viewer with:
      | Value      | AGE |
      | Category 1 | SEX |

  Scenario: Context menus as property paths
    Then "Show Inside Values" property of box plot viewer should be "true"
    When user picks "Misc > Show Inside Values" from the context menu of box plot viewer
    Then "Show Inside Values" property of box plot viewer should be "false"
    And box plot viewer should have less ink than before
    When user picks "Misc > Show Outside Values" from the context menu of box plot viewer
    Then "Show Outside Values" property of box plot viewer should be "false"
    And box plot viewer should have less ink than before
    When user picks "Misc > Show Inside Values" from the context menu of box plot viewer
    And user picks "Misc > Show Outside Values" from the context menu of box plot viewer
    Then "Show Inside Values" property of box plot viewer should be "true"
    And "Show Outside Values" property of box plot viewer should be "true"
    And box plot viewer should have more ink than before
    When user sets "Marker Size Column" property of box plot viewer to "WEIGHT"
    And user opens the context menu of box plot viewer
    And user hovers over Markers menu item in context menu
    Then "Markers > Size" menu item in context menu should be disabled
    When user closes the context menu
    And user sets "Marker Size Column" property of box plot viewer to ""

  Scenario: Statistics and group-comparison menu regions
    Given user sets properties of box plot viewer:
      | Show Statistics       | true  |
      | Show Group Comparison | false |
      | Show P Value          | true  |
    When user right-clicks on the "stats" area of box plot viewer
    And user hovers over "Group Comparison" menu item in context menu
    Then "Show Assumption Checks" menu item in context menu should be disabled
    When user hovers over "Show Assumption Checks" menu item in context menu
    Then tooltip should contain text "Show Group Comparison"
    When user closes the context menu
    And user right-clicks on the "p value" area of box plot viewer
    Then context menu should contain text "Show P Value"
    And context menu should not contain text "Statistics Format"
    When user closes the context menu
    And user moves the pointer away from box plot viewer
    And user hovers over the "p value" area of box plot viewer
    And user clicks on show-group-stats icon in box plot viewer
    Then "Show Group Comparison" property of box plot viewer should be "true"
    When user right-clicks on the "group comparison" area of box plot viewer
    Then "Add t-Test Table" menu item in context menu should be visible
    And "Show Assumption Checks" menu item in context menu should be enabled
    When user closes the context menu
    And user sets "Show Group Comparison" property of box plot viewer to "false"

  Scenario: Resize and auto layout
    Given user sets "Auto Layout" property of box plot viewer to "true"
    When user hovers over box plot viewer
    Then "Marker Color" column input in box plot viewer should be visible
    When user resizes box plot viewer to 170 by 150
    Then "Marker Color" column input in box plot viewer should be hidden
    When user restores the size of box plot viewer
    And user hovers over box plot viewer
    Then "Marker Color" column input in box plot viewer should be visible
    When user sets "Marker Color Column" property of box plot viewer to "SEX"
    And user resizes box plot viewer to 120 wide
    And user restores the size of box plot viewer
    Then no errors should have been logged
    When user sets "Marker Color Column" property of box plot viewer to ""

  Scenario: Marker gate and size scaling
    When user sets "Show Markers" property of box plot viewer to "false"
    Then box plot viewer should have less ink than before
    When user clicks on settings icon of box plot viewer
    Then "Marker Type" property in context panel should be disabled
    When user sets "Show Markers" property of box plot viewer to "true"
    Then box plot viewer should have more ink than before
    And "Marker Type" property in context panel should be enabled
    When user sets "Marker Size Column" property of box plot viewer to "WEIGHT"
    And user sets "Marker Size Scaling" property of box plot viewer to "logarithmic"
    Then box plot viewer should have repainted
    When user sets "Marker Size Scaling" property of box plot viewer to "linear"
    Then box plot viewer should have repainted
    When user sets "Marker Size Column" property of box plot viewer to ""

  Scenario: Whisker and control-band style
    When user sets "Whisker Line Width" property of box plot viewer to "4"
    Then box plot viewer should have repainted
    When user sets "Whisker Width Ratio" property of box plot viewer to "0.3"
    Then box plot viewer should have repainted
    When user sets "Control Band Color" property of box plot viewer to "#00AA00"
    Then no errors should have been logged
    When user sets properties of box plot viewer:
      | Whisker Line Width  | 2   |
      | Whisker Width Ratio | 0.5 |

  Scenario: Controls visibility
    Then "Show Size Selector" property of box plot viewer should be "false"
    And "Marker Size" column input in box plot viewer should be hidden
    When user sets "Show Size Selector" property of box plot viewer to "true"
    And user hovers over box plot viewer
    Then "Marker Size" column input in box plot viewer should be visible
    When user sets "Show Value Selector" property of box plot viewer to "false"
    Then Value column input in box plot viewer should be hidden
    When user sets "Show Color Selector" property of box plot viewer to "false"
    Then "Marker Color" column input in box plot viewer should be hidden
    When user sets "Show Category Selector" property of box plot viewer to "false"
    Then box plot viewer should have repainted by at least 2000 pixels
    When user sets "Show Value Axis" property of box plot viewer to "false"
    Then box plot viewer should have repainted by at least 2000 pixels
    And box plot viewer should not have a "y axis" area
    When user sets "Show Category Axis" property of box plot viewer to "false"
    Then box plot viewer should have repainted by at least 2000 pixels
    And box plot viewer should not have an "x axis" area
    When user sets "Show Category Axis" property of box plot viewer to "true"
    Then box plot viewer should have repainted by at least 2000 pixels
    And box plot viewer should have an "x axis" area
    When user sets "Show Value Axis" property of box plot viewer to "true"
    Then box plot viewer should have repainted by at least 2000 pixels
    And box plot viewer should have a "y axis" area
    When user sets "Show Category Selector" property of box plot viewer to "true"
    Then box plot viewer should have repainted by at least 2000 pixels
    When user sets properties of box plot viewer:
      | Show Value Selector    | true  |
      | Show Color Selector    | true  |
      | Show Size Selector     | false |
    Then box plot viewer should have repainted by at least 2000 pixels
    When user hovers over box plot viewer
    Then "Marker Size" column input in box plot viewer should be hidden
    And Value column input in box plot viewer should be visible
    And "Marker Color" column input in box plot viewer should be visible

  Scenario: Title and description
    When user sets properties of box plot viewer:
      | Show Title | true        |
      | Title      | Age by Race |
    Then title of box plot viewer should have text "Age by Race"
    When user sets properties of box plot viewer:
      | Description                 | Box plot of patient ages |
      | Description Visibility Mode | Always                   |
    Then description of box plot viewer should have text "Box plot of patient ages"
    When user sets "Description Position" property of box plot viewer to "Bottom"
    Then description of box plot viewer should be visible
    When user sets "Description Visibility Mode" property of box plot viewer to "Never"
    Then description of box plot viewer should be absent
    When user sets properties of box plot viewer:
      | Show Title                  | false |
      | Title                       |       |
      | Description                 |       |
      | Description Visibility Mode | Auto  |
      | Description Position        | Top   |

  Scenario: Axis font
    When user sets "Axis Font" property of box plot viewer to "normal normal 16px \"Roboto\""
    Then "Axis Font" property of box plot viewer should be "normal normal 16px \"Roboto\""
    And box plot viewer should have repainted
    When user sets "Axis Font" property of box plot viewer to "normal normal 10px \"Roboto\""
    Then no errors should have been logged

  Scenario: Date category mapping
    When user sets "Category 1" property of box plot viewer to "STARTED"
    Then "Category 1" property of box plot viewer should be "STARTED"
    When user sets "Category 1 Map" property of box plot viewer to "month"
    Then box plot viewer should have repainted
    And "Category 1 Map" property of box plot viewer should be "month"
    When user sets "Category 1 Map" property of box plot viewer to "quarter"
    Then box plot viewer should have repainted
    When user sets "Category 1" property of box plot viewer to "RACE"
    Then "Category 1" property of box plot viewer should be "RACE"
    And no errors should have been logged
    When user sets "Category 1" property of box plot viewer to "SEX"

  Scenario: Custom tooltip
    When user sets properties of box plot viewer:
      | Category 1   | RACE                |
      | Marker Size  | 10                  |
      | Row Tooltip  | AGE\nSEX\nWEIGHT    |
      | Show Tooltip | show custom tooltip |
    Then "Row Tooltip" property of box plot viewer should be "AGE\nSEX\nWEIGHT"
    When user hovers over the "marker" area of box plot viewer
    Then the tooltip should show columns "AGE, SEX, WEIGHT"
    When user moves the pointer away from box plot viewer
    Then tooltip should be hidden
    When user sets properties of box plot viewer:
      | Show Tooltip | inherit from table |
      | Row Tooltip  |                    |
      | Category 1   | SEX                |
    And user hovers over the "marker" area of box plot viewer
    Then tooltip should be visible
    And the tooltip should show some columns
    And the tooltip should not show columns "AGE, SEX, WEIGHT"
    When user moves the pointer away from box plot viewer

  Scenario: Table switching resets Category 2
    Given user opens spgi dataset
    And user switches to the "demog-1000" table view
    When user sets "Category 2" property of box plot viewer to "RACE"
    Then "Category 2" property of box plot viewer should be "RACE"
    When user sets "Table" property of box plot viewer to "spgi-100"
    Then "Category 2" property of box plot viewer should be ""
    And no errors should have been logged
    When user sets properties of box plot viewer:
      | Value      | Average Mass |
      | Category 1 | Series       |
    Then "Table" property of box plot viewer should be "spgi-100"
    And box plot viewer should be bound to table "spgi-100"
    And box plot viewer should have a "category Triazoles" area
    When user sets "Table" property of box plot viewer to "demog-1000"
    And user sets properties of box plot viewer:
      | Value      | AGE |
      | Category 1 | SEX |
      | Category 2 |     |
    Then "Table" property of box plot viewer should be "demog-1000"
    And box plot viewer should be bound to table "demog-1000"

  Scenario: Coloring keeps the render valid
    When user sets "Marker Color Column" property of box plot viewer to "RACE"
    Then "Marker Color Column" property of box plot viewer should be "RACE"
    And box plot viewer should be painted
    And no errors should have been logged
    When user sets "Marker Color Column" property of box plot viewer to ""

  Scenario: Double-click resets the view
    Given user listens for "d4-boxplot-reset-view" event on box plot viewer
    And user remembers the value range of box plot viewer
    When user zooms into the value axis of box plot viewer
    Then box plot viewer should show a narrower value range than before
    When user double-clicks on empty plot space of box plot viewer
    Then "d4-boxplot-reset-view" event should have fired on box plot viewer
    And box plot viewer should show the remembered value range
