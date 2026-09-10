@journey @viewers @realizes:viewers.pie-chart
Feature: Pie chart property surface
  Every setting of the pie chart with a claim about the disc it drew, never a property read back on
  its own: the wedges and their shares, the two sort orders, the rotation Start Angle applies, the
  radius Max Radius caps, the explosion Shift produces, where the labels go and how many are drawn,
  the outline, the on-chart column selector, the margins, the size below which the labels are
  dropped, and Donut mode with its hole and centre label. Selector and legend, aggregations,
  clicks, nulls and persistence have features of their own, and Row Source is the cross-viewer
  `viewers/row-source` journey. One journey on demog-1000 with a pie chart of RACE — Caucasian 896
  rows (89.6 %), Other 62, Black 27, Asian 15 (1.5 %) — and every scenario puts back what it
  changed.
  Not translated: Outline Line Width 0, because the canvas ignores a line width of 0 and keeps the
  previous one, so there is no honest paint claim for it; and the property echoes of the old
  pie-chart-spec (Appearance, Labels, Outline, Column selector, Legend, Title and description,
  Selection, Auto layout), which are folded into the geometric claims below.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a pie chart viewer with:
      | Category | RACE |
    Then the "slices" reading of pie chart viewer should be 4
    And pie chart viewer should show 1000 rows

  Scenario: A wedge per category, sized by its share of the rows
    Then pie chart viewer should be painted
    And pie chart viewer should have a "slice Caucasian" area
    And pie chart viewer should have a "slice Asian" area
    And pie chart viewer should have a "slice Black" area
    And pie chart viewer should have a "slice Other" area
    And the "share of Caucasian" reading of pie chart viewer should be 89.6
    And the "share of Asian" reading of pie chart viewer should be 1.5
    And the "angle value of Caucasian" reading of pie chart viewer should be 896
    And the "angle value of Asian" reading of pie chart viewer should be 15
    And the "pie" area of pie chart viewer should be painted in at least 4 colors
    And no errors should have been logged

  Scenario: Sorting by value puts the smallest wedge first and the order reverses
    Then "Pie Sort Type" property of pie chart viewer should be "by value"
    And the slices of pie chart viewer should be ordered by share ascending
    And the "start angle of Asian" reading of pie chart viewer should be 0
    When user sets "Pie Sort Order" property of pie chart viewer to "desc"
    Then the slices of pie chart viewer should be ordered by share descending
    And the "start angle of Caucasian" reading of pie chart viewer should be 0
    And pie chart viewer should have repainted by at least 500 pixels
    When user sets "Pie Sort Order" property of pie chart viewer to "asc"
    Then the slices of pie chart viewer should be ordered by share ascending
    And the "start angle of Asian" reading of pie chart viewer should be 0
    And no errors should have been logged

  Scenario: Sorting by category draws the wedges alphabetically
    When user sets "Pie Sort Type" property of pie chart viewer to "by category"
    Then the slices of pie chart viewer should be ordered by category
    And the "start angle of Asian" reading of pie chart viewer should be 0
    And pie chart viewer should have repainted by at least 500 pixels
    When user sets "Pie Sort Type" property of pie chart viewer to "by value"
    Then the slices of pie chart viewer should be ordered by share ascending
    And pie chart viewer should have repainted by at least 500 pixels
    And no errors should have been logged

  Scenario: Start Angle rotates the whole disc
    Then the "start angle of Asian" reading of pie chart viewer should be 0
    When user sets "Start Angle" property of pie chart viewer to "90"
    Then the "start angle of Asian" reading of pie chart viewer should be 90
    And pie chart viewer should have repainted by at least 500 pixels
    When user sets "Start Angle" property of pie chart viewer to "180"
    Then the "start angle of Asian" reading of pie chart viewer should be 180
    And pie chart viewer should have repainted by at least 500 pixels
    When user sets "Start Angle" property of pie chart viewer to "0"
    Then the "start angle of Asian" reading of pie chart viewer should be 0
    And pie chart viewer should have repainted by at least 500 pixels
    And no errors should have been logged

  Scenario: Max Radius caps the disc and gives it back
    Then the "pie radius" reading of pie chart viewer should be 150
    When user sets "Max Radius" property of pie chart viewer to "100"
    Then the "pie radius" reading of pie chart viewer should be 100
    And the "outer radius of Caucasian" reading of pie chart viewer should be 100
    And the "pie" area of pie chart viewer should be narrower than before
    And pie chart viewer should have less ink than before
    When user sets "Max Radius" property of pie chart viewer to "150"
    Then the "pie radius" reading of pie chart viewer should be 150
    And the "pie" area of pie chart viewer should be wider than before
    And pie chart viewer should have more ink than before
    And no errors should have been logged

  Scenario: Shift explodes the wedges out of the centre
    Then the slices of pie chart viewer should sit 0 pixels off the centre
    When user sets "Shift" property of pie chart viewer to "20"
    Then the slices of pie chart viewer should sit 20 pixels off the centre
    And pie chart viewer should have repainted by at least 500 pixels
    When user sets "Shift" property of pie chart viewer to "0"
    Then the slices of pie chart viewer should sit 0 pixels off the centre
    And pie chart viewer should have repainted by at least 500 pixels
    And no errors should have been logged

  Scenario: Label Position moves the labels out of the wedges and back in
    Then pie chart viewer should have a "label of Caucasian" area
    And pie chart viewer should have a "label of Asian" area
    And the "labels shown" reading of pie chart viewer should be 4
    When user sets "Label Position" property of pie chart viewer to "Inside"
    Then pie chart viewer should not have a "label of Caucasian" area
    And the "labels shown" reading of pie chart viewer should be 3
    And pie chart viewer should have repainted by at least 300 pixels
    When user sets "Label Position" property of pie chart viewer to "Outside"
    Then pie chart viewer should have a "label of Caucasian" area
    And pie chart viewer should have a "label of Asian" area
    And the "labels shown" reading of pie chart viewer should be 4
    And pie chart viewer should have repainted by at least 300 pixels
    When user sets "Label Position" property of pie chart viewer to "Auto"
    Then pie chart viewer should have a "label of Caucasian" area
    And the "labels shown" reading of pie chart viewer should be 4
    And no errors should have been logged

  Scenario: With nothing left to say a wedge gets no label at all
    Then the "labels shown" reading of pie chart viewer should be 4
    When user sets properties of pie chart viewer:
      | Show Label      | false |
      | Show Percentage | false |
    Then the "labels shown" reading of pie chart viewer should be 0
    And pie chart viewer should not have a "label of Caucasian" area
    And pie chart viewer should have repainted by at least 300 pixels
    When user sets "Show Value" property of pie chart viewer to "true"
    Then the "labels shown" reading of pie chart viewer should be higher than before
    And pie chart viewer should have repainted by at least 300 pixels
    When user sets properties of pie chart viewer:
      | Show Label      | true  |
      | Show Percentage | true  |
      | Show Value      | false |
    Then the "labels shown" reading of pie chart viewer should be 4
    And pie chart viewer should have repainted by at least 300 pixels
    And no errors should have been logged

  Scenario: Outline Line Width thickens the border between the wedges
    When user sets "Outline Line Width" property of pie chart viewer to "5"
    Then "Outline Line Width" property of pie chart viewer should be "5"
    And pie chart viewer should have repainted by at least 300 pixels
    When user sets "Outline Line Width" property of pie chart viewer to "1"
    Then "Outline Line Width" property of pie chart viewer should be "1"
    And pie chart viewer should have repainted by at least 300 pixels
    And no errors should have been logged

  Scenario: Show Column Selector takes the on-chart selector off the canvas
    Then pie chart viewer should have a "column selector" area
    When user sets "Show Column Selector" property of pie chart viewer to "false"
    Then pie chart viewer should not have a "column selector" area
    And pie chart viewer should have repainted by at least 300 pixels
    When user sets "Show Column Selector" property of pie chart viewer to "true"
    Then pie chart viewer should have a "column selector" area
    And pie chart viewer should have repainted by at least 300 pixels
    And no errors should have been logged

  Scenario: Auto Layout off hands the margins to the chart box
    Then the "pie radius" reading of pie chart viewer should be 150
    When user sets properties of pie chart viewer:
      | Auto Layout | false |
      | Margin Left | 100   |
      | Margin Top  | 100   |
    Then the "pie radius" reading of pie chart viewer should be lower than before
    And the "view" area of pie chart viewer should be narrower than before
    And pie chart viewer should have repainted by at least 500 pixels
    When user sets properties of pie chart viewer:
      | Auto Layout | true |
      | Margin Left | 10   |
      | Margin Top  | 10   |
    Then the "pie radius" reading of pie chart viewer should be 150
    And the "view" area of pie chart viewer should be wider than before
    And no errors should have been logged

  Scenario: A viewer too small for its labels drops them and the selector
    Then the "labels shown" reading of pie chart viewer should be 4
    And pie chart viewer should have a "column selector" area
    When user resizes pie chart viewer to 150 by 150
    Then the "labels shown" reading of pie chart viewer should be 0
    And pie chart viewer should not have a "column selector" area
    And the "pie radius" reading of pie chart viewer should be lower than before
    And the "slices" reading of pie chart viewer should be 4
    When user restores the size of pie chart viewer
    Then the "labels shown" reading of pie chart viewer should be 4
    And pie chart viewer should have a "column selector" area
    And the "pie radius" reading of pie chart viewer should be 150
    And no errors should have been logged

  Scenario: The Table property rebinds the chart to another table
    Given user opens spgi dataset
    And user switches to the "demog-1000" table view
    When user sets "Table" property of pie chart viewer to "spgi-100"
    Then pie chart viewer should be bound to table "spgi-100"
    When user sets "Category" property of pie chart viewer to "Primary Series Name"
    Then the "slices" reading of pie chart viewer should be 5
    And pie chart viewer should show 100 rows
    And pie chart viewer should be painted
    When user sets properties of pie chart viewer:
      | Table    | demog-1000 |
      | Category | RACE       |
    Then pie chart viewer should be bound to table "demog-1000"
    And the "slices" reading of pie chart viewer should be 4
    And pie chart viewer should show 1000 rows
    And the "share of Caucasian" reading of pie chart viewer should be 89.6
    And no errors should have been logged

  Scenario: Donut mode opens a hole and puts the column name in it
    Then pie chart viewer should not have a "donut hole" area
    And pie chart viewer should not have a "centre label" area
    When user sets "Mode" property of pie chart viewer to "Donut"
    Then pie chart viewer should have a "donut hole" area
    And pie chart viewer should have a "centre label" area
    And the "pie radius" reading of pie chart viewer should be the same as before
    And pie chart viewer should have less ink than before
    And pie chart viewer should have repainted by at least 1000 pixels
    When user sets "Center Label" property of pie chart viewer to "Race mix"
    Then the "centre label" area of pie chart viewer should have repainted
    When user sets "Show Center Label" property of pie chart viewer to "false"
    Then pie chart viewer should not have a "centre label" area
    And pie chart viewer should have a "donut hole" area
    And pie chart viewer should have repainted
    When user sets properties of pie chart viewer:
      | Show Center Label | true |
      | Center Label      |      |
      | Mode              | Pie  |
    Then pie chart viewer should not have a "donut hole" area
    And pie chart viewer should not have a "centre label" area
    And pie chart viewer should have more ink than before
    And no errors should have been logged
