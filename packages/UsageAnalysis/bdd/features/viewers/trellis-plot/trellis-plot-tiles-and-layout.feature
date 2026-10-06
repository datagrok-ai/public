@journey @viewers @realizes:viewers.trellis-plot
Feature: Trellis plot tiled view and auto layout
  A trellis split by one column lays its cells out as tiles: a rectangle Tiles Per Row wide, padded
  to a full row, against the single band it draws with tiling off. Then what auto layout drops as
  the viewer shrinks — the control panel, the column selectors, each label strip on its own
  threshold — what it keeps when it is switched off, the title and the description in each of the
  four slots, and the full-screen icon that only exists in the cell the pointer is in. One journey
  on demog-1000; every scenario puts back what it changed.

  Packing is off, so all six DIS_POP categories take part: with packing on the trellis drops UC,
  whose every row is missing the HEIGHT or the WEIGHT the inner scatter plot draws. UC is the cell
  that is drawn and stays blank in the two rungs wide enough to show it.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a trellis plot viewer with:
      | X Column Names  | DIS_POP      |
      | Y Column Names  |              |
      | Viewer Type     | Scatter plot |
      | Pack Categories | false        |
    Then the "x categories" reading of trellis plot viewer should be 6
    And "Tiles" property of trellis plot viewer should be "true"

  Scenario: Tiled view rebuilds the band into a padded rectangle
    Then "Tiles Per Row" property of trellis plot viewer should be "4"
    And the cells of trellis plot viewer should be 4 wide and 2 tall
    And the "cells" reading of trellis plot viewer should be 8
    And the "cells drawn" reading of trellis plot viewer should be 6
    And the "blank cells" reading of trellis plot viewer should be 1
    And trellis plot viewer should have a "cell RA" area
    And trellis plot viewer should have a "cell UC" area
    And the "distinct cell signatures" reading of trellis plot viewer should be 6
    And no errors should have been logged

  Scenario: Tiles Per Row 1 gives a single-column strip
    When user sets "Tiles Per Row" property of trellis plot viewer to "1"
    Then the cells of trellis plot viewer should be 1 wide and 5 tall
    And the "cells" reading of trellis plot viewer should be 5
    And the "cells drawn" reading of trellis plot viewer should be 5
    And the "blank cells" reading of trellis plot viewer should be 0
    And trellis plot viewer should have a "cell RA" area
    And trellis plot viewer should not have a "cell UC" area
    And no errors should have been logged

  Scenario: Tiles Per Row 3 gives three tiles a row
    When user sets "Tiles Per Row" property of trellis plot viewer to "3"
    Then the cells of trellis plot viewer should be 3 wide and 2 tall
    And the "cells" reading of trellis plot viewer should be 6
    And the "cells drawn" reading of trellis plot viewer should be 6
    And trellis plot viewer should have a "cell UC" area
    And no errors should have been logged

  Scenario: Tiles off returns a single band
    When user sets "Tiles" property of trellis plot viewer to "false"
    Then the cells of trellis plot viewer should be 5 wide and 1 tall
    And the "cells" reading of trellis plot viewer should be 5
    And the "cells drawn" reading of trellis plot viewer should be 5
    And the "blank cells" reading of trellis plot viewer should be 0
    And trellis plot viewer should not have a "cell UC" area
    When user sets properties of trellis plot viewer:
      | Tiles         | true |
      | Tiles Per Row | 4    |
    Then the cells of trellis plot viewer should be 4 wide and 2 tall
    And no errors should have been logged

  Scenario: Auto layout drops the control panel when the viewer gets short
    When user resizes trellis plot viewer to 900 by 700
    Then trellis plot viewer should have a "control panel" area
    When user resizes trellis plot viewer to 900 by 280
    Then trellis plot viewer should not have a "control panel" area
    And "Show Control Panel" property of trellis plot viewer should be "true"
    When user restores the size of trellis plot viewer
    Then trellis plot viewer should have a "control panel" area
    And no errors should have been logged

  Scenario: With auto layout off the control panel stays at any size
    When user sets "Auto Layout" property of trellis plot viewer to "false"
    And user resizes trellis plot viewer to 900 by 280
    Then trellis plot viewer should have a "control panel" area
    When user sets "Show Control Panel" property of trellis plot viewer to "false"
    Then trellis plot viewer should not have a "control panel" area
    When user sets properties of trellis plot viewer:
      | Show Control Panel | true |
      | Auto Layout        | true |
    And user restores the size of trellis plot viewer
    Then trellis plot viewer should have a "control panel" area
    And no errors should have been logged

  Scenario: The two label strips and the two selector strips drop on their own thresholds
    When user sets properties of trellis plot viewer:
      | X Column Names | SEX  |
      | Y Column Names | RACE |
    Then the "x labels shown" reading of trellis plot viewer should be 2
    And the "y labels shown" reading of trellis plot viewer should be 4
    And trellis plot viewer should have an "x selectors" area
    And trellis plot viewer should have a "y selectors" area
    When user resizes trellis plot viewer to 240 by 400
    Then trellis plot viewer should not have an "x selectors" area
    And trellis plot viewer should have a "y selectors" area
    When user resizes trellis plot viewer to 900 by 200
    Then the "x labels shown" reading of trellis plot viewer should be 2
    And the "y labels shown" reading of trellis plot viewer should be 0
    And trellis plot viewer should have an "x selectors" area
    And trellis plot viewer should not have a "y selectors" area
    When user resizes trellis plot viewer to 400 by 180
    Then the "x labels shown" reading of trellis plot viewer should be 0
    And the "y labels shown" reading of trellis plot viewer should be 0
    And the "cells" reading of trellis plot viewer should be 8
    When user restores the size of trellis plot viewer
    Then the "x labels shown" reading of trellis plot viewer should be 2
    And the "y labels shown" reading of trellis plot viewer should be 4
    And no errors should have been logged
    When user sets properties of trellis plot viewer:
      | X Column Names | DIS_POP |
      | Y Column Names |         |
    Then the cells of trellis plot viewer should be 4 wide and 2 tall

  Scenario: Title and description
    When user sets properties of trellis plot viewer:
      | Show Title | true         |
      | Title      | Demographics |
    Then title of trellis plot viewer should have text "Demographics"
    When user sets properties of trellis plot viewer:
      | Description                 | By race and sex |
      | Description Visibility Mode | Always          |
      | Description Position        | Top             |
    Then description of trellis plot viewer should have text "By race and sex"
    And the "description slot" reading of trellis plot viewer should be "top"
    When user sets "Description Position" property of trellis plot viewer to "Bottom"
    Then the "description slot" reading of trellis plot viewer should be "bottom"
    When user sets "Description Position" property of trellis plot viewer to "Left"
    Then the "description slot" reading of trellis plot viewer should be "left"
    When user sets "Description Position" property of trellis plot viewer to "Right"
    Then the "description slot" reading of trellis plot viewer should be "right"
    When user sets "Description Visibility Mode" property of trellis plot viewer to "Never"
    Then the "description slot" reading of trellis plot viewer should be ""
    And description of trellis plot viewer should be absent
    When user sets properties of trellis plot viewer:
      | Title                       |      |
      | Description                 |      |
      | Description Visibility Mode | Auto |
      | Description Position        | Top  |
    Then no errors should have been logged

  Scenario: The full screen icon lives in the cell the pointer is in
    When user sets properties of trellis plot viewer:
      | X Column Names | SEX  |
      | Y Column Names | RACE |
    Then the cells of trellis plot viewer should be 2 wide and 4 tall
    When user moves the pointer away from trellis plot viewer
    Then trellis plot viewer should not have a "full screen icon" area
    When user hovers over the "cell body F | Caucasian" area of trellis plot viewer
    Then trellis plot viewer should have a "full screen icon" area
    When user moves the pointer away from trellis plot viewer
    Then trellis plot viewer should not have a "full screen icon" area
    When user hovers over the "cell body M | Asian" area of trellis plot viewer
    Then trellis plot viewer should have a "full screen icon" area
    When user sets "Allow Viewer Full Screen" property of trellis plot viewer to "false"
    And user hovers over the "cell body F | Caucasian" area of trellis plot viewer
    Then trellis plot viewer should not have a "full screen icon" area
    When user sets "Allow Viewer Full Screen" property of trellis plot viewer to "true"
    And user moves the pointer away from trellis plot viewer
    Then no errors should have been logged
