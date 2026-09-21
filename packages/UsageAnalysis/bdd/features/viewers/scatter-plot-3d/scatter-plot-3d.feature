@journey @viewers @realizes:viewers.scatter-plot3d
Feature: 3D scatter plot
  The 3D scatter plot through its own readings — the WebGL frame's signature, the camera, the rows
  drawn, the current row, the rows a mouse-over group highlights — and its hit areas: the axes are
  assigned automatically and follow the column properties, a categorical color draws a legend of
  the values and a numeric one takes it away, marker type and opacity, axes and a logarithmic axis
  redraw the scene, a drag rotates the camera and Reset View brings it home, the wheel zooms, a
  click on a point sets the current row and Shift-click selects it, filtered-out points come back
  as ghosts, a hover on a bar chart highlights the matching points unless Show Mouse Over Row
  Group is off, and the legend docks where Legend Position says. One journey on demog-1000.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a 3d scatter plot viewer
    Then the "rows shown" reading of 3d scatter plot viewer should be 872

  Scenario: The axes are assigned automatically
    Then X column input in 3d scatter plot viewer should contain text "AGE"
    And Y column input in 3d scatter plot viewer should contain text "HEIGHT"
    And Z column input in 3d scatter plot viewer should contain text "WEIGHT"
    And properties of 3d scatter plot viewer should be:
      | X | AGE    |
      | Y | HEIGHT |
      | Z | WEIGHT |
    And no errors should have been logged

  Scenario: Reassigning X and Z moves the selectors and redraws the scene
    When user sets properties of 3d scatter plot viewer:
      | X | WEIGHT |
      | Z | AGE    |
    Then X column input in 3d scatter plot viewer should contain text "WEIGHT"
    And Z column input in 3d scatter plot viewer should contain text "AGE"
    And the "scene signature" reading of 3d scatter plot viewer should differ from before
    When user sets properties of 3d scatter plot viewer:
      | X | AGE    |
      | Z | WEIGHT |
    Then the "scene signature" reading of 3d scatter plot viewer should differ from before
    And no errors should have been logged

  Scenario: Color by a category draws a legend of its values
    Then legend of 3d scatter plot viewer should be hidden
    When user sets "Color" property of 3d scatter plot viewer to "SEX"
    Then legend of 3d scatter plot viewer should have 2 items
    And the "scene signature" reading of 3d scatter plot viewer should differ from before
    And no errors should have been logged

  Scenario: Color by a number takes the legend away
    When user sets "Color" property of 3d scatter plot viewer to "AGE"
    Then legend of 3d scatter plot viewer should be hidden
    And the "scene signature" reading of 3d scatter plot viewer should differ from before
    And no errors should have been logged

  Scenario: Marker type redraws the markers
    When user sets "Marker Type" property of 3d scatter plot viewer to "box"
    Then the "scene signature" reading of 3d scatter plot viewer should differ from before
    When user sets "Marker Type" property of 3d scatter plot viewer to "sphere"
    Then the "scene signature" reading of 3d scatter plot viewer should differ from before
    When user sets "Marker Type" property of 3d scatter plot viewer to "cylinder"
    Then the "scene signature" reading of 3d scatter plot viewer should differ from before
    When user sets "Marker Type" property of 3d scatter plot viewer to "octahedron"
    Then the "scene signature" reading of 3d scatter plot viewer should differ from before
    And no errors should have been logged

  Scenario: Marker opacity redraws the markers
    When user sets "Marker Opacity" property of 3d scatter plot viewer to "25"
    Then the "scene signature" reading of 3d scatter plot viewer should differ from before
    When user sets "Marker Opacity" property of 3d scatter plot viewer to "69"
    Then the "scene signature" reading of 3d scatter plot viewer should differ from before
    And no errors should have been logged

  Scenario: Show Axes hides and restores the axes
    When user sets "Show Axes" property of 3d scatter plot viewer to "false"
    Then the "scene signature" reading of 3d scatter plot viewer should differ from before
    When user sets "Show Axes" property of 3d scatter plot viewer to "true"
    Then the "scene signature" reading of 3d scatter plot viewer should differ from before
    And no errors should have been logged

  Scenario: A logarithmic X axis re-scales the scene without error
    When user sets "X Axis Type" property of 3d scatter plot viewer to "logarithmic"
    Then the "scene signature" reading of 3d scatter plot viewer should differ from before
    And no errors should have been logged
    And no error or warning balloon should have been shown
    When user sets "X Axis Type" property of 3d scatter plot viewer to "linear"
    Then the "scene signature" reading of 3d scatter plot viewer should differ from before
    And no errors should have been logged

  Scenario: A drag rotates the scene and Reset View brings the camera home
    When user drags across the "view" area of 3d scatter plot viewer
    Then the "camera x" reading of 3d scatter plot viewer should differ from before
    And the "scene signature" reading of 3d scatter plot viewer should differ from before
    When user picks "Reset View" from the context menu of 3d scatter plot viewer
    Then the "camera x" reading of 3d scatter plot viewer should be 0
    And the "camera y" reading of 3d scatter plot viewer should be 0
    And the "camera distance" reading of 3d scatter plot viewer should be 4
    And the "scene signature" reading of 3d scatter plot viewer should differ from before
    And no errors should have been logged

  Scenario: The mouse wheel zooms the scene
    When user scrolls the mouse wheel up over the "view" area of 3d scatter plot viewer
    Then the "camera distance" reading of 3d scatter plot viewer should be lower than before
    And the "scene signature" reading of 3d scatter plot viewer should differ from before
    When user scrolls the mouse wheel down over the "view" area of 3d scatter plot viewer
    Then the "camera distance" reading of 3d scatter plot viewer should be higher than before
    And no errors should have been logged

  Scenario: A click makes a row current, Shift-click selects it
    When user clicks on the "point" area of 3d scatter plot viewer
    Then the "current row" reading of 3d scatter plot viewer should differ from before
    And the table should have a current row
    When user clears the row selection
    And user clicks on the "point" area of 3d scatter plot viewer holding Shift
    Then some rows should be selected
    And the "scene signature" reading of 3d scatter plot viewer should differ from before
    When user clears the row selection
    Then no errors should have been logged

  Scenario: Show Filtered Out Points brings the filtered-away rows back
    When user filters rows where "SEX" is "F"
    Then 553 rows should pass the filter
    And the "rows shown" reading of 3d scatter plot viewer should be lower than before
    When user sets "Show Filtered Out Points" property of 3d scatter plot viewer to "true"
    Then the "rows shown" reading of 3d scatter plot viewer should be 872
    And the "scene signature" reading of 3d scatter plot viewer should differ from before
    When user sets "Show Filtered Out Points" property of 3d scatter plot viewer to "false"
    Then the "rows shown" reading of 3d scatter plot viewer should be lower than before
    When user resets the filter
    Then the "rows shown" reading of 3d scatter plot viewer should be 872
    And no errors should have been logged

  Scenario: A hover on a bar chart highlights the matching points
    When user adds a bar chart viewer with:
      | Split | SEX |
    Then bar chart viewer should have a "bar F" area
    And "Show Mouse Over Row Group" property of 3d scatter plot viewer should be "true"
    And the "highlighted rows" reading of 3d scatter plot viewer should be 0
    When user takes a snapshot of 3d scatter plot viewer
    And user hovers over the "bar F" area of bar chart viewer
    Then the "highlighted rows" reading of 3d scatter plot viewer should be higher than before
    And the "scene signature" reading of 3d scatter plot viewer should differ from before
    When user moves the pointer away from bar chart viewer
    And user sets "Show Mouse Over Row Group" property of 3d scatter plot viewer to "false"
    And user hovers over the "bar F" area of bar chart viewer
    Then the "highlighted rows" reading of 3d scatter plot viewer should be 0
    When user moves the pointer away from bar chart viewer
    And user sets "Show Mouse Over Row Group" property of 3d scatter plot viewer to "true"
    And user clicks on close icon of bar chart viewer
    Then bar chart viewer should be absent
    And no errors should have been logged

  Scenario: Legend Position docks the legend on the side it names
    When user sets properties of 3d scatter plot viewer:
      | Color             | SEX    |
      | Legend Visibility | Always |
    Then legend of 3d scatter plot viewer should have 2 items
    When user sets "Legend Position" property of 3d scatter plot viewer to "Left"
    Then the legend of 3d scatter plot viewer should be on the left
    When user sets "Legend Position" property of 3d scatter plot viewer to "Right"
    Then the legend of 3d scatter plot viewer should be on the right
    When user sets properties of 3d scatter plot viewer:
      | Legend Position   | Auto |
      | Legend Visibility | Auto |
      | Color             |      |
    Then legend of 3d scatter plot viewer should be hidden
    And no errors should have been logged
