@journey @viewers @realizes:viewers.scatter-plot3d
Feature: 3D scatter plot look
  The look properties the TestTrack "Manual scenarios" list and scatter-plot-3d.feature leaves out,
  each claimed on the frame's `scene signature` (the WebGL picture has no readable pixels): Y and Z
  go logarithmic and come back to the very frame they left, Size by a number and Marker Size Scaling
  redraw the markers, a label column draws labels and clearing it gives the unlabelled frame back,
  Marker Random Rotation turns the markers and switching it off puts them back, a Back Color picked
  in the context panel repaints the scene and white gives the first frame back, and hiding an
  on-viewer column selector takes it off the plot and gives its strip to the scene. That scenario
  first hides and shows the X selector once: until the first layout change the canvas keeps the size
  it was created with, taller than its area (MISSING.md, Answered), so a size compared with the
  first frame says nothing. A logarithmic
  axis over zero values cannot be drawn and says why in an error balloon (cars: diesel, turbo and
  two.doors are 0/1 columns). Not here, see MISSING.md: Axis Line Color and the grid line toggles,
  which the renderer never receives (the frame does not change, core#929), and the log-scale error as a
  reading. Dynamic Camera Movement is not claimed (an endless rotation, nothing to settle on).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a 3d scatter plot viewer
    Then the "rows shown" reading of 3d scatter plot viewer should be 872

  Scenario: Y and Z axes go logarithmic and come back
    When user remembers the "scene signature" reading of 3d scatter plot viewer
    And user sets "Y Axis Type" property of 3d scatter plot viewer to "logarithmic"
    Then the "scene signature" reading of 3d scatter plot viewer should differ from before
    When user sets "Z Axis Type" property of 3d scatter plot viewer to "logarithmic"
    Then the "scene signature" reading of 3d scatter plot viewer should differ from before
    And the "rows shown" reading of 3d scatter plot viewer should be 872
    And no error or warning balloon should have been shown
    When user sets "Y Axis Type" property of 3d scatter plot viewer to "linear"
    Then the "scene signature" reading of 3d scatter plot viewer should differ from before
    When user sets "Z Axis Type" property of 3d scatter plot viewer to "linear"
    Then the "scene signature" reading of 3d scatter plot viewer should be as remembered
    And no errors should have been logged

  Scenario: Size by a number scales the markers and Marker Size Scaling redraws them
    When user remembers the "scene signature" reading of 3d scatter plot viewer
    And user sets "Size" property of 3d scatter plot viewer to "WEIGHT"
    Then the "scene signature" reading of 3d scatter plot viewer should differ from before
    When user sets "Marker Size Scaling" property of 3d scatter plot viewer to "logarithmic"
    Then the "scene signature" reading of 3d scatter plot viewer should differ from before
    When user sets "Marker Size Scaling" property of 3d scatter plot viewer to "linear"
    Then the "scene signature" reading of 3d scatter plot viewer should differ from before
    When user sets "Size" property of 3d scatter plot viewer to ""
    Then the "scene signature" reading of 3d scatter plot viewer should be as remembered
    And no errors should have been logged

  Scenario: A label column draws labels and clearing it takes them away
    When user remembers the "scene signature" reading of 3d scatter plot viewer
    And user sets "Label" property of 3d scatter plot viewer to "SEX"
    Then the "scene signature" reading of 3d scatter plot viewer should differ from before
    When user sets "Label" property of 3d scatter plot viewer to ""
    Then the "scene signature" reading of 3d scatter plot viewer should be as remembered
    And no errors should have been logged

  Scenario: Random marker rotation turns the markers and switching it off puts them back
    When user remembers the "scene signature" reading of 3d scatter plot viewer
    And user sets "Marker Random Rotation" property of 3d scatter plot viewer to "true"
    Then the "scene signature" reading of 3d scatter plot viewer should differ from before
    When user sets "Marker Random Rotation" property of 3d scatter plot viewer to "false"
    Then the "scene signature" reading of 3d scatter plot viewer should be as remembered
    And no errors should have been logged

  Scenario: Back Color picked in the context panel repaints the scene behind the markers
    When user remembers the "scene signature" reading of 3d scatter plot viewer
    And user picks "Properties..." from the context menu of 3d scatter plot viewer
    Then the context panel should show "3d scatter plot"
    Given "Style" category in context panel is expanded
    When user clicks on editor of "Back Color" property in context panel
    And user picks the color "#FFA500" in the color picker
    Then "Back Color" property of 3d scatter plot viewer should be "#FFA500"
    And the "scene signature" reading of 3d scatter plot viewer should not be as remembered
    When user sets "Back Color" property of 3d scatter plot viewer to "#FFFFFF"
    Then the "scene signature" reading of 3d scatter plot viewer should be as remembered
    And no errors should have been logged

  Scenario: Hidden column selectors leave the plot and give their strip to the scene
    When user sets "Show X Selector" property of 3d scatter plot viewer to "false"
    And user sets "Show X Selector" property of 3d scatter plot viewer to "true"
    Then X column input in 3d scatter plot viewer should be visible
    When user sets "Show X Selector" property of 3d scatter plot viewer to "false"
    Then X column input in 3d scatter plot viewer should be hidden
    And the "view" area of 3d scatter plot viewer should be taller than before
    When user sets "Show Color Selector" property of 3d scatter plot viewer to "false"
    Then Color column input in 3d scatter plot viewer should be hidden
    And the "view" area of 3d scatter plot viewer should be taller than before
    When user sets "Show Y Selector" property of 3d scatter plot viewer to "false"
    Then Y column input in 3d scatter plot viewer should be hidden
    And the "view" area of 3d scatter plot viewer should be wider than before
    When user sets "Show Z Selector" property of 3d scatter plot viewer to "false"
    Then Z column input in 3d scatter plot viewer should be hidden
    And the "view" area of 3d scatter plot viewer should be wider than before
    When user sets properties of 3d scatter plot viewer:
      | Show X Selector     | true |
      | Show Y Selector     | true |
      | Show Z Selector     | true |
      | Show Color Selector | true |
    Then X column input in 3d scatter plot viewer should be visible
    And Y column input in 3d scatter plot viewer should be visible
    And Z column input in 3d scatter plot viewer should be visible
    And Color column input in 3d scatter plot viewer should be visible
    And the "view" area of 3d scatter plot viewer should be narrower than before
    And no errors should have been logged

  Scenario: A logarithmic axis over zero values says why it cannot be drawn
    Given user opens cars dataset
    And user adds a 3d scatter plot viewer
    Then Y column input in 3d scatter plot viewer should contain text "turbo"
    When user sets "Y Axis Type" property of 3d scatter plot viewer to "logarithmic"
    Then an error balloon containing "'turbo' contains negative or zero values and can't be visualized in the log scale" should have been shown
    When user sets "Y Axis Type" property of 3d scatter plot viewer to "linear"
    Then the "rows shown" reading of 3d scatter plot viewer should be 30
    And no errors should have been logged
