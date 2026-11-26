@guide @tutorial
Feature: How-to · Annotations
  For the maintainer: an annotation region on a scatter plot of iris, Petal.Length against Petal.Width, colored by
  Species. The plot's context menu, Tools > Draw Annotation Region, puts it in drawing mode; a drag from (4.45, 1.85)
  to (5.15, 1.35) in the units of the axes draws a rectangle around the points where versicolor and virginica mix
  (18 versicolor and 8 virginica, 26 rows: the half-step margins keep every point of that 0.1 grid inside), and the
  Formula Lines dialog opens on it. Its Title ("Overlap") and Description are typed there; OK keeps it on the viewer
  (its Annotation Regions property, in the browser only). Pointing at a marker-free spot inside it (4.82, 1.45,
  below the title drawn at its middle) shows the title, the description and "26 rows"; a click there selects those 26 rows. iris is opened and the plot
  added off camera; nothing is made on the server.

  Scenario: Mark a group of points with an annotation region
    Draw a region around points on a scatter plot, give it a name and a note, then use it to select them.

    Given user is logged in
    And simple mode is off
    And user opens iris dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    And user adds a scatter plot viewer with:
      | X     | Petal.Length |
      | Y     | Petal.Width  |
      | Color | Species      |
    And scatter plot viewer takes 65% of the space it shares with the grid
    When the video says "An annotation region marks and names a part of a chart"
    # caption: Right-click the plot: Tools > Draw Annotation Region
    And user picks "Tools > Draw Annotation Region" from the context menu of the "view" area of scatter plot viewer
    Then the "region drawing mode" reading of scatter plot viewer should be "true"
    # caption: Drag a box around the points to mark
    When user drags across scatter plot viewer from 4.45, 1.85 to 5.15, 1.35 on its axes
    Then "Formula Lines" dialog should be visible
    # caption: Give it a name and a note
    When user enters "Overlap" into "Title" text input in "Formula Lines" dialog
    And user enters "Versicolor and virginica mix here" into "Description" text area in "Formula Lines" dialog
    And user clicks on OK button in "Formula Lines" dialog
    Then scatter plot viewer should have a 'region "Overlap"' area
    And the "regions shown" reading of scatter plot viewer should be 1
    # caption: Point at it: the name, the note and how many rows it holds
    When user points at 4.82, 1.45 inside the 'region "Overlap"' area of scatter plot viewer
    Then tooltip should contain text "26 rows"
    # caption: A click selects the rows inside
    When user clicks at 4.82, 1.45 inside the 'region "Overlap"' area of scatter plot viewer
    Then 26 rows should be selected
    And some rows where "Species" is "versicolor" should be selected
    And some rows where "Species" is "virginica" should be selected
    And no rows where "Species" is "setosa" should be selected
    # caption: The status bar counts them: 26 selected
    When user hovers over "Selected: 26" text in status bar
    Then no errors should have been logged
