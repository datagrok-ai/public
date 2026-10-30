@journey @viewers @realizes:viewers.scatter-plot3d
Feature: 3D scatter plot interaction
  What a person does on the plot that scatter-plot-3d.feature leaves out. Hovering a point shows
  its row in a tooltip (the `point` area is the drawn marker nearest to the camera: at the home
  camera on demog-1000 that is row 779, X0273T44004500022, AGE 52, WEIGHT 163 near the top of the Z axis;
  the plot holds its camera on a test page). Ctrl-click toggles that row's selection and
  Ctrl+Shift-click takes it out, while neither moves the current row; the selection is claimed on
  the table, not on the frame, since the mouse-over box drawn round the hovered marker hides its
  new color. The
  context menu's Filter > Show Filtered Out Points brings the filtered-away rows back and a second
  pick takes them away. A double-click on empty space brings a rotated camera home. Row Source
  Selected draws only the selected rows that have all three coordinates (499 of the 553 women,
  373 of the 447 men). A string column on an axis draws: the frame changes and every row stays
  on the plot (public#3594); where its categories land is not read (MISSING.md). Two regressions on their own tables close the journey: a boolean label column on
  the 5850 rows of demog draws (GROK-19048, core#129: it used to crash the tab), and a SMILES
  column on an axis does not freeze the plot (public#398; needs Chem for the column to be a
  molecule column). The current-row highlight and Save as PNG are in MISSING.md.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a 3d scatter plot viewer
    Then the "rows shown" reading of 3d scatter plot viewer should be 872
    And "USUBJID" of the current row should be "X0273T21000300003"

  Scenario: Hovering a point shows its row in a tooltip
    When user hovers over the "point" area of 3d scatter plot viewer
    Then the mouse-over row of the table should be 779
    And exactly one tooltip should be shown
    And the tooltip should show "USUBJID" as "X0273T44004500022"
    And the tooltip should show "AGE" as "52"
    And no errors should have been logged

  Scenario: Ctrl-click toggles a point's selection and Ctrl+Shift-click removes it
    Given user clears the row selection
    When user clicks on the "point" area of 3d scatter plot viewer holding Control
    Then 1 row should be selected
    And only rows where "USUBJID" is "X0273T44004500022" should be selected
    And "USUBJID" of the current row should be "X0273T21000300003"
    When user clicks on the "point" area of 3d scatter plot viewer holding Control
    Then no rows should be selected
    When user clicks on the "point" area of 3d scatter plot viewer holding Shift
    Then only rows where "USUBJID" is "X0273T44004500022" should be selected
    When user clicks on the "point" area of 3d scatter plot viewer holding Control+Shift
    Then no rows should be selected
    When user clicks on the "point" area of 3d scatter plot viewer holding Control+Shift
    Then no rows should be selected
    And "USUBJID" of the current row should be "X0273T21000300003"
    And no errors should have been logged

  Scenario: Filter > Show Filtered Out Points in the context menu brings the filtered-out points back
    When user filters rows where "SEX" is "F"
    Then 553 rows should pass the filter
    And the "rows shown" reading of 3d scatter plot viewer should be 499
    When user picks "Filter > Show Filtered Out Points" from the context menu of 3d scatter plot viewer
    And user closes the context menu
    Then "Show Filtered Out Points" property of 3d scatter plot viewer should be "true"
    And the "rows shown" reading of 3d scatter plot viewer should be 872
    And the "scene signature" reading of 3d scatter plot viewer should differ from before
    When user picks "Filter > Show Filtered Out Points" from the context menu of 3d scatter plot viewer
    And user closes the context menu
    Then the "rows shown" reading of 3d scatter plot viewer should be 499
    When user resets the filter
    Then the "rows shown" reading of 3d scatter plot viewer should be 872
    And no errors should have been logged

  Scenario: A double-click on empty space brings the camera home
    When user drags across the "view" area of 3d scatter plot viewer
    Then the "camera x" reading of 3d scatter plot viewer should differ from before
    When user double-clicks on empty plot space of 3d scatter plot viewer
    Then the "camera x" reading of 3d scatter plot viewer should be 0
    And the "camera y" reading of 3d scatter plot viewer should be 0
    And the "camera distance" reading of 3d scatter plot viewer should be 4
    And the "scene signature" reading of 3d scatter plot viewer should differ from before
    And no errors should have been logged

  Scenario: Row Source Selected draws only the selected rows
    Given user clears the row selection
    When user selects rows where "SEX" is "F"
    Then 553 rows should be selected
    When user sets "Row Source" property of 3d scatter plot viewer to "Selected"
    Then the "rows shown" reading of 3d scatter plot viewer should be 499
    When user selects rows where "SEX" is "M"
    Then 447 rows should be selected
    And the "rows shown" reading of 3d scatter plot viewer should be 373
    When user clears the row selection
    Then the "rows shown" reading of 3d scatter plot viewer should be 0
    When user sets "Row Source" property of 3d scatter plot viewer to "Filtered"
    Then the "rows shown" reading of 3d scatter plot viewer should be 872
    And no errors should have been logged

  Scenario: A string column on an axis draws every row
    When user takes a snapshot of 3d scatter plot viewer
    And user picks "RACE" in the "X" column selector of 3d scatter plot viewer
    Then X column input in 3d scatter plot viewer should contain text "RACE"
    And the "rows shown" reading of 3d scatter plot viewer should be 872
    And the "scene signature" reading of 3d scatter plot viewer should differ from before
    When user picks "AGE" in the "X" column selector of 3d scatter plot viewer
    Then X column input in 3d scatter plot viewer should contain text "AGE"
    And no errors should have been logged

  Scenario: A boolean label column on 5850 rows draws without freezing (GROK-19048)
    Given user opens demog dataset
    And user adds a 3d scatter plot viewer
    Then the "rows shown" reading of 3d scatter plot viewer should be 5097
    When user sets "Label" property of 3d scatter plot viewer to "CONTROL"
    Then the "scene signature" reading of 3d scatter plot viewer should differ from before
    And no errors should have been logged

  Scenario: A SMILES column on an axis does not freeze the plot (public#398)
    Given the "Chem" package is installed
    And user opens smiles dataset
    And user adds a 3d scatter plot viewer
    Then "canonical_smiles" column should have semantic type "Molecule"
    When user takes a snapshot of 3d scatter plot viewer
    And user picks "canonical_smiles" in the "X" column selector of 3d scatter plot viewer
    Then X column input in 3d scatter plot viewer should contain text "canonical_smiles"
    And the "rows shown" reading of 3d scatter plot viewer should be 1000
    And the "scene signature" reading of 3d scatter plot viewer should differ from before
    When user takes a snapshot of 3d scatter plot viewer
    And user drags across the "view" area of 3d scatter plot viewer
    Then the "camera x" reading of 3d scatter plot viewer should differ from before
    And the "scene signature" reading of 3d scatter plot viewer should differ from before
    And no errors should have been logged
