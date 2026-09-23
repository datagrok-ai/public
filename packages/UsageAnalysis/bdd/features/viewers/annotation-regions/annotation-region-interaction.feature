@journey @viewers @realizes:viewers.scatter-plot
Feature: Annotation region interaction
  What a region does once it is on the plot: hovering counts it, clicking selects the rows it
  holds (the intersection where regions overlap, never past the filter), the modifier keys compose
  the selection the way every viewer's click does, a dataframe region drawn from the table's tag
  shows and hides on its own switch, a hidden region is kept in the look but not drawn, the axes
  may be swapped but not replaced, a column rename follows into the regions, the regions survive a
  layout round trip (that one in `annotation-regions-persistence.feature`, on a fresh view), and
  Escape leaves the drawing mode.
  One journey on demog-1000 (AGE 18..89, WEIGHT 41.6..165, no blanks) with two nested regions:
  "Outer" spans AGE 20.5..60.5 and "Inner" AGE 30.5..38.5, both the whole WEIGHT range, so the
  centre of Inner lies in both and the centre of Outer in Outer alone. AGE is an integer column, so
  the half-value bounds put no row on an edge.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | xColumnName       | AGE    |
      | yColumnName       | WEIGHT |
      | lassoTool         | false  |
      | annotationRegions | [{"type":"area","x":"AGE","y":"WEIGHT","header":"Outer","area":[[20.5,30],[60.5,30],[60.5,180],[20.5,180]]},{"type":"area","x":"AGE","y":"WEIGHT","header":"Inner","area":[[30.5,30],[38.5,30],[38.5,180],[30.5,180]]}] |
    And user resizes scatter plot viewer to 800 by 500
    Then the "viewer regions" reading of scatter plot viewer should be 2
    And the "regions shown" reading of scatter plot viewer should be 2
    And scatter plot viewer should have a "region Outer" area
    And scatter plot viewer should have a "region Inner" area

  Scenario: Hovering counts the regions under the pointer
    When user hovers over the "region Outer" area of scatter plot viewer
    Then the "regions hovered" reading of scatter plot viewer should be 1
    When user hovers over the "region Inner" area of scatter plot viewer
    Then the "regions hovered" reading of scatter plot viewer should be 2
    When user moves the pointer away from scatter plot viewer
    Then the "regions hovered" reading of scatter plot viewer should be 0
    And no errors should have been logged

  Scenario: A click selects the rows the region holds
    When user clicks on the "region Outer" area of scatter plot viewer
    Then only rows where "AGE" is between 21 and 60 should be selected
    And no errors should have been logged

  Scenario: A click where regions overlap selects the rows they share
    When user clicks on the "region Inner" area of scatter plot viewer
    Then only rows where "AGE" is between 31 and 38 should be selected
    And no errors should have been logged

  Scenario: Control drops a selected region's rows and Shift adds a region's rows
    When user clicks on the "region Inner" area of scatter plot viewer
    Then only rows where "AGE" is between 31 and 38 should be selected
    When user clicks on the "region Inner" area of scatter plot viewer holding Control
    Then no rows should be selected
    When user clicks on the "region Outer" area of scatter plot viewer
    And user clicks on the "region Inner" area of scatter plot viewer holding Shift
    Then only rows where "AGE" is between 21 and 60 should be selected
    When user clears the row selection
    Then no errors should have been logged

  Scenario: A click selects nothing the filter has taken out
    When user moves the pointer away from scatter plot viewer
    And user filters rows where "SEX" is "F"
    Then scatter plot viewer should show fewer rows than before
    When user hovers over the "region Outer" area of scatter plot viewer
    Then the "regions hovered" reading of scatter plot viewer should be 1
    When user clicks on the "region Outer" area of scatter plot viewer
    Then some rows should be selected
    And every selected row should pass the filter
    And no rows where "SEX" is "M" should be selected
    When user resets the filter
    And user clears the row selection
    Then no errors should have been logged

  Scenario: A dataframe region is drawn next to the viewer's and hidden on its own switch
    When user sets the ".annotation-regions" tag of the table to:
      """
      [{"type":"area","x":"AGE","y":"WEIGHT","header":"Shared","area":[[65.5,30],[85.5,30],[85.5,180],[65.5,180]]}]
      """
    Then the "dataframe regions" reading of scatter plot viewer should be 1
    And the "regions shown" reading of scatter plot viewer should be 3
    And scatter plot viewer should have a "region Shared" area
    When user sets "showDataframeAnnotationRegions" property of scatter plot viewer to "false"
    Then the "regions shown" reading of scatter plot viewer should be 2
    And scatter plot viewer should not have a "region Shared" area
    And scatter plot viewer should have a "region Outer" area
    When user sets "showDataframeAnnotationRegions" property of scatter plot viewer to "true"
    Then scatter plot viewer should have a "region Shared" area
    When user sets the ".annotation-regions" tag of the table to ""
    Then the "dataframe regions" reading of scatter plot viewer should be 0
    And the "regions shown" reading of scatter plot viewer should be 2
    And scatter plot viewer should not have a "region Shared" area
    And no errors should have been logged

  Scenario: A hidden region is kept in the look but not drawn
    When user sets properties of scatter plot viewer:
      | annotationRegions | [{"type":"area","x":"AGE","y":"WEIGHT","header":"Outer","area":[[20.5,30],[60.5,30],[60.5,180],[20.5,180]]},{"type":"area","x":"AGE","y":"WEIGHT","header":"Inner","hidden":true,"area":[[30.5,30],[38.5,30],[38.5,180],[30.5,180]]}] |
    Then the "viewer regions" reading of scatter plot viewer should be 2
    And the "regions shown" reading of scatter plot viewer should be 1
    And scatter plot viewer should have a "region Outer" area
    And scatter plot viewer should not have a "region Inner" area
    When user sets properties of scatter plot viewer:
      | annotationRegions | [{"type":"area","x":"AGE","y":"WEIGHT","header":"Outer","area":[[20.5,30],[60.5,30],[60.5,180],[20.5,180]]},{"type":"area","x":"AGE","y":"WEIGHT","header":"Inner","area":[[30.5,30],[38.5,30],[38.5,180],[30.5,180]]}] |
    Then the "regions shown" reading of scatter plot viewer should be 2
    And no errors should have been logged

  Scenario: Swapped axes keep the regions and another column drops them
    When user sets properties of scatter plot viewer:
      | xColumnName | WEIGHT |
      | yColumnName | AGE    |
    Then scatter plot viewer should have a "region Outer" area
    And scatter plot viewer should have a "region Inner" area
    When user sets "xColumnName" property of scatter plot viewer to "HEIGHT"
    Then scatter plot viewer should not have a "region Outer" area
    And scatter plot viewer should not have a "region Inner" area
    And the "region titles shown" reading of scatter plot viewer should be 0
    And the "viewer regions" reading of scatter plot viewer should be 2
    When user sets properties of scatter plot viewer:
      | xColumnName | AGE    |
      | yColumnName | WEIGHT |
    Then scatter plot viewer should have a "region Outer" area
    And scatter plot viewer should have a "region Inner" area
    And no errors should have been logged

  Scenario: A column rename follows into the regions
    When user renames "AGE" column to "AGE (years)"
    Then "annotationRegions" property of scatter plot viewer should contain "\"x\":\"AGE (years)\""
    And scatter plot viewer should have a "region Outer" area
    And scatter plot viewer should have a "region Inner" area
    When user renames "AGE (years)" column to "AGE"
    Then "annotationRegions" property of scatter plot viewer should contain "\"x\":\"AGE\""
    And scatter plot viewer should have a "region Outer" area
    And scatter plot viewer should have a "region Inner" area
    And no errors should have been logged

  Scenario: Escape leaves the drawing mode with the regions untouched
    When user sets "showViewerAnnotationRegions" property of scatter plot viewer to "false"
    And user moves the pointer away from scatter plot viewer
    And user picks "Tools > Draw Annotation Region" from the context menu of scatter plot viewer
    Then the "region drawing mode" reading of scatter plot viewer should be "true"
    When user presses Escape
    Then the "region drawing mode" reading of scatter plot viewer should be "false"
    And the "viewer regions" reading of scatter plot viewer should be 2
    When user sets "showViewerAnnotationRegions" property of scatter plot viewer to "true"
    Then scatter plot viewer should have a "region Outer" area
    And scatter plot viewer should have a "region Inner" area
    And no errors should have been logged
