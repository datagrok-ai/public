@journey @viewers @realizes:viewers.scatter-plot @realizes:viewers.box-plot
Feature: Annotation region titles
  Where a region's title is drawn and what must not move it. A title lives either in the data, at
  an anchor the viewer keeps in world coordinates, or — for a band that is constant on an axis
  column — in a strip the layout reserves above (X band) or to the right of (Y band) the plot, and
  only when the title fits the band; a title that does not fit is meant to render in the data and
  reserve nothing (GROK-20388). The viewer publishes every drawn title as a `region <name> title`
  hit area, the strips as the `title strip top` / `title strip right` readings, and the count as
  `region titles shown`, so a claim about placement is a rectangle, not a screenshot.
  Two defects are regression-guarded here: clicking a title relaid it out (GROK-20157), and a
  size or style change made it jump (GROK-20158), and a band title the strip cannot take — too wide
  for it, or Auto Layout off — was not drawn at all: its anchor lay level with a vertex of the
  band's sampled edge, which the polygon's ray cast counted twice and put outside the band.
  One journey on demog-1000: AGE runs 18..89 and WEIGHT 41.6..165, with no blanks in either; the
  area region spans the whole WEIGHT range, so a click on it selects exactly the rows of its AGE band.
  The scatter plot is held at 800 by 500 so the band widths the strip claims depend on are known,
  and its markers are small so a gesture at a title lands on the region, not on a marker (a marker
  under the pointer takes precedence over the regions).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | xColumnName       | AGE    |
      | yColumnName       | WEIGHT |
      | lassoTool         | false  |
      | markerDefaultSize | 2      |
      | annotationRegions | [{"type":"area","x":"AGE","y":"WEIGHT","header":"Adults","area":[[30.5,30],[60.5,30],[60.5,180],[30.5,180]]}] |
    And user resizes scatter plot viewer to 800 by 500
    Then the "regions shown" reading of scatter plot viewer should be 1
    And the "region titles shown" reading of scatter plot viewer should be 1

  Scenario: An in-data title sits inside its region and reserves no strip
    Then scatter plot viewer should have a "region Adults" area
    And scatter plot viewer should have a "region Adults title" area
    And the "region Adults title" area of scatter plot viewer should lie inside the "region Adults" area
    And the "title strip top" reading of scatter plot viewer should be 0
    And the "title strip right" reading of scatter plot viewer should be 0
    And no errors should have been logged

  Scenario: Clicking the title selects the region's rows and leaves the title where it was
    When user remembers the place of the "region Adults title" area of scatter plot viewer
    And user hovers over the "region Adults title" area of scatter plot viewer
    Then the "regions hovered" reading of scatter plot viewer should be 1
    When user clicks on the "region Adults title" area of scatter plot viewer
    Then only rows where "AGE" is between 31 and 60 should be selected
    And the "region Adults title" area of scatter plot viewer should be placed as remembered
    When user clears the row selection
    And user moves the pointer away from scatter plot viewer
    Then the "region Adults title" area of scatter plot viewer should be placed as remembered
    And no errors should have been logged

  Scenario: Growing the viewer and restoring it puts the title back where it was
    When user remembers the place of the "region Adults title" area of scatter plot viewer
    And user resizes scatter plot viewer to 1000 by 600
    Then the "region Adults title" area of scatter plot viewer should lie inside the "region Adults" area
    When user resizes scatter plot viewer to 800 by 500
    Then the "region Adults title" area of scatter plot viewer should be placed as remembered
    And no errors should have been logged

  Scenario: Zooming keeps the title inside the region
    When user scrolls the mouse wheel up over the "region Adults" area of scatter plot viewer
    Then scatter plot viewer should show a narrower value range than before
    And the "region Adults title" area of scatter plot viewer should lie inside the "region Adults" area
    When user scrolls the mouse wheel down over the "region Adults" area of scatter plot viewer
    Then scatter plot viewer should show a wider value range than before
    And the "region Adults title" area of scatter plot viewer should lie inside the "region Adults" area
    And no errors should have been logged

  Scenario: A larger annotation font draws a taller title
    When user sets "annotationFont" property of scatter plot viewer to "normal normal 20px \"Roboto\""
    Then the "region Adults title" area of scatter plot viewer should be taller than before
    And the "region Adults title" area of scatter plot viewer should lie inside the "region Adults" area
    When user sets "annotationFont" property of scatter plot viewer to "normal normal 10px \"Roboto\""
    Then the "region Adults title" area of scatter plot viewer should be shorter than before
    And no errors should have been logged

  Scenario: A band on the X column whose title fits takes the strip above the plot
    When user sets properties of scatter plot viewer:
      | annotationRegions | [{"type":"formula","header":"Adults","formula1":"${AGE} = 30","formula2":"${AGE} = 60"}] |
    Then the "regions shown" reading of scatter plot viewer should be 1
    And the "region titles shown" reading of scatter plot viewer should be 1
    And the "title strip top" reading of scatter plot viewer should be at least 1
    And the "title strip right" reading of scatter plot viewer should be 0
    And the "region Adults title" area of scatter plot viewer should lie above the "view" area
    And no errors should have been logged

  Scenario: A band on the Y column takes the strip to the right of the plot
    When user sets properties of scatter plot viewer:
      | annotationRegions | [{"type":"formula","header":"Middle","formula1":"${WEIGHT} = 80","formula2":"${WEIGHT} = 120"}] |
    Then the "regions shown" reading of scatter plot viewer should be 1
    And the "title strip right" reading of scatter plot viewer should be at least 1
    And the "title strip top" reading of scatter plot viewer should be 0
    And the "region Middle title" area of scatter plot viewer should lie to the right of the "view" area
    And no errors should have been logged

  Scenario: Region titles stay drawn after the formula lines go empty
    When user sets properties of scatter plot viewer:
      | annotationRegions | [{"type":"area","x":"AGE","y":"WEIGHT","header":"Adults","area":[[30.5,30],[60.5,30],[60.5,180],[30.5,180]]}] |
      | formulaLines      | [{"type":"band","title":"Young","formula":"${AGE} in (18, 28)","orientation":"Vertical","column2":"WEIGHT"}] |
    Then the "formula lines" reading of scatter plot viewer should be 1
    And the "title strip top" reading of scatter plot viewer should be at least 1
    And the "region titles shown" reading of scatter plot viewer should be 1
    When user sets "formulaLines" property of scatter plot viewer to "[]"
    Then the "formula lines" reading of scatter plot viewer should be 0
    And the "title strip top" reading of scatter plot viewer should be 0
    And the "region titles shown" reading of scatter plot viewer should be 1
    And scatter plot viewer should have a "region Adults title" area
    And the "region Adults title" area of scatter plot viewer should lie inside the "region Adults" area
    And no errors should have been logged

  Scenario: A box plot keeps a band title in place through a plot style and a bin count change
    Given user adds a box plot viewer with:
      | category1ColumnName | RACE |
      | valueColumnName     | AGE  |
      | annotationRegions   | [{"type":"formula","header":"Prime","formula1":"${AGE} = 30","formula2":"${AGE} = 50"}] |
    Then the "regions shown" reading of box plot viewer should be 1
    And box plot viewer should have a "region Prime" area
    And box plot viewer should have a "region Prime title" area
    And the "title strip right" reading of box plot viewer should be at least 1
    When user remembers the place of the "region Prime title" area of box plot viewer
    And user sets "plotStyle" property of box plot viewer to "violin"
    Then the "region Prime title" area of box plot viewer should be placed as remembered
    When user sets "bins" property of box plot viewer to "20"
    Then the "region Prime title" area of box plot viewer should be placed as remembered
    When user sets "plotStyle" property of box plot viewer to "box"
    Then the "region Prime title" area of box plot viewer should be placed as remembered
    And no errors should have been logged

  Scenario: A band title that does not fit renders in the data and reserves nothing
    When user sets properties of scatter plot viewer:
      | annotationRegions | [{"type":"formula","header":"Adults of the study population","formula1":"${AGE} = 40","formula2":"${AGE} = 45"}] |
    Then the "regions shown" reading of scatter plot viewer should be 1
    And the "region titles shown" reading of scatter plot viewer should be 1
    And the "title strip top" reading of scatter plot viewer should be 0
    And the "region Adults of the study population title" area of scatter plot viewer should lie inside the "view" area
    And no errors should have been logged

  Scenario: Auto Layout off drops the strip and the title moves into the data
    When user sets properties of scatter plot viewer:
      | annotationRegions | [{"type":"formula","header":"Adults","formula1":"${AGE} = 30","formula2":"${AGE} = 60"}] |
      | autoLayout        | true |
    Then the "title strip top" reading of scatter plot viewer should be at least 1
    When user sets "autoLayout" property of scatter plot viewer to "false"
    Then the "title strip top" reading of scatter plot viewer should be 0
    And the "region Adults title" area of scatter plot viewer should lie inside the "view" area
    When user sets "autoLayout" property of scatter plot viewer to "true"
    Then the "title strip top" reading of scatter plot viewer should be at least 1
    And the "region Adults title" area of scatter plot viewer should lie above the "view" area
    And no errors should have been logged
