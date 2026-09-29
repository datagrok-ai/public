@viewers @realizes:viewers.line-chart @realizes:viewers.density-plot
Feature: Annotation regions on the other two-dimensional viewers
  The line chart and the density plot share the scatter plot's region machinery: an area region
  keyed by its X and Y columns is drawn, published as a `region <name>` hit area with its title as
  `region <name> title`, and counted. The line chart aggregates its Y columns per X value (AGE is
  an integer column, so every value repeats) and binds a region to the aggregated column — the
  caption a region drawn on the chart stores, `avg(WEIGHT)`, not the plain `WEIGHT`, which the
  chart does not match. A region belongs to the chart of its Y column, so it stays when a second Y
  column stacks another chart under it; a formula region constant on X spans every chart and its
  title takes the strip above the first chart (inside the first chart's box, above the band's
  own rectangle). On demog-1000. PowerPack is installed on this stand, so the Formula Lines dialog
  opens after a region is drawn and the scenario accepts it with OK.
  The band title after a resize is claimed only after the resize: at creation the title lands a
  few seconds late (the formula-lines feature initialises asynchronously and relays the strip out).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset

  Scenario: A region drawn on a line chart is keyed to the aggregated Y column
    Given user adds a line chart viewer with:
      | xColumnName  | AGE    |
      | yColumnNames | WEIGHT |
    And user resizes line chart viewer to 800 by 500
    Then the "aggregated" reading of line chart viewer should be "true"
    When user picks "Tools > Draw Annotation Region" from the context menu of line chart viewer
    And user drags across the "plot" area of line chart viewer
    Then the "viewer regions" reading of line chart viewer should be 1
    And "annotationRegions" property of line chart viewer should contain "avg(WEIGHT)"
    And line chart viewer should have a "region 1" area
    When user clicks OK button in "Formula Lines" dialog
    Then no errors should have been logged

  Scenario: A line chart draws a titled area region on the chart of its Y column
    Given user adds a line chart viewer with:
      | xColumnName       | AGE    |
      | yColumnNames      | WEIGHT |
      | annotationRegions | [{"type":"area","x":"AGE","y":"avg(WEIGHT)","header":"Adults","area":[[30.5,60],[60.5,60],[60.5,140],[30.5,140]]}] |
    And user resizes line chart viewer to 800 by 500
    Then the "viewer regions" reading of line chart viewer should be 1
    And the "regions shown" reading of line chart viewer should be 1
    And line chart viewer should have a "region Adults" area
    And line chart viewer should have a "region Adults title" area
    And the "region Adults title" area of line chart viewer should lie inside the "region Adults" area
    When user sets "yColumnNames" property of line chart viewer to "WEIGHT, HEIGHT"
    Then the "charts" reading of line chart viewer should be 2
    And line chart viewer should have a "region Adults" area
    And the "region Adults title" area of line chart viewer should lie inside the "chart WEIGHT" area
    And no errors should have been logged

  Scenario: A line chart band constant on X reserves the strip above its charts
    Given user adds a line chart viewer with:
      | xColumnName       | AGE            |
      | yColumnNames      | WEIGHT, HEIGHT |
      | annotationRegions | [{"type":"formula","header":"Adults","formula1":"${AGE} = 30","formula2":"${AGE} = 60"}] |
    Then the "regions shown" reading of line chart viewer should be 1
    And line chart viewer should have a "region Adults" area
    And line chart viewer should have a "region Adults title" area
    And the "title strip top" reading of line chart viewer should be at least 1
    And the "region Adults title" area of line chart viewer should lie above the "region Adults" area
    And the "region Adults title" area of line chart viewer should lie inside the "chart 1" area
    And no errors should have been logged

  Scenario: A resized line chart keeps its band title
    Given user adds a line chart viewer with:
      | xColumnName       | AGE            |
      | yColumnNames      | WEIGHT, HEIGHT |
      | annotationRegions | [{"type":"formula","header":"Adults","formula1":"${AGE} = 30","formula2":"${AGE} = 60"}] |
    Then the "regions shown" reading of line chart viewer should be 1
    And line chart viewer should have a "region Adults" area
    When user resizes line chart viewer to 800 by 500
    Then line chart viewer should have a "region Adults" area
    And line chart viewer should have a "region Adults title" area
    And the "title strip top" reading of line chart viewer should be at least 1
    And no errors should have been logged

  Scenario: A density plot draws a titled area region inside its bins
    Given user adds a density plot viewer with:
      | xColumnName       | AGE    |
      | yColumnName       | WEIGHT |
      | annotationRegions | [{"type":"area","x":"AGE","y":"WEIGHT","header":"Adults","area":[[30.5,60],[60.5,60],[60.5,140],[30.5,140]]}] |
    And user resizes density plot viewer to 800 by 500
    Then the "viewer regions" reading of density plot viewer should be 1
    And the "regions shown" reading of density plot viewer should be 1
    And density plot viewer should have a "region Adults" area
    And density plot viewer should have a "region Adults title" area
    And the "region Adults title" area of density plot viewer should lie inside the "region Adults" area
    And the "region Adults" area of density plot viewer should lie inside the "view" area
    When user sets "showViewerAnnotationRegions" property of density plot viewer to "false"
    Then density plot viewer should not have a "region Adults" area
    And density plot viewer should not have a "region Adults title" area
    When user sets "showViewerAnnotationRegions" property of density plot viewer to "true"
    Then density plot viewer should have a "region Adults" area
    And density plot viewer should have a "region Adults title" area
    And no errors should have been logged
