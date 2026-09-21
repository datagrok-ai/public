@viewers @realizes:viewers.scatter-plot @realizes:viewers.density-plot @realizes:viewers.histogram @realizes:viewers.bar-chart
Feature: Hovering and clicking regions on every viewer
  What regions do under the pointer on the viewers other than the scatter plot journey covers:
  pointing at a region counts it and shows its title, pointing at an overlap counts both, a click
  selects the region's rows (the rows both regions share, at an overlap) and Control-click on a
  fully selected region removes them. A region reacts only where nothing else lies under the
  pointer — a marker, a bin, a bar, the histogram's bins slider at the top of its plot — so every
  gesture aims at a part of the region the data leaves empty: the left end of the "Tall" band,
  the left end of the overlap of "Tall" and "Heavy" (tall and heavy people are rare, and the
  plot's Color and Size selectors sit over its top-right corner), the top-left corner of the
  histogram's "Older" band (short bars far from the slider), the middle of the bar chart's band
  between two bars. On demog-1000: 159 rows have a HEIGHT of 180 to 200, 21 of them
  weigh 110 to 170, 171 are aged 60 to 90, and the bar chart's "Mid" band at avg(AGE) 45.7 to
  47.1 holds the two bars whose average lies in it — Black (46.7) and Other (47.1), 89 rows: a
  band on an aggregated axis selects by the bar, not by the row's own value. From
  `annotation-regions-interaction.md`.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset

  Scenario: Two overlapping bands on the scatter plot
    Given user adds a scatter plot viewer with:
      | xColumnName       | WEIGHT |
      | yColumnName       | HEIGHT |
      | markerDefaultSize | 2      |
      | annotationRegions | [{"type":"formula","header":"Tall","formula1":"${HEIGHT} = 180","formula2":"${HEIGHT} = 200"},{"type":"formula","header":"Heavy","formula1":"${WEIGHT} = 110","formula2":"${WEIGHT} = 170"}] |
    And user resizes scatter plot viewer to 800 by 500
    Then the "regions shown" reading of scatter plot viewer should be 2
    And the "region titles shown" reading of scatter plot viewer should be 2
    When user hovers over the "left edge of region Tall" area of scatter plot viewer
    Then the "regions hovered" reading of scatter plot viewer should be 1
    And tooltip should contain text "Tall"
    And tooltip should contain text "159 rows"
    When user hovers over the "left edge of overlap of region Tall and region Heavy" area of scatter plot viewer
    Then the "regions hovered" reading of scatter plot viewer should be 2
    And tooltip should contain text "Tall"
    And tooltip should contain text "Heavy"
    And tooltip should contain text "21 rows"
    When user moves the pointer away from scatter plot viewer
    Then the "regions hovered" reading of scatter plot viewer should be 0
    When user hovers over the "left edge of region Tall" area of scatter plot viewer
    And user clicks on the "left edge of region Tall" area of scatter plot viewer
    Then only rows where "HEIGHT" is between 180 and 200 should be selected
    When user clicks on the "left edge of region Tall" area of scatter plot viewer holding Control
    Then no rows should be selected
    When user hovers over the "left edge of overlap of region Tall and region Heavy" area of scatter plot viewer
    And user clicks on the "left edge of overlap of region Tall and region Heavy" area of scatter plot viewer
    Then 21 rows should be selected
    And every selected row should pass the filter
    When user clears the row selection
    And user moves the pointer away from scatter plot viewer
    Then no errors should have been logged

  Scenario: Two overlapping bands on the density plot
    Given user adds a density plot viewer with:
      | xColumnName       | WEIGHT |
      | yColumnName       | HEIGHT |
      | annotationRegions | [{"type":"formula","header":"Tall","formula1":"${HEIGHT} = 180","formula2":"${HEIGHT} = 200"},{"type":"formula","header":"Heavy","formula1":"${WEIGHT} = 110","formula2":"${WEIGHT} = 170"}] |
    And user resizes density plot viewer to 800 by 500
    Then the "regions shown" reading of density plot viewer should be 2
    When user hovers over the "left edge of region Tall" area of density plot viewer
    Then the "regions hovered" reading of density plot viewer should be 1
    When user hovers over the "left edge of overlap of region Tall and region Heavy" area of density plot viewer
    Then the "regions hovered" reading of density plot viewer should be 2
    When user clicks on the "left edge of overlap of region Tall and region Heavy" area of density plot viewer
    Then 21 rows should be selected
    When user hovers over the "left edge of region Tall" area of density plot viewer
    And user clicks on the "left edge of region Tall" area of density plot viewer
    Then only rows where "HEIGHT" is between 180 and 200 should be selected
    When user clicks on the "left edge of region Tall" area of density plot viewer holding Control
    Then no rows should be selected
    When user moves the pointer away from density plot viewer
    Then no errors should have been logged

  Scenario: A value band on the histogram reacts where no bar is under it
    Given user adds a histogram viewer with:
      | valueColumnName   | AGE |
      | annotationRegions | [{"type":"formula","header":"Older","formula1":"${AGE} = 60","formula2":"${AGE} = 90"}] |
    And user resizes histogram viewer to 800 by 500
    Then the "regions shown" reading of histogram viewer should be 1
    When user hovers over the "top left corner of region Older" area of histogram viewer
    Then the "regions hovered" reading of histogram viewer should be 1
    And tooltip should contain text "Older"
    When user clicks on the "top left corner of region Older" area of histogram viewer
    Then only rows where "AGE" is between 60 and 90 should be selected
    When user clicks on the "top left corner of region Older" area of histogram viewer holding Control
    Then no rows should be selected
    When user moves the pointer away from histogram viewer
    Then the "regions hovered" reading of histogram viewer should be 0
    And no errors should have been logged

  Scenario: A value band on the bar chart reacts between the bars
    Given user adds a bar chart viewer with:
      | splitColumnName   | RACE     |
      | valueColumnName   | AGE      |
      | valueAggrType     | avg      |
      | orientation       | vertical |
      | annotationRegions | [{"type":"formula","header":"Mid","formula1":"${AGE} = 45.7","formula2":"${AGE} = 47.1"}] |
    And user resizes bar chart viewer to 800 by 500
    Then the "regions shown" reading of bar chart viewer should be 1
    And bar chart viewer should have a "region Mid" area
    And the "region Mid" and "view" areas of bar chart viewer should be the same width
    When user hovers over the "region Mid" area of bar chart viewer
    Then the "regions hovered" reading of bar chart viewer should be 1
    And tooltip should contain text "Mid"
    When user clicks on the "region Mid" area of bar chart viewer
    Then only rows where "RACE" is one of "Black, Other" should be selected
    When user clicks on the "region Mid" area of bar chart viewer holding Control
    Then no rows should be selected
    When user moves the pointer away from bar chart viewer
    Then no errors should have been logged

  Scenario: A value band on the box plot reacts between the boxes
    Given user adds a box plot viewer with:
      | category1ColumnName | RACE |
      | valueColumnName     | AGE  |
      | annotationRegions   | [{"type":"formula","header":"Middle age","formula1":"${AGE} = 36.0","formula2":"${AGE} = 56.0"}] |
    And user resizes box plot viewer to 800 by 500
    Then the "regions shown" reading of box plot viewer should be 1
    And box plot viewer should have a "region Middle age" area
    And the "region Middle age" and "view" areas of box plot viewer should be the same width
    When user hovers over the "left edge of region Middle age" area of box plot viewer
    Then the "regions hovered" reading of box plot viewer should be 1
    And tooltip should contain text "Middle age"
    And tooltip should contain text "517 rows"
    When user clicks on the "left edge of region Middle age" area of box plot viewer
    Then only rows where "AGE" is between 36 and 56 should be selected
    When user clicks on the "left edge of region Middle age" area of box plot viewer holding Control
    Then no rows should be selected
    When user moves the pointer away from box plot viewer
    Then the "regions hovered" reading of box plot viewer should be 0
    And no errors should have been logged

  Scenario: Two overlapping bands on the line chart
    Given user adds a line chart viewer with:
      | xColumnName       | AGE    |
      | yColumnNames      | HEIGHT |
      | annotationRegions | [{"type":"formula","header":"Medium height","formula1":"${avg(HEIGHT)} = 165","formula2":"${avg(HEIGHT)} = 175"},{"type":"formula","header":"Older","formula1":"${AGE} = 60","formula2":"${AGE} = 90"}] |
    And user resizes line chart viewer to 800 by 500
    Then the "regions shown" reading of line chart viewer should be 2
    When user hovers over the "left edge of region Medium height" area of line chart viewer
    Then the "regions hovered" reading of line chart viewer should be 1
    And tooltip should contain text "Medium height"
    And tooltip should contain text "822 rows"
    When user hovers over the "left edge of overlap of region Medium height and region Older" area of line chart viewer
    Then the "regions hovered" reading of line chart viewer should be 2
    And tooltip should contain text "Medium height"
    And tooltip should contain text "Older"
    And tooltip should contain text "102 rows"
    When user clicks on the "left edge of region Medium height" area of line chart viewer
    Then 822 rows should be selected
    When user clicks on the "left edge of region Medium height" area of line chart viewer holding Control
    Then no rows should be selected
    When user clicks on the "left edge of overlap of region Medium height and region Older" area of line chart viewer
    Then 102 rows should be selected
    When user clears the row selection
    And user moves the pointer away from line chart viewer
    Then the "regions hovered" reading of line chart viewer should be 0
    And no errors should have been logged
