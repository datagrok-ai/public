@journey @viewers @realizes:viewers.scatter-plot @realizes:viewers.histogram @realizes:viewers.line-chart
Feature: Annotation regions
  Regions drawn on a plot and kept in the viewer's look: how one is created from the Tools menu,
  what the viewer reports about it, what the visibility switches do to it, and what the one-
  dimensional viewers persist instead (an axis-locked value band written as a formula region).
  Every claim is a reading or a hit area the viewer publishes — `viewer regions` and `dataframe
  regions` count what the look holds, `regions shown` counts what is drawn, and a `region N` hit
  area exists only while the region is actually on screen, so hiding one removes its area rather
  than leaving the rectangle of the frame that last drew it.
  One journey on demog-1000. PowerPack is installed on this stand, so the Formula Lines dialog
  opens after a region is drawn and each scenario accepts it with OK; the spec's "PowerPack absent"
  section cannot be translated here. Both the axis annotations group and the viewer's annotation
  property group are called "Annotations", so the claim about the count axis is made about
  its Add Line, Add Band and Add Region items. Picking one of them is not translated: the platform
  hides a group's flyout as soon as the pointer leaves the group row, and it takes the item away
  before the click lands, so no honest gesture reaches it. The count axis is
  not claimed about: a right-click there carries the viewer's whole menu, where the viewer's own
  "Annotations" property group has the same name, so the difference cannot be stated honestly at
  the label level. Region hover and click selection are
  the manual checks of `annotation-regions-ui.md` and stay manual.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | xColumnName | AGE    |
      | yColumnName | WEIGHT |
      | lassoTool   | false  |
    Then the "viewer regions" reading of scatter plot viewer should be 0
    And the "regions shown" reading of scatter plot viewer should be 0
    And scatter plot viewer should not have a "region 1" area

  Scenario: Drawing a rectangle adds a region the viewer paints and reports
    When user picks "Tools > Draw Annotation Region" from the context menu of scatter plot viewer
    Then the "region drawing mode" reading of scatter plot viewer should be "true"
    When user drags across the "view" area of scatter plot viewer
    Then the "viewer regions" reading of scatter plot viewer should be 1
    And the "regions shown" reading of scatter plot viewer should be 1
    And the "dataframe regions" reading of scatter plot viewer should be 0
    And scatter plot viewer should have a "region 1" area
    And the "region 1" area of scatter plot viewer should be painted
    And scatter plot viewer should have repainted
    And "annotationRegions" property of scatter plot viewer should contain "area"
    When user clicks OK button in "Formula Lines" dialog
    Then the "viewer regions" reading of scatter plot viewer should be 1
    And the "regions shown" reading of scatter plot viewer should be 1
    And no errors should have been logged

  Scenario: Hiding the viewer's regions takes the region off the plot but not out of the look
    When user takes a snapshot of scatter plot viewer
    And user sets "showViewerAnnotationRegions" property of scatter plot viewer to "false"
    Then the "viewer regions" reading of scatter plot viewer should be 1
    And the "regions shown" reading of scatter plot viewer should be 0
    And scatter plot viewer should not have a "region 1" area
    And scatter plot viewer should have repainted
    When user sets "showViewerAnnotationRegions" property of scatter plot viewer to "true"
    Then the "regions shown" reading of scatter plot viewer should be 1
    And scatter plot viewer should have a "region 1" area
    And no errors should have been logged

  Scenario: Show Annotation Regions in the Tools menu switches both kinds at once
    When user picks "Tools > Show Annotation Regions" from the context menu of the "empty space" area of scatter plot viewer
    Then "showViewerAnnotationRegions" property of scatter plot viewer should be "false"
    And "showDataframeAnnotationRegions" property of scatter plot viewer should be "false"
    And the "regions shown" reading of scatter plot viewer should be 0
    When user picks "Tools > Show Annotation Regions" from the context menu of the "empty space" area of scatter plot viewer
    Then "showViewerAnnotationRegions" property of scatter plot viewer should be "true"
    And the "regions shown" reading of scatter plot viewer should be 1
    And no errors should have been logged

  Scenario: A histogram locks the categorical axis and stores the band as a formula region
    Given user adds a histogram viewer with:
      | valueColumnName | AGE |
    Then the "viewer regions" reading of histogram viewer should be 0
    When user picks "Tools > Draw Annotation Region" from the context menu of histogram viewer
    And user drags across the "view" area of histogram viewer
    Then the "viewer regions" reading of histogram viewer should be 1
    And the "regions shown" reading of histogram viewer should be 1
    And histogram viewer should have a "region 1" area
    And "annotationRegions" property of histogram viewer should contain "formula"
    And "annotationRegions" property of histogram viewer should contain "${AGE}"
    When user clicks OK button in "Formula Lines" dialog
    Then no errors should have been logged

  Scenario: The value axis offers Add Line, Add Band and Add Region
    When user right-clicks on the "x axis" area of histogram viewer
    Then the open menu should list "Annotations > Add Line"
    And the open menu should list "Annotations > Add Band"
    And the open menu should list "Annotations > Add Region"
    When user closes the context menu
    Then no errors should have been logged

  Scenario: A multi-axis line chart offers no Draw Annotation Region
    Given user adds a line chart viewer with:
      | xColumnName  | AGE            |
      | yColumnNames | WEIGHT, HEIGHT |
      | multiAxis    | false          |
    When user opens the context menu of line chart viewer
    Then the open menu should list "Tools > Draw Annotation Region"
    When user closes the context menu
    And user sets "multiAxis" property of line chart viewer to "true"
    And user opens the context menu of line chart viewer
    Then the open menu should not list "Tools > Draw Annotation Region"
    When user closes the context menu
    Then no errors should have been logged
