@viewers @realizes:viewers.scatter-plot @realizes:powerpack.dialogs.formula-lines
Feature: A lasso drawn on the scatter plot becomes an area region
  `annotation-regions-drawing.feature` draws a lasso region on the density plot and a rectangle on
  the scatter plot; this is the scatter plot's own lasso. With the Lasso Tool on, Tools > Draw
  Annotation Region turns drawing on, and the lasso drawn over the plot becomes an area region
  keyed to the X and Y columns and opens PowerPack's Formula Lines dialog. On demog-1000 with small
  markers; the region is removed through the dialog's Delete at the end. The feature claims only
  that an area region keyed to WEIGHT and HEIGHT is drawn; that it keeps the free-form shape the
  lasso went round has no claim yet (MISSING.md).

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | xColumnName       | WEIGHT |
      | yColumnName       | HEIGHT |
      | lassoTool         | true   |
      | markerDefaultSize | 2      |
    And user resizes scatter plot viewer to 800 by 500

  Scenario: The lasso becomes a titled area region, removed again from the dialog
    Then the "viewer regions" reading of scatter plot viewer should be 0
    When user picks "Tools > Draw Annotation Region" from the context menu of scatter plot viewer
    Then the "region drawing mode" reading of scatter plot viewer should be "true"
    When user drags a lasso over the "view" area of scatter plot viewer
    Then the "region drawing mode" reading of scatter plot viewer should be "false"
    And the "viewer regions" reading of scatter plot viewer should be 1
    And scatter plot viewer should have a "region 1" area
    And "annotationRegions" property of scatter plot viewer should contain "\"type\":\"area\""
    And "annotationRegions" property of scatter plot viewer should contain "\"x\":\"WEIGHT\",\"y\":\"HEIGHT\""
    And "Formula Lines" dialog should be visible
    And editor of "X column" input in "Formula Lines" dialog should have text "WEIGHT"
    And editor of "Y column" input in "Formula Lines" dialog should have text "HEIGHT"
    When user enters "Lasso A" into Title input in "Formula Lines" dialog
    And user clicks OK button in "Formula Lines" dialog
    Then the "Formula Lines" dialog should close
    And scatter plot viewer should have a "region Lasso A" area
    When user picks "Tools > Formula Lines..." from the context menu of the "bottom edge of view" area of scatter plot viewer
    Then "Formula Lines" dialog should be visible
    When user clicks on Delete button in "Formula Lines" dialog
    And user clicks OK button in "Formula Lines" dialog
    Then the "viewer regions" reading of scatter plot viewer should be 0
    And scatter plot viewer should not have a "region Lasso A" area
    And no errors should have been logged
