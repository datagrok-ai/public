@viewers @realizes:viewers.scatter-plot @realizes:viewers.density-plot @realizes:viewers.bar-chart @realizes:powerpack.dialogs.formula-lines
Feature: Drawing and editing annotation regions
  How a region is drawn on a viewer and changed afterwards. On a two-axis viewer the drawn
  rectangle (or lasso polygon) becomes an area region keyed to the X and Y columns; on a one-axis
  viewer the rectangle is locked to the value axis and saved as a two-formula region. Every new
  region opens PowerPack's Formula Lines dialog, whose Title and Description reach the region's
  title and tooltip, and a region is edited again from its own context menu. On demog-1000, with
  small markers so a gesture at a region lands on the region; from `annotation-regions.md`
  scenarios 2 to 4. Rows count what the viewer holds through the `viewer regions` reading and
  the `region <title>` hit areas; the dialog is read through its inputs.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset

  Scenario: A rectangle drawn on the scatter plot becomes a titled area region, edited from its menu
    Given user adds a scatter plot viewer with:
      | xColumnName       | WEIGHT |
      | yColumnName       | HEIGHT |
      | lassoTool         | false  |
      | markerDefaultSize | 2      |
    And user resizes scatter plot viewer to 800 by 500
    Then the "viewer regions" reading of scatter plot viewer should be 0
    When user picks "Tools > Draw Annotation Region" from the context menu of scatter plot viewer
    Then an info balloon containing "Click and drag to draw region" should have been shown
    And the "region drawing mode" reading of scatter plot viewer should be "true"
    When user drags the "top edge of view" area of scatter plot viewer to the "left edge of view" area
    Then the "region drawing mode" reading of scatter plot viewer should be "false"
    And the "viewer regions" reading of scatter plot viewer should be 1
    And scatter plot viewer should have a "region 1" area
    And "annotationRegions" property of scatter plot viewer should contain "\"type\":\"area\""
    And "Formula Lines" dialog should be visible
    And editor of "X column" input in "Formula Lines" dialog should have text "WEIGHT"
    And editor of "Y column" input in "Formula Lines" dialog should have text "HEIGHT"
    When user enters "Band A" into Title input in "Formula Lines" dialog
    And user clicks OK button in "Formula Lines" dialog
    Then the "Formula Lines" dialog should close
    And scatter plot viewer should have a "region Band A" area
    And scatter plot viewer should have a "region Band A title" area
    When user hovers over the "top edge of region Band A" area of scatter plot viewer
    Then the "regions hovered" reading of scatter plot viewer should be 1
    And tooltip should contain text "Band A"
    When user right-clicks on the "top edge of region Band A" area of scatter plot viewer
    Then the open menu should list "Edit..."
    And the open menu should list "Show Annotation Regions"
    When user picks "Edit..." from the open menu
    Then "Formula Lines" dialog should be visible
    And Title input in "Formula Lines" dialog should have value "Band A"
    When user enters "Band B" into Title input in "Formula Lines" dialog
    And user enters "Edited band" into Description input in "Formula Lines" dialog
    And user enters "3" into Width input in "Formula Lines" dialog
    And user clicks OK button in "Formula Lines" dialog
    Then scatter plot viewer should have a "region Band B" area
    And scatter plot viewer should not have a "region Band A" area
    And "annotationRegions" property of scatter plot viewer should contain "\"description\":\"Edited band\""
    And "annotationRegions" property of scatter plot viewer should contain "\"outlineWidth\":3"
    When user hovers over the "top edge of region Band B" area of scatter plot viewer
    Then tooltip should contain text "Band B"
    And tooltip should contain text "Edited band"
    When user moves the pointer away from scatter plot viewer
    And user picks "Tools > Formula Lines..." from the context menu of the "bottom edge of view" area of scatter plot viewer
    Then "Formula Lines" dialog should be visible
    When user clicks on Delete button in "Formula Lines" dialog
    And user clicks OK button in "Formula Lines" dialog
    Then the "viewer regions" reading of scatter plot viewer should be 0
    And scatter plot viewer should not have a "region Band B" area
    And no errors should have been logged

  Scenario: A lasso drawn on the density plot becomes a polygon region, and the dialog adds a formula region
    Given user adds a density plot viewer with:
      | xColumnName | WEIGHT |
      | yColumnName | HEIGHT |
      | lassoTool   | true   |
    And user resizes density plot viewer to 800 by 500
    When user picks "Tools > Draw Annotation Region" from the context menu of density plot viewer
    Then the "region drawing mode" reading of density plot viewer should be "true"
    When user drags a lasso over the "view" area of density plot viewer
    Then the "region drawing mode" reading of density plot viewer should be "false"
    And the "viewer regions" reading of density plot viewer should be 1
    And density plot viewer should have a "region 1" area
    And "annotationRegions" property of density plot viewer should contain "\"type\":\"area\""
    When user clicks OK button in "Formula Lines" dialog
    And user sets properties of density plot viewer:
      | lassoTool                    | false |
      | showViewerAnnotationRegions  | false |
    And user picks "Tools > Formula Lines..." from the context menu of density plot viewer
    Then "Formula Lines" dialog should be visible
    When user clicks on "ADD NEW" button in "Formula Lines" dialog
    And user picks "Region - Formula Lines" from the open menu
    And user clicks OK button in "Formula Lines" dialog
    Then the "viewer regions" reading of density plot viewer should be 2
    And "annotationRegions" property of density plot viewer should contain "\"type\":\"formula\""
    And "annotationRegions" property of density plot viewer should contain "\"formula1\":\"${HEIGHT} = ${WEIGHT}"
    And "annotationRegions" property of density plot viewer should contain "\"formula2\":\"${HEIGHT} = ${WEIGHT}"
    When user sets "showViewerAnnotationRegions" property of density plot viewer to "true"
    Then the "regions shown" reading of density plot viewer should be 2
    When user sets "annotationRegions" property of density plot viewer to "[]"
    Then the "viewer regions" reading of density plot viewer should be 0
    And no errors should have been logged

  Scenario: A vertical bar chart locks the drawn region to its full width
    Given user adds a bar chart viewer with:
      | splitColumnName | RACE     |
      | valueColumnName | AGE      |
      | valueAggrType   | avg      |
      | orientation     | vertical |
    And user resizes bar chart viewer to 800 by 500
    When user picks "Tools > Draw Annotation Region" from the context menu of the "top edge of view" area of bar chart viewer
    Then the "region drawing mode" reading of bar chart viewer should be "true"
    When user drags across the "top edge of view" area of bar chart viewer
    Then the "region drawing mode" reading of bar chart viewer should be "false"
    And the "viewer regions" reading of bar chart viewer should be 1
    And "annotationRegions" property of bar chart viewer should contain "\"type\":\"formula\""
    And "annotationRegions" property of bar chart viewer should contain "\"formula1\":\"${AGE} = "
    And "annotationRegions" property of bar chart viewer should contain "\"formula2\":\"${AGE} = "
    And bar chart viewer should have a "region 1" area
    And the "region 1" and "view" areas of bar chart viewer should be the same width
    And Title input in "Formula Lines" dialog should have value ""
    When user clicks OK button in "Formula Lines" dialog
    Then no errors should have been logged

  Scenario: A horizontal bar chart locks the drawn region to its full height
    Given user adds a bar chart viewer with:
      | splitColumnName | RACE       |
      | valueColumnName | AGE        |
      | valueAggrType   | avg        |
      | orientation     | horizontal |
    And user resizes bar chart viewer to 800 by 500
    When user picks "Tools > Draw Annotation Region" from the context menu of the "top edge of view" area of bar chart viewer
    Then the "region drawing mode" reading of bar chart viewer should be "true"
    When user drags across the "top edge of view" area of bar chart viewer
    Then the "region drawing mode" reading of bar chart viewer should be "false"
    And the "viewer regions" reading of bar chart viewer should be 1
    And "annotationRegions" property of bar chart viewer should contain "\"formula1\":\"${AGE} = "
    And bar chart viewer should have a "region 1" area
    And the "region 1" and "view" areas of bar chart viewer should be the same height
    When user clicks OK button in "Formula Lines" dialog
    Then no errors should have been logged

  Scenario: A box plot locks the drawn region to its full width
    Given user adds a box plot viewer with:
      | category1ColumnName | RACE |
      | valueColumnName     | AGE  |
    And user resizes box plot viewer to 800 by 500
    Then the "viewer regions" reading of box plot viewer should be 0
    When user picks "Tools > Draw Annotation Region" from the context menu of the "top edge of view" area of box plot viewer
    Then the "region drawing mode" reading of box plot viewer should be "true"
    When user drags across the "top edge of view" area of box plot viewer
    Then the "region drawing mode" reading of box plot viewer should be "false"
    And the "viewer regions" reading of box plot viewer should be 1
    And "annotationRegions" property of box plot viewer should contain "\"type\":\"formula\""
    And "annotationRegions" property of box plot viewer should contain "\"formula1\":\"${AGE} = "
    And "annotationRegions" property of box plot viewer should contain "\"formula2\":\"${AGE} = "
    And box plot viewer should have a "region 1" area
    And the "region 1" and "view" areas of box plot viewer should be the same width
    And Title input in "Formula Lines" dialog should have value ""
    When user clicks OK button in "Formula Lines" dialog
    Then no errors should have been logged

  Scenario: The dialog's colors and opacity reach the region
    Given user adds a scatter plot viewer with:
      | xColumnName       | AGE    |
      | yColumnName       | WEIGHT |
      | markerDefaultSize | 2      |
      | annotationRegions | [{"type":"area","x":"AGE","y":"WEIGHT","header":"Adults","area":[[30.5,30],[60.5,30],[60.5,180],[30.5,180]]}] |
    And user resizes scatter plot viewer to 800 by 500
    Then scatter plot viewer should have a "region Adults" area
    When user takes a snapshot of scatter plot viewer
    And user picks "Edit..." from the context menu of the "left edge of region Adults" area of scatter plot viewer
    Then "Formula Lines" dialog should be visible
    When user enters "#ff8800" into Region Color input in "Formula Lines" dialog
    And user enters "#003366" into Outline Color input in "Formula Lines" dialog
    And user enters "3" into Width input in "Formula Lines" dialog
    And user drags the slider of Opacity input in "Formula Lines" dialog to 60
    And user enters "#ff0000" into Color input in "Formula Lines" dialog
    And user clicks OK button in "Formula Lines" dialog
    Then the "Formula Lines" dialog should close
    And "annotationRegions" property of scatter plot viewer should contain "\"fillColor\""
    And "annotationRegions" property of scatter plot viewer should contain "\"outlineColor\""
    And "annotationRegions" property of scatter plot viewer should contain "\"outlineWidth\":3"
    And "annotationRegions" property of scatter plot viewer should contain "\"headerColor\""
    And "annotationRegions" property of scatter plot viewer should contain "\"opacity\""
    And scatter plot viewer should have repainted
    And scatter plot viewer should have a "region Adults" area
    And no errors should have been logged

  Scenario: A dataframe region added from the DataFrame tab is counted on its own
    Given user adds a scatter plot viewer with:
      | xColumnName | AGE    |
      | yColumnName | WEIGHT |
    And user resizes scatter plot viewer to 800 by 500
    Then the "viewer regions" reading of scatter plot viewer should be 0
    And the "dataframe regions" reading of scatter plot viewer should be 0
    When user picks "Tools > Formula Lines..." from the context menu of scatter plot viewer
    Then "Formula Lines" dialog should be visible
    When user clicks on DataFrame tab in "Formula Lines" dialog
    And user clicks on "ADD NEW" button in "Formula Lines" dialog
    And user picks "Region - Formula Lines" from the open menu
    And user enters "Medium weight" into Title input in "Formula Lines" dialog
    And user clicks OK button in "Formula Lines" dialog
    Then the "dataframe regions" reading of scatter plot viewer should be 1
    And the "viewer regions" reading of scatter plot viewer should be 0
    And scatter plot viewer should have a "region Medium weight" area
    When user sets "showDataframeAnnotationRegions" property of scatter plot viewer to "false"
    Then the "regions shown" reading of scatter plot viewer should be 0
    And the "dataframe regions" reading of scatter plot viewer should be 1
    When user sets the ".annotation-regions" tag of the table to ""
    Then the "dataframe regions" reading of scatter plot viewer should be 0
    And no errors should have been logged
