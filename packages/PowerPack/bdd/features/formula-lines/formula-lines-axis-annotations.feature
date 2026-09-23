@viewers @realizes:viewers.scatter-plot @realizes:viewers.density-plot @realizes:viewers.line-chart @realizes:viewers.histogram @realizes:viewers.bar-chart @realizes:powerpack.dialogs.formula-lines
Feature: The Annotations group on a numeric axis
  The context menu of a numeric axis carries an Annotations group with Add Line, Add Band and Add
  Region. Each item binds a new formula item to the axis column at the column's quartiles — the
  line at the median, the band and the region between the first and third quartile — writes it
  into the viewer's look and opens PowerPack's Formula Lines dialog on it. A categorical axis and
  a derived one (the histogram's counts, the bar chart's categories) have no Annotations group.
  On demog-1000: WEIGHT's median is 77.4 and its quartiles 64.2 and 91.0; HEIGHT's 168.5,
  160.9 and 177.6; AGE's 45, 36 and 56; the bar chart's avg(AGE) axis runs 45.7 to 47.1 around
  46.2. The axis is right-clicked at its end, away from the column selector in its middle.
  From `formula-lines-axis-annotations.md`; every scenario empties what it added.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset

  Scenario: The scatter plot's X axis adds a line, a band and a region on WEIGHT
    Given user adds a scatter plot viewer with:
      | xColumnName | WEIGHT |
      | yColumnName | HEIGHT |
    And user resizes scatter plot viewer to 800 by 500
    When user picks "Annotations > Add Line" from the context menu of the "left edge of x axis" area of scatter plot viewer
    Then "Formula Lines" dialog should be visible
    And editor of Column input in "Formula Lines" dialog should have text "WEIGHT"
    And Value input in "Formula Lines" dialog should have value "77.4"
    When user clicks OK button in "Formula Lines" dialog
    Then "formulaLines" property of scatter plot viewer should contain "${WEIGHT} = 77.4"
    And the "formula lines" reading of scatter plot viewer should be 1
    When user sets "formulaLines" property of scatter plot viewer to ""
    Then the "formula lines" reading of scatter plot viewer should be 0
    When user picks "Annotations > Add Band" from the context menu of the "left edge of x axis" area of scatter plot viewer
    And user clicks OK button in "Formula Lines" dialog
    Then "formulaLines" property of scatter plot viewer should contain "${WEIGHT} in (64.2, 91.0)"
    And the "formula lines" reading of scatter plot viewer should be 1
    When user sets "formulaLines" property of scatter plot viewer to ""
    And user picks "Annotations > Add Region" from the context menu of the "left edge of x axis" area of scatter plot viewer
    Then "Formula Lines" dialog should be visible
    And Title input in "Formula Lines" dialog should have value ""
    When user clicks OK button in "Formula Lines" dialog
    Then "annotationRegions" property of scatter plot viewer should contain "${WEIGHT} = 64.2"
    And "annotationRegions" property of scatter plot viewer should contain "${WEIGHT} = 91.0"
    And the "viewer regions" reading of scatter plot viewer should be 1
    And scatter plot viewer should have a "region 1" area
    And the "region 1" and "view" areas of scatter plot viewer should be the same height
    When user sets "annotationRegions" property of scatter plot viewer to "[]"
    Then the "viewer regions" reading of scatter plot viewer should be 0
    And no errors should have been logged

  Scenario: The scatter plot's Y axis adds a line, a band and a region on HEIGHT
    Given user adds a scatter plot viewer with:
      | xColumnName | WEIGHT |
      | yColumnName | HEIGHT |
    And user resizes scatter plot viewer to 800 by 500
    When user picks "Annotations > Add Line" from the context menu of the "bottom edge of y axis" area of scatter plot viewer
    Then editor of Column input in "Formula Lines" dialog should have text "HEIGHT"
    And Value input in "Formula Lines" dialog should have value "168.5"
    When user clicks OK button in "Formula Lines" dialog
    Then "formulaLines" property of scatter plot viewer should contain "${HEIGHT} = 168.5"
    And the "formula lines" reading of scatter plot viewer should be 1
    When user sets "formulaLines" property of scatter plot viewer to ""
    And user picks "Annotations > Add Band" from the context menu of the "bottom edge of y axis" area of scatter plot viewer
    And user clicks OK button in "Formula Lines" dialog
    Then "formulaLines" property of scatter plot viewer should contain "${HEIGHT} in (160.9, 177.6)"
    When user sets "formulaLines" property of scatter plot viewer to ""
    And user picks "Annotations > Add Region" from the context menu of the "bottom edge of y axis" area of scatter plot viewer
    And user clicks OK button in "Formula Lines" dialog
    Then "annotationRegions" property of scatter plot viewer should contain "${HEIGHT} = 160.9"
    And "annotationRegions" property of scatter plot viewer should contain "${HEIGHT} = 177.6"
    And scatter plot viewer should have a "region 1" area
    And the "region 1" and "view" areas of scatter plot viewer should be the same width
    When user sets "annotationRegions" property of scatter plot viewer to "[]"
    Then no errors should have been logged

  Scenario: The density plot's axes add items on their columns
    Given user adds a density plot viewer with:
      | xColumnName | WEIGHT |
      | yColumnName | HEIGHT |
    And user resizes density plot viewer to 800 by 500
    When user picks "Annotations > Add Line" from the context menu of the "left edge of x axis" area of density plot viewer
    And user clicks OK button in "Formula Lines" dialog
    Then "formulaLines" property of density plot viewer should contain "${WEIGHT} = 77.4"
    When user picks "Annotations > Add Band" from the context menu of the "bottom edge of y axis" area of density plot viewer
    And user clicks OK button in "Formula Lines" dialog
    Then "formulaLines" property of density plot viewer should contain "${HEIGHT} in (160.9, 177.6)"
    When user picks "Annotations > Add Region" from the context menu of the "left edge of x axis" area of density plot viewer
    And user clicks OK button in "Formula Lines" dialog
    Then "annotationRegions" property of density plot viewer should contain "${WEIGHT} = 64.2"
    And density plot viewer should have a "region 1" area
    When user sets properties of density plot viewer:
      | formulaLines      |    |
      | annotationRegions | [] |
    Then the "viewer regions" reading of density plot viewer should be 0
    And no errors should have been logged

  Scenario: The line chart's axes add items on AGE and on the aggregated HEIGHT
    Given user adds a line chart viewer with:
      | xColumnName  | AGE    |
      | yColumnNames | HEIGHT |
    And user resizes line chart viewer to 800 by 500
    When user picks "Annotations > Add Line" from the context menu of the "left edge of x axis" area of line chart viewer
    And user clicks OK button in "Formula Lines" dialog
    Then "formulaLines" property of line chart viewer should contain "${AGE} = 45.0"
    And line chart viewer should have a "formula line AGE = 45.0" area
    When user sets "formulaLines" property of line chart viewer to ""
    And user picks "Annotations > Add Line" from the context menu of the "bottom edge of y axis" area of line chart viewer
    And user clicks OK button in "Formula Lines" dialog
    Then "formulaLines" property of line chart viewer should contain "${avg(HEIGHT)} = 168.8"
    And line chart viewer should have a "formula line avg(HEIGHT) = 168.8" area
    When user sets "formulaLines" property of line chart viewer to ""
    And user picks "Annotations > Add Band" from the context menu of the "left edge of x axis" area of line chart viewer
    And user clicks OK button in "Formula Lines" dialog
    Then "formulaLines" property of line chart viewer should contain "${AGE} in (36.0, 56.0)"
    And line chart viewer should have a "formula band AGE in (36.0, 56.0)" area
    When user sets "formulaLines" property of line chart viewer to ""
    And user picks "Annotations > Add Region" from the context menu of the "bottom edge of y axis" area of line chart viewer
    And user clicks OK button in "Formula Lines" dialog
    Then "annotationRegions" property of line chart viewer should contain "${avg(HEIGHT)} = 164.9"
    And "annotationRegions" property of line chart viewer should contain "${avg(HEIGHT)} = 171.8"
    And line chart viewer should have a "region 1" area
    When user sets "annotationRegions" property of line chart viewer to "[]"
    Then no errors should have been logged

  Scenario: The histogram's value axis has the group and its count axis does not
    Given user adds a histogram viewer with:
      | valueColumnName | AGE |
    And user resizes histogram viewer to 800 by 500
    When user picks "Annotations > Add Line" from the context menu of the "left edge of x axis" area of histogram viewer
    And user clicks OK button in "Formula Lines" dialog
    Then "formulaLines" property of histogram viewer should contain "${AGE} = 45.0"
    When user sets "formulaLines" property of histogram viewer to ""
    And user picks "Annotations > Add Region" from the context menu of the "left edge of x axis" area of histogram viewer
    And user clicks OK button in "Formula Lines" dialog
    Then "annotationRegions" property of histogram viewer should contain "${AGE} = 36.0"
    And "annotationRegions" property of histogram viewer should contain "${AGE} = 56.0"
    And histogram viewer should have a "region 1" area
    And the "region 1" and "view" areas of histogram viewer should be the same height
    When user sets "annotationRegions" property of histogram viewer to "[]"
    And user right-clicks on the "y axis" area of histogram viewer
    Then the open menu should not list "Annotations"
    And the open menu should not list "Add Line"
    When user closes the context menu
    Then no errors should have been logged

  Scenario: The bar chart's value axis has the group and its category axis does not
    Given user adds a bar chart viewer with:
      | splitColumnName | RACE     |
      | valueColumnName | AGE      |
      | valueAggrType   | avg      |
      | orientation     | vertical |
    And user resizes bar chart viewer to 800 by 500
    When user picks "Annotations > Add Line" from the context menu of the "bottom edge of x axis" area of bar chart viewer
    And user clicks OK button in "Formula Lines" dialog
    Then "formulaLines" property of bar chart viewer should contain "${AGE} = 46.2"
    When user sets "formulaLines" property of bar chart viewer to ""
    And user picks "Annotations > Add Band" from the context menu of the "bottom edge of x axis" area of bar chart viewer
    And user clicks OK button in "Formula Lines" dialog
    Then "formulaLines" property of bar chart viewer should contain "${AGE} in (45.7, 47.1)"
    When user sets "formulaLines" property of bar chart viewer to ""
    And user picks "Annotations > Add Region" from the context menu of the "bottom edge of x axis" area of bar chart viewer
    And user clicks OK button in "Formula Lines" dialog
    Then "annotationRegions" property of bar chart viewer should contain "${AGE} = 45.7"
    And "annotationRegions" property of bar chart viewer should contain "${AGE} = 47.1"
    And bar chart viewer should have a "region 1" area
    When user sets "annotationRegions" property of bar chart viewer to "[]"
    And user right-clicks on the "y axis" area of bar chart viewer
    Then the open menu should not list "Annotations"
    When user closes the context menu
    Then no errors should have been logged
