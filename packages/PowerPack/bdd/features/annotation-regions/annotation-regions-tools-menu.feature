@viewers @realizes:viewers.scatter-plot @realizes:viewers.density-plot @realizes:viewers.line-chart @realizes:viewers.histogram @realizes:viewers.bar-chart
Feature: The Tools menu offers the region items
  Every viewer that supports annotation regions lists Show Annotation Regions, Draw Annotation
  Region and Formula Lines... under Tools in its context menu. The Lasso Tool switch is a
  top-level item on the scatter plot, a Tools item on the density plot and the line chart, and
  absent on the one-axis viewers, where a drawn region is locked to the value axis. On demog-1000;
  from `annotation-regions.md` scenario 1.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset

  Scenario: The scatter plot keeps the Lasso Tool at the top level
    Given user adds a scatter plot viewer with:
      | xColumnName | WEIGHT |
      | yColumnName | HEIGHT |
    When user opens the context menu of scatter plot viewer
    Then the open menu should list "Tools > Show Annotation Regions"
    And the open menu should list "Tools > Draw Annotation Region"
    And the open menu should list "Tools > Formula Lines..."
    And the open menu should list "Lasso Tool"
    When user closes the context menu
    Then no errors should have been logged

  Scenario: The density plot has the Lasso Tool under Tools
    Given user adds a density plot viewer with:
      | xColumnName | WEIGHT |
      | yColumnName | HEIGHT |
    When user opens the context menu of density plot viewer
    Then the open menu should list "Tools > Show Annotation Regions"
    And the open menu should list "Tools > Draw Annotation Region"
    And the open menu should list "Tools > Formula Lines..."
    And the open menu should list "Tools > Lasso Tool"
    When user closes the context menu
    Then no errors should have been logged

  Scenario: The line chart has the Lasso Tool under Tools
    Given user adds a line chart viewer with:
      | xColumnName  | AGE    |
      | yColumnNames | HEIGHT |
    When user opens the context menu of line chart viewer
    Then the open menu should list "Tools > Show Annotation Regions"
    And the open menu should list "Tools > Draw Annotation Region"
    And the open menu should list "Tools > Formula Lines..."
    And the open menu should list "Tools > Lasso Tool"
    When user closes the context menu
    Then no errors should have been logged

  Scenario: The histogram offers no lasso
    Given user adds a histogram viewer with:
      | valueColumnName | AGE |
    When user opens the context menu of histogram viewer
    Then the open menu should list "Tools > Show Annotation Regions"
    And the open menu should list "Tools > Draw Annotation Region"
    And the open menu should list "Tools > Formula Lines..."
    And the open menu should not list "Tools > Lasso Tool"
    And the open menu should not list "Lasso Tool"
    When user closes the context menu
    Then no errors should have been logged

  Scenario: The bar chart offers no lasso
    Given user adds a bar chart viewer with:
      | splitColumnName | RACE     |
      | valueColumnName | AGE      |
      | valueAggrType   | avg      |
      | orientation     | vertical |
    When user opens the context menu of bar chart viewer
    Then the open menu should list "Tools > Show Annotation Regions"
    And the open menu should list "Tools > Draw Annotation Region"
    And the open menu should list "Tools > Formula Lines..."
    And the open menu should not list "Tools > Lasso Tool"
    And the open menu should not list "Lasso Tool"
    When user closes the context menu
    Then no errors should have been logged

  Scenario: The horizontal bar chart offers no lasso either
    Given user adds a bar chart viewer with:
      | splitColumnName | RACE       |
      | valueColumnName | AGE        |
      | valueAggrType   | avg        |
      | orientation     | horizontal |
    When user opens the context menu of bar chart viewer
    Then the open menu should list "Tools > Show Annotation Regions"
    And the open menu should list "Tools > Draw Annotation Region"
    And the open menu should list "Tools > Formula Lines..."
    And the open menu should not list "Tools > Lasso Tool"
    And the open menu should not list "Lasso Tool"
    When user closes the context menu
    Then no errors should have been logged

  Scenario: The box plot offers no lasso
    Given user adds a box plot viewer with:
      | category1ColumnName | RACE |
      | valueColumnName     | AGE  |
    When user opens the context menu of box plot viewer
    Then the open menu should list "Tools > Show Annotation Regions"
    And the open menu should list "Tools > Draw Annotation Region"
    And the open menu should list "Tools > Formula Lines..."
    And the open menu should not list "Tools > Lasso Tool"
    And the open menu should not list "Lasso Tool"
    When user closes the context menu
    Then no errors should have been logged
