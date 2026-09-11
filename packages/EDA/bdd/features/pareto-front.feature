@eda @realizes:eda.viewer.pareto-front
Feature: Pareto front viewer
  ML | Pareto Front... adds the viewer to the table. Translated from
  files/TestTrack/EDA/pareto-front-viewer.md and the package's playwright/pareto-front-viewer.test.ts.

  The viewer reports no readings of its own (no getWidgetStatus), so the axes scenario reads the
  repaint of the scatter plot inside it and the Auto Axes Selection the viewer turns off itself when
  an axis is chosen by hand. The objectives, chosen in the context panel, are
  pareto-front-objectives.feature.

  Background:
    Given user is logged in
    And user opens cars dataset

  Scenario: The viewer picks the column of unique values as its label
    When user picks "ML > Pareto Front..." from the top menu
    Then pareto front viewer should be visible
    And "Label Columns" property of pareto front viewer should be "model"
    And "Minimize" property of pareto front viewer should be "highway.mpg, price"
    And no errors should have been logged

  Scenario: On demog the unique subject id is the label
    Given user opens demog dataset
    When user picks "ML > Pareto Front..." from the top menu
    Then pareto front viewer should be visible
    And "Label Columns" property of pareto front viewer should be "USUBJID"

  Scenario: The axes of the viewer follow its properties
    When user picks "ML > Pareto Front..." from the top menu
    And user sets "X Axis" property of pareto front viewer to "horsepower"
    Then pareto front viewer should have repainted
    And "Auto Axes Selection" property of pareto front viewer should not be "true"
    And no errors should have been logged
