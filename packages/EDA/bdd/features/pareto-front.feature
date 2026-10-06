@eda @realizes:eda.viewer.pareto-front
Feature: Pareto front viewer
  ML | Pareto Front... adds the viewer to the table. Translated from
  files/TestTrack/EDA/pareto-front-viewer.md and the package's playwright/pareto-front-viewer.test.ts.

  iris has no category of unique values (Species repeats), so its label stays empty. The properties
  chosen in the context panel, the objectives, the axes and the labels, are
  pareto-front-objectives.feature.

  Background:
    Given user is logged in
    And user opens cars dataset

  Scenario: The viewer picks the column of unique values as its label
    When user picks "ML > Pareto Front..." from the top menu
    Then pareto front viewer should be visible
    And scatter plot viewer in pareto front viewer should be painted
    And "Label Columns" property of pareto front viewer should be "model"
    And "Label Columns" property of scatter plot viewer in pareto front viewer should be "model"
    And the "labels shown" reading of scatter plot viewer in pareto front viewer should be at least 1
    And "Minimize" property of pareto front viewer should be "highway.mpg, price"
    And no errors should have been logged

  Scenario: On demog the unique subject id is the label
    Given user opens demog dataset
    When user picks "ML > Pareto Front..." from the top menu
    Then pareto front viewer should be visible
    And scatter plot viewer in pareto front viewer should be painted
    And "Label Columns" property of pareto front viewer should be "USUBJID"
    And "Label Columns" property of scatter plot viewer in pareto front viewer should be "USUBJID"
    And the "labels shown" reading of scatter plot viewer in pareto front viewer should be at least 1

  Scenario: Without a column of unique values the label stays empty
    Given user opens iris dataset
    When user picks "ML > Pareto Front..." from the top menu
    Then pareto front viewer should be visible
    And scatter plot viewer in pareto front viewer should be painted
    And "Minimize" property of pareto front viewer should be "Petal.Length, Petal.Width"
    And "Label Columns" property of pareto front viewer should be ""
    And "Label Columns" property of scatter plot viewer in pareto front viewer should be ""
    And the "labels shown" reading of scatter plot viewer in pareto front viewer should be 0
    And no errors should have been logged
