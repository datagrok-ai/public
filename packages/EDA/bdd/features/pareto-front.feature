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
    And "Label Columns" property of pareto front viewer should be "model"
    And "Minimize" property of pareto front viewer should be "highway.mpg, price"
    And no errors should have been logged

  Scenario: On demog the unique subject id is the label
    Given user opens demog dataset
    When user picks "ML > Pareto Front..." from the top menu
    Then pareto front viewer should be visible
    And "Label Columns" property of pareto front viewer should be "USUBJID"

  Scenario: Without a column of unique values the label stays empty
    Given user opens iris dataset
    When user picks "ML > Pareto Front..." from the top menu
    Then pareto front viewer should be visible
    And "Label Columns" property of pareto front viewer should be ""
    And no errors should have been logged
