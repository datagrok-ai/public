@journey @diffstudio @realizes:diffstudio.model.bioreactor
Feature: Fitting a model to data
  The Fit command of the ribbon opening the fitting view over the Bioreactor model, and a fit that
  runs: a parameter varied within the bounds the model carries, the model's output as the target, and
  the experiment table as the data it is fitted to. Translated from
  files/TestTrack/DiffStudio/fitting.md and the spec beside it.

  What the form takes before Run does anything is the switch beside each input — one parameter to
  vary, one output to aim at — and the Run icon stays grey until every varied parameter has both
  bounds. FKox, which the manual case names, carries none in the model, so varying it alone leaves
  the fit unrunnable; FFox carries 0.15 to 0.25 and is varied here instead.

  The target table is read from the stand into the workspace rather than dragged out of the Browse
  tree: a drag into a Vue-rendered table input has no phrase, and opening the file in a view of its
  own takes the form off screen and rebuilds it.

  Background:
    Given user is logged in
    And user opens the "Bioreactor" model of the Diff Studio library

  Scenario: Process mode cascades into the parameters the fit would use
    When user clicks on Multiaxis tab
    And user takes a snapshot of line chart viewer
    And user selects "Mode 1" in "Process mode" input
    Then "Process mode" input should have value "Mode 1"
    And line chart viewer should have repainted

  Scenario: Fit opens a view of its own
    When user clicks on Fit ribbon item
    Then the "Bioreactor - fitting" view should be current
    And no errors should have been logged

  Scenario: A parameter is varied, and its bounds appear with it
    When user selects "Default" in "Process mode" input
    And user switches on "FFox" input
    Then "FFox (min)" input should be switched on
    And "FFox (min)" input should have value "0.15"
    And "FFox (max)" input should have value "0.25"
    When user enters "1.0" into "FFox (max)" input
    Then "FFox (max)" input should have value "1.0"

  Scenario: The fit needs an output to aim at and a table to aim with
    Given the "System:AppData/DiffStudio/library/bioreactor-experiment.csv" file is loaded as a table
    When user switches on "Bioreactor" input
    And user selects "bioreactor-experiment" in "Bioreactor" input
    Then "Bioreactor" input should have value "bioreactor-experiment"
    And "argument" input should have value "t"

  Scenario: Running the fit lowers the loss it reports, iteration by iteration
    When user clicks on "Run" icon
    Then the table should have a column "RMSE by iterations"
    And the table should have 1 row
    And every value of "FFox" column should lie between 0.15 and 1.0
    And the "RMSE by iterations" table should have at least 2 rows
    And the "Loss" column of the "RMSE by iterations" table should never increase
    And no errors should have been logged
