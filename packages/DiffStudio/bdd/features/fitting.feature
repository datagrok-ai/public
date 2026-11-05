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
  own takes the form off screen and rebuilds it. That the fit's loss falls and its result stays within
  the bounds is the fitting library's claim (LibTests compute-utils/fitting), not this feature's; Process
  mode cascading into the inputs is open-model's.

  Background:
    Given user is logged in
    And user opens the "Bioreactor" model of the Diff Studio library

  Scenario: Fit opens a view of its own
    When user clicks on Fit ribbon item
    Then the "Bioreactor - fitting" view should be current
    And no errors should have been logged

  Scenario: A parameter is varied, and its bounds appear with it
    When user switches on "FFox" input
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

  Scenario: Running the fit gives its result and the loss of each iteration
    When user clicks on "Run" icon
    Then the table should have a column "RMSE by iterations"
    And the table should have 1 row
    And the "RMSE by iterations" table should have at least 2 rows
    And no errors should have been logged
