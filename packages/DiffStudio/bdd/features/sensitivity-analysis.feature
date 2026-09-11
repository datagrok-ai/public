@journey @diffstudio @realizes:diffstudio.model.bioreactor
Feature: Sensitivity analysis over a model
  The Sensitivity command of the ribbon opening its own view over the Bioreactor model. Translated
  from files/TestTrack/DiffStudio/sensitivity-analysis.md and the spec beside it.

  Choosing which parameters to vary is the switch each input carries: the analysis has nothing to
  vary until one is on, and a run started with none of them on ends with a single viewer. The switch
  is not part of the input's own control, so the library resolves it as the switch that governs the
  input — inside its host here, beside it in the fitting form. A parameter that is being varied is
  shown as a min and a max rather than as itself, and the switch moves to the min, which is why the
  round trip below names three different inputs for one parameter.

  Background:
    Given user is logged in
    And user opens the "Bioreactor" model of the Diff Studio library

  Scenario: The model is the one the ribbon acts on
    Then the "Bioreactor" view should be current
    And Sensitivity ribbon item should be visible
    And Fit ribbon item should be visible

  Scenario: The Edit toggle swaps the form for the equations, and back
    When user clicks on Edit ribbon item
    Then code editor should be visible
    And "Process mode" input should be absent
    When user clicks on Edit ribbon item
    Then code editor should be absent
    And "Process mode" input should be visible

  Scenario: Sensitivity opens a view of its own
    When user clicks on Sensitivity ribbon item
    Then the "Bioreactor - comparison" view should be current

  Scenario: A parameter is chosen by the switch beside its input
    When user switches on "FFox" input
    Then "FFox min" input should be switched on
    And "FFox max" input should be visible
    When user switches off "FFox min" input
    Then "FFox" input should be visible

  Scenario: Running the analysis over three parameters puts its viewers on screen
    When user switches on "FFox" input
    And user switches on "FKox" input
    And user switches on "FFred" input
    And user clicks on "Run" icon
    Then the current view should hold at least 4 viewers
    And no errors should have been logged
