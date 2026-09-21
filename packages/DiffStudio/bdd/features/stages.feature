@journey @diffstudio @realizes:diffstudio.model.acid-production
Feature: A staged model and its inputs
  Acid production from the library: its two plots, a stage duration typed straight into its input,
  and the tooltips the model declares. Translated from files/TestTrack/DiffStudio/stages.md and the
  spec beside it.

  The tooltip texts come from the model's own .ivp file (files/library/ga-production.ivp), the same
  source the old spec read at run time.

  Background:
    Given user is logged in
    And user opens the "Acid production" model of the Diff Studio library

  Scenario: The model arrives with its inputs and its plots
    Then the "Acid production" view should be current
    And "1-st stage" input should be visible
    And Multiaxis tab should be visible
    And Facet tab should be visible

  Scenario: Changing a stage duration redraws the solution
    When user clicks on Multiaxis tab
    And user takes a snapshot of line chart viewer
    And user enters "50" into "1-st stage" input
    Then "1-st stage" input should have value "50"
    And line chart viewer should have repainted

  Scenario: The inputs explain themselves on hover
    When user hovers over "1-st stage" input
    Then tooltip should contain text "Duration of the 1-st stage"
    When user hovers over biomass input
    Then tooltip should contain text "Aspergillus niger biomass"
    When user hovers over glucose input
    Then tooltip should contain text "Glucose"
    And no errors should have been logged
