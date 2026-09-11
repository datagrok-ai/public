@journey @diffstudio @realizes:diffstudio.model.pk-pd
Feature: A cyclic model and its dosing inputs
  PK-PD from the library: its two plots, the Count input driven by its clickers, and the tooltips
  the model declares for its dosing inputs. Translated from
  files/TestTrack/DiffStudio/cyclic-models.md and the spec beside it.

  The tooltip texts are the ones the model's own .ivp file carries in square brackets
  (files/library/pk-pd.ivp); the old spec read that file at run time and fell back to these when it
  could not, so they are written here as what the product is expected to show.

  Background:
    Given user is logged in
    And user opens the "PK-PD" model of the Diff Studio library

  Scenario: The model arrives with its inputs and its plots
    Then the "PK-PD" view should be current
    And count input should be visible
    And Multiaxis tab should be visible
    And Facet tab should be visible

  Scenario: The clickers move Count and the solution follows
    When user clicks on Multiaxis tab
    And user takes a snapshot of line chart viewer
    And user hovers over count input
    And user clicks on plus icon in count input
    Then count input should not have value "1"
    And line chart viewer should have repainted

  Scenario: The dosing inputs explain themselves on hover
    When user hovers over begin input
    Then tooltip should contain text "Begin of dosing interval"
    When user hovers over end input
    Then tooltip should contain text "End of dosing interval"
    When user hovers over step input
    Then tooltip should contain text "Time step of simulation"
    And no errors should have been logged
