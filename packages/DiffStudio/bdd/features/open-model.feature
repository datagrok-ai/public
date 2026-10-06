@journey @diffstudio @realizes:diffstudio.app.diff-studio @realizes:diffstudio.model.bioreactor
Feature: Opening a model from the library
  The Bioreactor model as the library serves it: its charts, its inputs, and the two things the
  manual case asks to watch — that a change to an input reaches the table and the chart at once,
  and that Process mode cascades into the parameters below it. Translated from
  files/TestTrack/DiffStudio/open-model.md and the spec beside it.

  The model's charts are real viewers of a table view named after the model, so what the old spec
  did by hashing a canvas is asked of the platform instead: a repaint is the viewer's own word.

  The Facet plot draws twelve small multiples that are not viewers of the table view, so the tier's
  per-viewer colour count does not reach them; the colours are counted over every canvas of the view
  instead, the reduction the old spec used (alpha, paper, ink and greys dropped, four bits a
  channel) as a step of its own.

  Background:
    Given user is logged in
    And user opens the "Bioreactor" model of the Diff Studio library

  Scenario: The model arrives with its table and its inputs
    Then the "Bioreactor" view should be current
    And grid should be visible
    And "Process mode" input should be visible
    And "switch at" input should have value "135"
    And FFox input should have value "0.20"

  Scenario: The line chart offers the Multiaxis and Facet tabs
    Then Multiaxis tab should be visible
    And Facet tab should be visible
    And Grid tab should be present

  Scenario: The Facet tab draws the model as small multiples
    When user clicks on Facet tab
    Then Facet tab should be selected
    And Multiaxis tab should not be selected
    And the canvases of open tableview should be painted in at least 10 colors

  Scenario: Changing "switch at" redraws the table and the chart
    When user clicks on Multiaxis tab
    And user takes a snapshot of line chart viewer
    And user enters "150" into "switch at" input
    Then "switch at" input should have value "150"
    And the page address should contain "switchat=150"
    And line chart viewer should have repainted

  Scenario: The slider moves "switch at" and the chart follows
    When user takes a snapshot of line chart viewer
    And user drags the slider of "switch at" input to 100
    Then "switch at" input should have a value between 95 and 105
    And line chart viewer should have repainted

  Scenario: Process mode cascades into the parameters below it
    When user takes a snapshot of line chart viewer
    And user selects "Mode 1" in "Process mode" input
    Then "Process mode" input should have value "Mode 1"
    And FFox input should not have value "0.20"
    And line chart viewer should have repainted
    And no errors should have been logged
