@guide @tutorial
Feature: How-to · Compute
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). Diff Studio from Browse > Apps: its hub (templates, the library),
  the Lotka-Volterra predator-prey model opened from a library card, the prey's starting value typed in, a
  parameter's tooltip and its slider dragged, and the end time raised, each recomputing the solution in the browser
  (diff-grok, no script). The browse tab opens Browse, which stays open (Diff Studio keeps it shown while a model is current). Opening a model adds it to
  Diff Studio's recent models, a file in the account's home that is put back at feature end; nothing else is made
  on the server.

  Scenario: Run a model in Diff Studio
    Diff Studio solves differential equations right in the browser. Open a model, change its inputs and watch the solution follow.

    Given user is logged in
    And simple mode is off
    And the "DiffStudio" package is installed
    And the package autostarts have completed
    And the file "diff-studio-recent.d42" of the user's home is put back at feature end
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    When the video says "Diff Studio is one of the apps in Browse"
    And user clicks on browse tab
    And user clicks on Apps tree node inside browse tree
    And user double-clicks on Diff-Studio gallery card
    Then the "Diff Studio" view should be current
    # caption: Templates to start from, and a library of ready models
    When user hovers over "Lotka-Volterra" model card
    And user double-clicks on "Lotka-Volterra" model card
    Then the "Lotka-Volterra" view should be current
    And "Prey" input should be visible
    When the video says "Inputs on the left, the solution on the right"
    # caption: Change a starting value: the curves are recomputed at once
    And user enters "30" into "Prey" input
    Then "Prey" input should have the value "30"
    # caption: Point at an input to see what it stands for
    When user hovers over "alpha" input
    # caption: Drag its slider and watch the curves follow
    And user drags the slider of "alpha" input to 0.6
    Then "alpha" input should have a value between 0.5 and 0.7
    # caption: Finish: how far in time to solve
    When user enters "150" into "Finish" input
    Then "Finish" input should have the value "150"
    And no errors should have been logged
