@guide @tutorial
Feature: How-to · Compute
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). Diff Studio's Create button opens the Basic template in the
  equations editor; the clip turns it into logistic population growth, dN/dt = r * N * (1 - N / K): the name, the
  equation, the starting value with a slider range and a hint, and a #parameters block with r and K. Refresh solves
  it in the browser (diff-grok, no script) and names the view after the model; Edit turned off shows the inputs the
  model declares, and r's slider changes how fast the curve rises. Diff Studio is opened from Browse > Apps, and the
  editor is widened by dragging the dock divider beside it. Running a model adds it to Diff Studio's recent models, a file in the account's home that is put back at
  feature end; nothing else is made on the server (the model is never saved).

  Scenario: Write your own model
    Type the equations of your own model into Diff Studio, and it builds the inputs and the charts for you.

    Given user is logged in
    And simple mode is off
    And the "DiffStudio" package is installed
    And the package autostarts have completed
    And the file "diff-studio-recent.d42" of the user's home is put back at feature end
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    # caption: Diff Studio is in Browse > Apps
    When user clicks on browse tab
    And user clicks on Apps tree node inside browse tree
    And user double-clicks on Diff-Studio gallery card
    Then the "Diff Studio" view should be current
    # caption: Create starts a new model from a template
    When user clicks on Create button
    Then code editor should be visible
    # caption: Drag the divider to widen the editor
    When user drags equations editor divider by 300 and 0 pixels
    When the video says "A model is plain text: equations, time, starting values, parameters"
    # caption: Give it a name
    And user replaces the line starting with "#name" in code editor with "#name: Population growth"
    # caption: One line per equation: logistic growth
    And user replaces the line starting with "dy/dt" in code editor with "dN/dt = r * N * (1 - N / K)"
    # caption: The starting value, with a range for its slider and a hint
    And user replaces the line starting with "y = 0" in code editor with "N = 10 {min: 1; max: 100} [Starting population]"
    # caption: Parameters get inputs of their own
    And user puts "#parameters:" on a new line after the line starting with "N = 10" in code editor
    And user puts "  r = 0.5 {min: 0.1; max: 2} [Growth rate]" on a new line after the line starting with "#parameters" in code editor
    And user puts "  K = 1000 {min: 100; max: 2000} [Carrying capacity]" on a new line after the line starting with "r = 0.5" in code editor
    # caption: Refresh solves it
    And user clicks on Refresh ribbon item
    Then the "Population growth" view should be current
    # caption: Turn Edit off to get back to the inputs
    When user clicks on Edit ribbon item
    Then code editor should be absent
    And "r" input should be visible
    # caption: Every parameter is now an input with a slider
    When user drags the slider of "r" input to 1.5
    Then "r" input should have a value between 1.4 and 1.6
    And no errors should have been logged
