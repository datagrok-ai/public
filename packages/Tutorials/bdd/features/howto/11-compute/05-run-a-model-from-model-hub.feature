@guide @tutorial
Feature: How-to · Compute
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). Model Hub (Browse > Apps > Model Hub, Compute2's catalog of the
  functions tagged as models): a card pointed at, then Ball Flight Simulation (DiffStudio's ballFlight, JavaScript,
  solved in the browser) run with a double-click; its form takes the throw's velocity and angle, and every change
  runs the model again, redrawing the trajectory and the maximum distance and height; the ribbon's Fit inputs icon
  is pointed at. Runs go to the model's history only through its save icon, so nothing is made on the server.

  Scenario: Run a model from Model Hub
    Model Hub lists the models you can run on your server. Open one and it runs with a form for its inputs.

    Given user is logged in
    And simple mode is off
    And the "DiffStudio" package is installed
    And the "Compute2" package is installed
    And the package autostarts have completed
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    # caption: Model Hub is in Browse > Apps
    When user clicks on browse tab
    And user clicks on Apps tree node inside browse tree
    And user double-clicks on Model-Hub gallery card
    Then "Ball Flight Simulation" link should be visible
    # caption: Each card is a model; point at one to read about it
    When user hovers over "Ball Flight Simulation" link
    # caption: Double-click to run it
    And user double-clicks on "Ball Flight Simulation" link
    Then "Angle" input should be visible
    When the video says "Inputs on the left; the results update as you change them"
    And user enters "30" into "Angle" input
    Then "Angle" input should have the value "30"
    When user drags the slider of "Velocity" input to 58
    Then "Velocity" input should have a value between 57 and 59
    When the video says "Max distance and height, at the bottom, follow each change" at the top
    And the video says "Fitting and sensitivity analysis are on its ribbon too" at the top
    And user hovers over "Fit inputs" icon
    Then no errors should have been logged
