@guide @tutorial
Feature: How-to · Compute
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). Sensitivity analysis of the PK model of Diff Studio's library
  (two compartments, one 500 mg dose), opened from Browse > Apps > Diff Studio: Sensitivity opens its own view;
  Method is set to Sobol with 1000 samples, and the five PK parameters are switched on to vary within the ranges the
  model gives them. With Browse open, the analysis view (a table view) opened with the toolbox in its place on the
  left, which squeezed the form's labels and the result viewers, so the clip closes Browse with the x on its title bar
  before Sensitivity: then no side panel comes. The model's table output is read at its last row, the drug levels 12 hours after the dose. Run
  solves the model for every sample in the browser (no script) and draws Sobol's bar charts, how much of the spread
  of each output each parameter explains. Opening the model adds it to Diff Studio's recent models, a file in the
  account's home that is put back at feature end; nothing else is made on the server.

  Scenario: Find which inputs matter with sensitivity analysis
    Sensitivity analysis varies a model's inputs together and shows which of them change the result the most.

    Given user is logged in
    And simple mode is off
    And the "DiffStudio" package is installed
    And the "Compute2" package is installed
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
    # caption: A two-compartment PK model from the library
    When user double-clicks on "PK" model card
    Then Sensitivity ribbon item should be visible
    When user hovers over browse panel
    # caption: Close Browse for more room
    And user clicks on browse panel close icon
    Then browse panel should be hidden
    # caption: Sensitivity opens a view of its own
    When user clicks on Sensitivity ribbon item and steps off the ribbon
    Then "Method" input should be visible
    # caption: Sobol: how much of the spread each input explains
    When user selects "Sobol" in "Method" input
    # caption: Samples: more give steadier answers
    And user enters "1000" into "Samples" input
    Then "Samples" input should have the value "1000"
    # caption: Switch on the inputs to vary, each within its range
    When user switches on "Absorption rate" input
    And user switches on "Clearance" input
    And user switches on "Central volume" input
    And user switches on "Peripheral volume" input
    And user switches on "Intercompartmental clearance" input
    Then "Clearance min" input should be visible
    # caption: The output: drug levels at the last time point, 12 h
    When user hovers over "PK" input
    Then "PK" input should have the value "Last row"
    # caption: Run: the model is solved for every sample
    When user clicks on "Run" icon
    Then the open tableview should have 2 bar chart viewers
    When the video says "The longer the bar, the more it matters: here, clearance"
    # caption: Correlations give the direction: more clearance, lower level
    And user hovers over the "cell Clearance x Central conc., mg/L" area of correlation plot viewer
    Then no errors should have been logged
