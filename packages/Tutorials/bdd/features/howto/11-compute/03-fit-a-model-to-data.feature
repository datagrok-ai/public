@guide @tutorial
Feature: How-to · Compute
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). Fitting a model to measurements: the PK model of Diff Studio's
  library (two compartments, a dose every 12 hours), opened from Browse > Apps > Diff Studio, and Fit, which opens
  the fitting view. With Browse open, the fitting view (a table view) opened with the toolbox in its place on the
  left, which cut the result grid's charts in half (the status bar's Toolbox toggle only gave the dock back to
  Browse), so the clip closes Browse with the x on its title bar before Fit (the sidebar's compass,
  clicked with the pointer over the sidebar, does not close it): then no side panel comes. The target is the model's table output aimed at drug concentrations measured over three doses,
  picked with the table input's "Add file from Files" icon (Files > App Data > DiffStudio > library > pk.csv, read
  into the form without a view of its own); Count is set to the three doses of the data; the five PK parameters are
  switched on to be fitted within the ranges the model gives them. Run fits in the browser (Nelder-Mead in a web
  worker, no script): the grid lists the fits best first, with the fitted values and the model drawn over the data.
  A fit of all five recovers clearance 15.8 against 16.1 in an offline check (the data was not made with the
  model's defaults; fewer parameters fit poorly). Opening the model adds it to Diff Studio's recent models, a file in
  the account's home that is put back at feature end; nothing else is made on the server.

  Scenario: Fit a model to data
    Fitting finds the parameter values that make a model match what you measured.

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
    Then Fit ribbon item should be visible
    When user hovers over browse panel
    # caption: Close Browse for more room
    And user clicks on browse panel close icon
    Then browse panel should be hidden
    # caption: Fit opens a view of its own
    When user clicks on Fit ribbon item and steps off the ribbon
    Then the current view should show the text "Target"
    When the video says "Target: what the model should match"
    # caption: Here, a table of drug levels over time, from Files
    And user clicks on "Add file from Files" icon
    Then "Select a file" dialog should be visible
    When user expands "Files > App Data" tree node inside "Select a file" dialog
    And user expands "Files > App Data > DiffStudio" tree node inside "Select a file" dialog
    And user expands "Files > App Data > DiffStudio > library" tree node inside "Select a file" dialog
    And user scrolls "Files > App Data > DiffStudio > library > pk.csv" tree node inside "Select a file" dialog into view
    And user double-clicks on "Files > App Data > DiffStudio > library > pk.csv" tree node inside "Select a file" dialog
    Then "PK" input should have the value "pk"
    And "argument" input should have the value "Time, h"
    # caption: The data covers three doses: so must the model
    When user enters "3" into "Count" input
    # caption: Switch on the parameters to fit; each gets a range
    And user switches on the switch beside "Absorption rate" input
    And user switches on the switch beside "Clearance" input
    And user switches on the switch beside "Central volume" input
    And user switches on the switch beside "Peripheral volume" input
    And user switches on the switch beside "Intercompartmental clearance" input
    Then "Clearance (min)" input should be visible
    # caption: Run the fit
    When user clicks on "Run" icon
    Then the table should have a column "RMSE by iterations"
    When the video says "The best fit comes first: its values, and the model over the data"
    Then no errors should have been logged
