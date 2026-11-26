@guide @tutorial
Feature: How-to · More viewers
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). The Radar (Charts package) from the Add viewer gallery on
  cars: the settings icon, five columns picked as its axes in the Values dialog (Value), the values of the current
  row shown on the axes (Misc), then a click on porsche and on volvo1 in the grid draws each car's shape over
  the band where most cars fall (the current row's green line shows grey while the pointer is still on that
  grid row, its mouse-over line drawn over it). cars is opened off camera; nothing is made on the server.

  Scenario: Compare rows across many columns with a radar
    A radar draws a row as a shape over several axes, so you can see where it stands against the rest.

    Given user is logged in
    And simple mode is off
    And the package autostarts have completed
    And user opens cars dataset
    And the context panel is closed
    And the browse panel is closed
    And the toolbox pane is hidden
    When the video says "A radar shows a row as a shape across several columns"
    And user clicks on "Add viewer" icon in toolbar
    And user types "radar" into viewer gallery search in "Add Viewer" dialog
    And user clicks on first "Radar" card in "Add Viewer" dialog
    Then radar viewer should be visible
    # caption: Pick the columns for the axes in the settings
    When user hovers over radar viewer
    And user clicks on settings icon of radar viewer
    And user expands "Value" category in context panel
    And user clicks on "..." button in "Values" property in context panel
    Then "Select columns..." dialog should be visible
    When user clicks on "None" link in "Select columns..." dialog
    And user toggles the "curb.weight" column in the column list of "Select columns..." dialog
    And user toggles the "eng.size" column in the column list of "Select columns..." dialog
    And user toggles the "horsepower" column in the column list of "Select columns..." dialog
    And user toggles the "city.mpg" column in the column list of "Select columns..." dialog
    And user toggles the "price" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then "Values" property of radar viewer should be "curb.weight, eng.size, horsepower, city.mpg, price"
    # caption: Misc: write the values of the current row on the axes
    When user expands "Misc" category in context panel
    And user checks "Show Values" property in context panel
    When the video says "The shaded band is where most cars fall" at the top
    # caption: Click a car to draw its shape over the rest
    And user clicks on the "cell 19 of model" area of grid
    Then row 19 should be current
    # caption: Another car, another shape
    When user clicks on the "cell 28 of model" area of grid
    Then row 28 should be current
    And no errors should have been logged
