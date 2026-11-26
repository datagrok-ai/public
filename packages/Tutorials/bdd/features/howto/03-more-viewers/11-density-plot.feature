@guide @tutorial
Feature: How-to · More viewers
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). The density plot from the toolbox on earthquakes: Longitude
  and Latitude picked on the plot (the plate boundaries appear), the densest bin hovered (its range and count),
  then the settings icon: more bins (Misc) and a logarithmic color scale (Color). earthquakes is opened off
  camera; nothing is made on the server.

  Scenario: See where points crowd with a density plot
    A density plot counts the points in small bins and colors each bin by its count, so dense areas stand out.

    Given user is logged in
    And simple mode is off
    And user opens earthquakes dataset
    And the context panel is closed
    And the browse panel is closed
    And the toolbox pane is shown
    When the video says "When points pile up, count them in bins"
    And user clicks on density plot icon on toolbox
    Then density plot viewer should be added to the open tableview
    # caption: Hide the toolbox to give the plot room
    When user clicks on toolbox tab
    Then toolbox should be hidden
    # caption: Pick the columns right on the plot
    When user picks "Longitude" in the "x" column selector of density plot viewer
    And user picks "Latitude" in the "y" column selector of density plot viewer
    Then the "x column" reading of density plot viewer should be "Longitude"
    # caption: Hover a bin for its range and count
    When user hovers over the "densest bin" area of density plot viewer
    Then tooltip should be visible
    # caption: The rest is in the settings
    When user hovers over density plot viewer
    And user clicks on settings icon of density plot viewer
    Then the context panel should show "Density plot"
    # caption: Misc: more bins, finer detail
    When user expands "Misc" category in context panel
    And user enters "100" into "Bins" property in context panel
    Then the "bins" reading of density plot viewer should be 100
    # caption: Color: a log scale brings out the sparse areas too
    When user expands "Color" category in context panel
    And user selects "logarithmic" in "Color Transform Type" property in context panel
    Then the "color transform" reading of density plot viewer should be "logarithmic"
    And no errors should have been logged
