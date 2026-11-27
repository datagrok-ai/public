@guide @tutorial
Feature: How-to · More viewers
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). The Sankey diagram (Charts package) from the Add viewer
  gallery on the UK energy flows demo table (Files > Demo > energy_uk.csv: value, source, target), which it reads
  on its own; a node clicked (the 18 rows of its flows selected), then the settings icon: the links
  colored by their source (Color). energy_uk is opened from Browse off camera; nothing is made on the server.

  Scenario: Show flows with a Sankey diagram
    A Sankey diagram shows how an amount flows from one category to the next; the wider the band, the bigger the flow.

    Given user is logged in
    And simple mode is off
    And the package autostarts have completed
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    And user double-clicks on Files---Demo---energy_uk.csv tree node inside browse tree
    And the "energy_uk" view should be current
    And the browse panel is closed
    And the toolbox pane is hidden
    And the context panel is closed
    When the video says "A Sankey diagram shows amounts flowing between categories"
    And user clicks on "Add viewer" icon in toolbar
    And user types "sankey" into viewer gallery search in "Add Viewer" dialog
    And user clicks on first "Sankey" card in "Add Viewer" dialog
    Then sankey viewer should be visible
    And "Source" property of sankey viewer should be "source"
    When the video says "It found the source, target and amount columns itself" at the top
    # caption: Click a node to select the rows of its flows
    And user clicks on the "node Electricity grid" area of sankey viewer
    Then 18 rows should be selected
    # caption: The rest is in the settings
    When user hovers over sankey viewer
    And user clicks on settings icon of sankey viewer
    Then the context panel should show "Sankey"
    # caption: Color: color the links by where they come from
    When user expands "Color" category in context panel
    And user selects "source" in "Link Color" property in context panel
    Then "Link Color Column Name" property of sankey viewer should be "source"
    And no errors should have been logged
