@guide @tutorial
Feature: How-to · More viewers
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). The network diagram from the Add viewer gallery on the Game
  of Thrones season 1 links (Files > Demo > got-s1-edges.csv: Source, Target, Weight, ...): Source and Target
  picked as its node columns on the viewer, then the settings icon: the Source nodes sized by the count of their Weight values (their links) and the edges widened by
  the sum of Weight (Data); last, a node clicked (its 33 rows selected; the click comes last because a selection
  moves the context panel off the viewer's settings). got-s1-edges is opened from Browse off camera; nothing is made on the server.

  Scenario: Show links between things with a network diagram
    A network diagram draws each value as a node and each row as a link between two nodes.

    Given user is logged in
    And simple mode is off
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    And user double-clicks on Files---Demo---got-s1-edges.csv tree node inside browse tree
    And the "got-s1-edges" view should be current
    And the browse panel is closed
    And the toolbox pane is hidden
    And the context panel is closed
    When the video says "A network diagram: who is linked to whom"
    And user clicks on "Add viewer" icon in toolbar
    And user types "network" into viewer gallery search in "Add Viewer" dialog
    And user clicks on first "Network diagram" card in "Add Viewer" dialog
    Then network diagram viewer should be visible
    # caption: Pick the two columns each link joins
    When user picks "Source" in the column selector at the "node 1 selector" area of network diagram viewer
    And user picks "Target" in the column selector at the "node 2 selector" area of network diagram viewer
    Then the "node 2 column" reading of network diagram viewer should be "Target"
    # caption: The rest is in the settings
    When user hovers over network diagram viewer
    And user clicks on settings icon of network diagram viewer
    Then the context panel should show "Network diagram"
    # caption: Data: bigger nodes for characters with more links
    When user selects "Weight" in "Node1 Size" property in context panel
    Then "node1SizeColumnName" property of network diagram viewer should be "Weight"
    # caption: and links as wide as the sum of their Weight
    When user selects "Weight" in "Edge Width" property in context panel
    And user selects "sum" in "Edge Width Aggr Type" property in context panel
    Then "edgeWidthAggrType" property of network diagram viewer should be "sum"
    # caption: Click a character to select the rows of their links
    When user clicks on the "node \"Eddard Stark\"" area of network diagram viewer
    Then 33 rows should be selected
    And no errors should have been logged
