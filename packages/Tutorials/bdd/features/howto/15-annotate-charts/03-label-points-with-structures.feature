@guide @tutorial
Feature: How-to · Annotations
  For the maintainer: structure labels on a scatter plot of spgi-100 (100 molecules) on its Chemical Space X and Y.
  The gear opens the plot's settings; in Labels, the "..." of Label Columns opens the column picker, where Structure
  and Id are ticked: a molecule column is drawn as the structure, with the ID under it, wherever there is room. The plot's
  context menu holds the same settings: Labels > Show Labels For > Selected (a radio item, which leaves the menu
  open: closed by a step of its own) keeps labels for the selected rows only, and a Shift-drag round the small cluster at
  X 6.5 to 7.2, Y 5.7 to 6.2 selects its four molecules (rows 59, 60, 87 and 89, all aminopiperidines), which get
  their labels. Label Content Size 1.6 draws the structures larger (picked in the panel, before Show Labels For: a pick from the
  panel's native list left that row's description tooltip over the plot until the next press). The Configure a scatter plot clip of the Viewers playlist labels points with a text column; this
  one is about molecules and about labelling only the points that matter. spgi-100 is opened and the plot added off
  camera; nothing is made on the server.

  Scenario: Label points with structures
    Show the structure and the ID next to the points of a chart, and only for the points you select.

    Given user is logged in
    And simple mode is off
    And user opens spgi dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    And user adds a scatter plot viewer with:
      | X | Chemical Space X |
      | Y | Chemical Space Y |
    And scatter plot viewer takes 65% of the space it shares with the grid
    When the video says "Labels put a column's values next to the points"
    # caption: The gear icon opens the plot's settings
    And user hovers over empty plot space of scatter plot viewer
    And user clicks on settings icon of scatter plot viewer
    Then the context panel should show "Scatter plot"
    # caption: Labels > Label Columns: pick what to show
    When user expands "Labels" category in context panel
    And user clicks on "..." button in "Label Columns" property in context panel
    Then "Select columns..." dialog should be visible
    When user toggles the "Structure" column in the column list of "Select columns..." dialog
    And user toggles the "Id" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    # caption: A molecule column is drawn as its structure
    Then "Label Columns" property of scatter plot viewer should contain "Structure"
    And "Label Columns" property of scatter plot viewer should contain "Id"
    And the "labels shown" reading of scatter plot viewer should be at least 5
    # caption: Label Content Size: bigger structures
    When user enters "1.6" into "Label Content Size" property in context panel
    Then "Label Content Size" property of scatter plot viewer should be "1.6"
    # caption: The plot's menu has them too. Show Labels For: Selected
    When user picks "Labels > Show Labels For > Selected" from the context menu of the "view" area of scatter plot viewer
    And user closes the context menu
    Then "Show Labels For" property of scatter plot viewer should be "Selected"
    And the "labels shown" reading of scatter plot viewer should be 0
    # caption: Shift+drag round some points: their labels appear
    When user drags a selection box on scatter plot viewer from 6.3, 6.5 to 7.5, 5.4 on its axes
    Then 4 rows should be selected
    And the "labels shown" reading of scatter plot viewer should be 4
    And no errors should have been logged
