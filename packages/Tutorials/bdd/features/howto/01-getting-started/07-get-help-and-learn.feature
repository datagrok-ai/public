@guide @tutorial
Feature: How-to · Basics
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). Help where you are, then learning step by step, in one clip (it
  absorbed the former "Learn with interactive tutorials" clip): the Context Help toggle of the status bar (F1) docks a
  panel that shows the help page of whatever was clicked last (the grid's column, then a scatter plot) and closes with
  the same toggle; the sidebar's Help icon opens a menu with the documentation, the community forum, videos and
  feedback (Escape closes it); the menu behind the sidebar's logo has Help > Tutorials, which opens the Tutorials app
  with its tracks and cards; a card starts its tutorial (with a demog table of its own), whose bar on top says what to
  do next and ticks a step off as soon as it is done (Grid Customization: go a few rows down), and the tutorial closes
  from its bar. demog and its scatter plot are set up off camera; the help pages come from the stand, which the gate
  checks first. A running tutorial records its progress in the account's settings: the "tutorials" and
  "achievement-badges" settings are put back at feature end, and the Tutorials app is closed then. Nothing is made on
  the server.

  Scenario: Get help and learn
    Context Help explains what you click, the Help menu leads to the docs and community, and tutorials teach step by step.

    Given user is logged in
    And simple mode is off
    And the stand serves the help pages
    And the "tutorials" user settings are put back at feature end
    And the "achievement-badges" user settings are put back at feature end
    And the "Grid Customization" tutorial is not completed yet
    And the Tutorials app is closed
    And user opens demog dataset
    And user adds scatter plot viewer
    And the toolbox pane is hidden
    And the context panel is closed
    # caption: Context Help: the help page of what you click
    When user clicks on "Context Help" status bar toggle
    Then help panel should be visible
    # caption: Click a column: help on tables and columns
    When user clicks on the "header AGE" area of grid
    Then help panel should contain text "Dataframe"
    # caption: Click a viewer: help on that viewer
    When user clicks on empty plot space of scatter plot viewer
    Then help panel should contain text "scatterplot"
    # caption: The same toggle closes it
    When user clicks on "Context Help" status bar toggle
    Then help panel should be hidden
    # caption: The sidebar's Help: docs, community, videos, feedback
    When user clicks on "Help" tab
    And user hovers over "Help" menu item inside context menu
    # caption: Each one opens in a new browser tab
    And user hovers over "Community" menu item inside context menu
    And user hovers over "YouTube" menu item inside context menu
    And user presses Escape
    When the video says "To learn by doing, take a tutorial"
    # caption: The menu behind the logo: Help > Tutorials
    And user clicks on "Menu" tab
    And user picks "Help > Tutorials" from the open menu
    Then Tutorials view should be visible
    # caption: Tutorials come in tracks, from data analysis to chemistry
    When user hovers over "Grid Customization" tutorial card
    # caption: Click one to start it
    And user clicks on "Grid Customization" tutorial card
    Then the tutorial progress should be 1 of 30
    # caption: The bar on top says what to do next
    When user hovers over tutorial bar
    # caption: This icon lists all the steps
    And user clicks on all steps icon
    Then Tutorials panel should be visible
    # caption: Do it, and the step is ticked off
    When user clicks on the "cell 1 of AGE" area of grid
    And user presses ArrowDown
    And user presses ArrowDown
    And user presses ArrowDown
    Then the tutorial step "Go a few rows down in the grid" should be done
    # caption: Close it any time from its bar
    When user closes the tutorial
    Then no errors should have been logged
