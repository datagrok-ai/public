@guide @tutorial
Feature: How-to · Basics
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). Arranging the workspace in one clip (it absorbed the former "Work
  with several tables side by side" clip): every open table has a tab and a click on a tab switches to it, a tab
  dragged to a side of another view docks the two side by side, and the splitter between them shares the width; then
  the toggles at the right end of the status bar show and hide the panels around the tables: Toolbox, Context Panel,
  Console (opened, then closed again), Tabs (the view tabs hidden and shown) and Presentation mode, which hides
  everything but the current view (cars; the docked iris goes too) and is left through its "back to design mode" link. cars and iris are opened off camera
  with the toolbox and the context panel hidden; the arrangement lives in the page only, nothing is made on the server.

  Once iris has left the document area, that area shows Home rather than cars (2026-11-22; the Link Tables how-to
  meets it too), so the clip clicks the cars tab to bring cars back on the left. Floating a view is left out: the
  view floats where the drag ended, half off the window.

  Scenario: Arrange your workspace
    Put two tables side by side, then show and hide the panels around them with the toggles at the bottom right.

    Given user is logged in
    And simple mode is off
    And user opens cars dataset
    And user opens iris dataset
    And the toolbox pane is hidden
    And the context panel is closed
    # caption: Each open table has its own tab: click one to switch
    When user clicks on the tab of the "cars" view
    Then the "cars" view should be current
    # caption: Drag a tab to the side of a view to see both
    When user drags the tab of the "iris" view to the right side of the "cars" view
    Then the "iris" view should be docked right-of the "cars" view
    # caption: Click cars to show it on the left
    When user clicks on the tab of the "cars" view
    Then the "cars" view should be current
    # caption: Drag the line between them to share the width
    When user resizes the "cars" view to 60% of the width it shares with the "iris" view
    Then the "cars" view should take about 60% of the width it shares with the "iris" view
    # caption: The toggles at the bottom right show and hide panels
    When user clicks on "Toolbox" status bar toggle
    Then toolbox should be visible
    And "Toolbox" status bar toggle should be selected
    # caption: Context Panel, also F4
    When user clicks on "Context Panel" status bar toggle
    Then context panel should be visible
    # caption: Console: commands and the log
    When user clicks on "Console" status bar toggle
    Then console should be visible
    # caption: Click again to hide it
    When user clicks on "Console" status bar toggle
    Then console should be hidden
    # caption: Tabs: hide the view tabs for more room
    When user clicks on "Tabs" status bar toggle
    Then "cars" view tab should be hidden
    When user clicks on "Tabs" status bar toggle
    Then "cars" view tab should be visible
    # caption: Presentation mode: just the current view, also F7
    When user clicks on "Presentation mode" status bar toggle
    Then status bar should be hidden
    And toolbox should be hidden
    # caption: This link brings the panels back
    When user clicks on "back to design mode" link
    Then status bar should be visible
    And toolbox should be visible
    And no errors should have been logged
