@guide @tutorial
Feature: How-to · Dashboards
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). The dashboard "Demographics overview" (demog with a scatter
  plot, a histogram and a bar chart) is saved off camera through the ribbon's Save dialog, the Share offer
  cancelled. On camera Save opens the dialog again, the Presentation mode switch's tooltip says what it does, the
  switch is turned on and the dashboard saved. Every view is closed off camera; opened again from Browse >
  Dashboards, the dashboard comes up in presentation mode (no sidebar, toolbox or status bar, the "Press F7"
  balloon and the "back to design mode" link), and the link brings the panels back. The view's own menu and ribbon
  stay in presentation mode on this build (see projects-presentation-mode).

  The project, with its table, views and picture, is deleted when the feature starts and ends; the harness turns
  presentation mode off at every feature's start.

  TODO: not filmed yet. On localhost (2026-11-22) the switch is on when OK is pressed, but the server keeps no
  ".presentation" meta param for the project, and the reopened dashboard comes up in design mode; the feature fails
  on its server claim until that is fixed.

  Scenario: Open a dashboard in presentation mode
    Presentation mode opens a dashboard with just its charts, for people who only need to look at them.

    Given user is logged in
    And simple mode is off
    And no project named "Demographics overview" is on the server
    And user opens demog dataset
    And user adds a scatter plot viewer with:
      | xColumnName     | AGE    |
      | yColumnName     | WEIGHT |
      | colorColumnName | SEX    |
    And user adds a histogram viewer with:
      | valueColumnName | HEIGHT |
    And user adds a bar chart viewer with:
      | splitColumnName | RACE |
    And user clicks on Save ribbon item
    And user enters "Demographics overview" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    And the "Save project" dialog should close
    And "Share Demographics overview" dialog should be visible
    And user clicks on CANCEL button in "Share Demographics overview" dialog
    And the "Share Demographics overview" dialog should close
    And user closes all views
    And the browse panel is open
    And user clicks on Dashboards tree node inside browse tree
    And user enters "Demographics overview" into gallery search
    And "Demographics overview" gallery card should be visible
    And user double-clicks on "Demographics overview" gallery card and the project opens
    And the "demog" view should be current
    And the notification balloons are closed
    And the toolbox pane is shown
    And the context panel is closed
    And the browse panel is closed
    When the video says "Some people only need to look at the charts"
    And user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    # caption: Presentation mode: open it with only the charts showing
    When user hovers over "Presentation mode" input in "Save project" dialog
    Then tooltip should contain text "Design mode"
    When user switches on "Presentation mode" input in "Save project" dialog
    Then "Presentation mode" input in "Save project" dialog should be switched on
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And the project named "Demographics overview" should be saved to open in presentation mode
    Given user closes all views
    And the browse panel is open
    When the video says "Open it again from Browse > Dashboards"
    And user clicks on Dashboards tree node inside browse tree
    And user enters "Demographics overview" into gallery search
    Then "Demographics overview" gallery card should be visible
    # caption: It opens without the side panels and the status bar
    When user double-clicks on "Demographics overview" gallery card and the project opens
    Then scatter plot viewer should be visible
    And "back to design mode" link should be visible
    And browse tab should be hidden
    And toolbox tab should be hidden
    And status bar should be hidden
    And an info balloon containing "Press F7 to go back to the design mode" should have been shown
    # caption: back to design mode, or F7, brings the panels back
    When user clicks on "back to design mode" link
    Then status bar should be visible
    And browse tab should be visible
    And no errors should have been logged
