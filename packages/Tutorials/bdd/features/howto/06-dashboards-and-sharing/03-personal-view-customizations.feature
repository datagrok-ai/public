@guide @tutorial
Feature: How-to · Dashboards
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). The dashboard "Demographics overview" (demog with a scatter
  plot colored by SEX and a bar chart) is saved off camera through the ribbon's Save dialog (a project saved through
  the API keeps views with no table link, and its personal version then opens beside the original view as a second
  one), the Share offer cancelled. On camera the scatter plot
  is colored by RACE and saved with the Save project dialog's "Save personal view customizations", which keeps the
  views in the account's own settings (`projectCustomViews`, `projectCustomLayout`) and leaves the project as it
  was. Then every move is filmed: the view's tab is closed (the platform shows Browse again by itself), the dashboard
  opened again from Browse > Dashboards comes with the personal version and the platform's balloon about it (the
  check claims the balloon; the film does not show balloons, so the line does not mention it), its card
  in the Dashboards gallery and the status bar's Context Panel toggle show its Custom views pane, whose Reset drops
  the personal version, and closing the tab and opening the card again shows the shared one.

  The project, with its table, views and picture, is deleted when the feature starts and ends; the two user-setting
  maps are put back as they were at feature end (Reset already removes the entry on camera).

  Scenario: Change a dashboard just for yourself
    Personal view customizations: your own version of a shared dashboard, while everyone else keeps the original.

    Given user is logged in
    And simple mode is off
    And the "projectCustomViews" user settings are put back at feature end
    And the "projectCustomLayout" user settings are put back at feature end
    And no project named "Demographics overview" is on the server
    And user opens demog dataset
    And user adds a scatter plot viewer with:
      | xColumnName     | AGE    |
      | yColumnName     | WEIGHT |
      | colorColumnName | SEX    |
    And user adds a bar chart viewer with:
      | splitColumnName | RACE |
    And user clicks on Save ribbon item
    And user enters "Demographics overview" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    And the "Save project" dialog should close
    And "Share Demographics overview" dialog should be visible
    And user clicks on CANCEL button in "Share Demographics overview" dialog
    And the "Share Demographics overview" dialog should close
    And the notification balloons are closed
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    When the video says "A team dashboard. You want it a bit different, just for you"
    # caption: Color the points by race instead of sex
    And user picks "RACE" in the "color" column selector of scatter plot viewer
    Then "colorColumnName" property of scatter plot viewer should be "RACE"
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    # caption: Save personal view customizations: for you only
    When user hovers over "Save personal view customizations" text in "Save project" dialog
    And user selects "Save personal view customizations" in radio input in "Save project" dialog
    Then Name text input in "Save project" dialog should be disabled
    # caption: The dashboard itself stays as it was, for everyone else
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    # caption: Close it, and open it again later from Browse > Dashboards
    When user closes "demog" view tab
    Then the project "Demographics overview" should not be open
    And browse panel should be visible
    When user clicks on Dashboards tree node inside browse tree
    And user enters "Demographics overview" into gallery search
    Then "Demographics overview" gallery card should be visible
    # caption: It opens with your changes: colored by race
    When user double-clicks on "Demographics overview" gallery card and the project opens
    Then the "demog" view should be current
    And a warning balloon containing "Project has personal view customizations" should have been shown
    And the open table views should be exactly "demog"
    And "colorColumnName" property of scatter plot viewer should be "RACE"
    # caption: To go back to the shared version, select it in Dashboards
    When user clicks on the tab of the "Projects" view
    And user clicks on "Demographics overview" gallery card
    # caption: Its details are in the Context Panel
    And user clicks on "Context Panel" status bar toggle
    Then the context panel should show "Demographics overview"
    # caption: Custom views: your version of this dashboard
    When user expands "Custom views" pane in context panel
    Then "Custom views" pane in context panel should contain text "Project has personal view customizations"
    # caption: Reset drops your version
    When user clicks on Reset button in "Custom views" pane in context panel
    Then "Custom views" pane in context panel should not contain text "Project has personal view customizations"
    # caption: Close it and open it again: the shared version, colored by sex
    When user closes "demog" view tab
    Then the project "Demographics overview" should not be open
    When user double-clicks on "Demographics overview" gallery card and the project opens
    Then the "demog" view should be current
    And "colorColumnName" property of scatter plot viewer should be "SEX"
    And no errors should have been logged
