@guide @tutorial
Feature: How-to · Dashboards
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). The dashboard "Demographics overview" (demog with a scatter
  plot and a bar chart) is saved through the API and opened off camera. On camera a histogram is added and saved
  with the Save project dialog's default, "Save original project"; then a box plot is added and saved with "Save a
  copy" as "Demographics by race". Then the view's tab is closed (the platform shows Browse again by itself),
  Browse > Dashboards and a search list both, and the original opens with the histogram and without the box plot.
  (Closing the copy's view closes its project but leaves its demog table in the workspace, with no view: the
  original then opens beside it, so "no table should be left" is not claimed there.) Both projects, with their tables, views and pictures,
  are deleted when the feature starts and ends.

  Scenario: Save changes to the original or as a copy
    Save a changed dashboard over the original for everyone, or as a copy that leaves the original as it was.

    Given user is logged in
    And simple mode is off
    And no project named "Demographics overview" is on the server
    And no project named "Demographics by race" is on the server
    And user opens demog dataset
    And user adds a scatter plot viewer with:
      | xColumnName     | AGE    |
      | yColumnName     | WEIGHT |
      | colorColumnName | SEX    |
    And user adds a bar chart viewer with:
      | splitColumnName | RACE |
    And user saves the current view as project "Demographics overview"
    And user closes all views
    And user opens the "Demographics overview" project and waits for its table
    And the toolbox pane is shown
    And the context panel is closed
    And the browse panel is closed
    When the video says "This dashboard is already saved. Now change it"
    # caption: Add a histogram
    And user clicks on "histogram" icon in toolbox
    Then histogram viewer should be visible
    # caption: SAVE now asks where the change should go
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    # caption: Save original project: everyone who opens it gets the change
    When user hovers over "Save original project" text in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "Demographics overview" uploaded' should have been shown
    When the video says "For a variant, save a copy instead"
    # caption: Add a box plot
    And user clicks on "box plot" icon in toolbox
    Then box plot viewer should be visible
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    # caption: Save a copy: a new dashboard, the original stays as it was
    When user hovers over "Save a copy" text in "Save project" dialog
    And user selects "Save a copy" in radio input in "Save project" dialog
    Then Name text input in "Save project" dialog should have value "Copy of Demographics overview"
    # caption: Give the copy its own name
    When user enters "Demographics by race" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And an info balloon containing 'Project "Demographics by race" uploaded' should have been shown
    # caption: Close it. Browse > Dashboards now lists both
    When user closes "demog" view tab
    Then the project "Demographics by race" should not be open
    And browse panel should be visible
    When user clicks on Dashboards tree node inside browse tree
    And user enters "Demographics" into gallery search
    Then "Demographics overview" gallery card should be visible
    And "Demographics by race" gallery card should be visible
    # caption: The copy and the original, side by side
    When user hovers over "Demographics by race" gallery card
    # caption: Open the original
    When user double-clicks on "Demographics overview" gallery card and the project opens
    Then the "demog" view should be current
    And histogram viewer should be visible
    And box plot viewer should be absent
    When the video says "It has the histogram, but not the box plot"
    And user hovers over histogram viewer
    Then no errors should have been logged
