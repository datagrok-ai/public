@journey @realizes:views.projects @realizes:core#239
Feature: A dashboard saves the arrangement of its views
  core#239: saving a dashboard records how its views are docked — side by side or stacked, with their
  proportions, and floating views as shares of the window — and reopening it puts them back. cars,
  iris and demog-1000 are arranged by dragging their tabs onto the dock compass: iris to the right of
  cars, cars narrowed to about 30% of the width the two share by its splitter, demog-1000 below cars.
  A markup view is dragged out to float. Views other than tables join a new dashboard only by choice:
  the Save project dialog lists them greyed out, its "+" includes one and a second is left out. The
  reopened dashboard brings back the arrangement (shares within 5 points, the float within 2% of the
  window) and only the markup view that was included. A markup view opened next to the open
  dashboard is then dragged onto it in the Dashboards panel, saved with it and reopened with it.

  Side panels (Browse, Toolbox, Context panel) are not part of the arrangement: the reopen happens
  with the Browse panel open, and the shares hold regardless. Markup views are opened through the API
  (no menu opens one). The dashboard is named with the run's time (letters and digits only: the gallery search misses "-") and removed with its tables, views
  and picture when the feature starts and ends.

  Background:
    Given user is logged in
    And simple mode is off
    And the user's own project "BDDLayout{time}" is removed now and at feature end

  Scenario: Tables docked side by side and below, a markup view floating
    Given user opens cars dataset
    And user opens iris dataset
    And user opens demog-1000 dataset
    When user drags the tab of the "iris" view to the right side of the "cars" view
    Then the "iris" view should be docked right-of the "cars" view
    When user resizes the "cars" view to 30% of the width it shares with the "iris" view
    And user drags the tab of the "demog-1000" view to the bottom side of the "cars" view
    Then the "demog-1000" view should be docked below the "cars" view
    And the "iris" view should be docked right-of the "demog-1000" view
    And the "cars" view should take about 30% of the width it shares with the "iris" view
    Given user opens a markup view named "BDD notes"
    When user drags the tab of the "BDD notes" view out to float
    Then the "BDD notes" view should be floating
    When user remembers where the "BDD notes" view floats
    Given user opens a markup view named "BDD scratch"
    Then no errors should have been logged

  Scenario: The Save dialog lists the markup views left out until one is included
    When user clicks on the tab of the "cars" view
    And user opens the Save project dialog from the ribbon
    Then the "BDD notes" view should be excluded in the Save project dialog
    And the "BDD scratch" view should be excluded in the Save project dialog
    When user includes the "BDD notes" view in the Save project dialog
    Then the "BDD scratch" view should be excluded in the Save project dialog
    When user types "BDDLayout{time}" into text input in "Save project" dialog
    And user clicks on OK in the Save project dialog and the project uploads
    Then the "Save project" dialog should close
    When user presses Escape
    Then 1 project named "BDDLayout{time}" should be on the server

  Scenario: The reopened dashboard brings the arrangement back
    When user closes all views
    Then no table should be left in the workspace
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDLayout{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDLayout{time} gallery card and the project opens
    Then the "iris" view should be docked right-of the "cars" view
    And the "cars" view should take about 30% of the width it shares with the "iris" view
    And the "demog-1000" view should be docked below the "cars" view
    And the "iris" view should be docked right-of the "demog-1000" view
    And the "BDD notes" view should float where it was remembered
    And "BDD scratch" view should be absent
    And no error or warning balloon should have been shown

  Scenario: A view dragged onto the open dashboard in the Dashboards panel is saved with it
    Given user opens a markup view named "BDD moved"
    When user clicks on Dashboards tab
    Then "New Dashboard > BDD moved" tree node inside browse tree should be visible
    When user drags "New Dashboard > BDD moved" tree node inside browse tree to BDDLayout{time} tree node inside browse tree
    Then Move entity dialog should be visible
    And Move entity dialog should contain text "BDD moved"
    When user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    Given BDDLayout{time} tree node inside browse tree is expanded
    Then "BDDLayout{time} > BDD moved" tree node inside browse tree should be visible
    When user clicks on Save button in BDDLayout{time} tree node inside browse tree
    Then "Save project" dialog should be visible
    And the "BDD moved" view should be included in the Save project dialog
    When user clicks on OK in the Save project dialog and the project uploads
    Then the "Save project" dialog should close
    When user closes all views
    Then no table should be left in the workspace
    When user opens the "BDDLayout{time}" project and waits until it has opened
    Then "BDD moved" view should be present
    And the "iris" view should be docked right-of the "cars" view
    And no errors should have been logged
