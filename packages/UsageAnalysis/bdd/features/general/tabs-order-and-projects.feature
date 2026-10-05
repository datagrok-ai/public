@serial
Feature: Dataset tabs reordered by drag and drop keep their order through a project
  Two files dropped onto the window and two demo tables open as four view tabs in that order; a tab
  dragged onto another moves next to it and can be dragged back, a dataset opened later goes to the end,
  and a project saved through the ribbon's Save dialog, reopened after a reload, brings the tabs back
  in the order they had when it was saved. Translated from the manual-only TestTrack case
  General/tabs-reordering-ui.md (tried at the user's request).

  The order is read from the tab strip itself: the shell is shown in full (simple mode off) with the
  Browse and Toolbox panes docked, whose tabs come first, then Home, so the first table's tab is the
  4th view tab. A tab dropped on the middle of another lands right after it (measured on the local
  stand). "Smooth", "snaps into place" and "no glitches" are not observable claims; the error
  floor stands for them.

  The project is named with the run's time and removed, with its tables, views and picture, at
  feature end and swept at its start.

  Background:
    Given user is logged in
    And simple mode is off
    And the browse panel is open
    And the toolbox pane is shown
    And the user's own project "BDD-Tabs-{time}" is removed now and at feature end
    When user drops the "fixtures/browse-import.csv" file of the project onto status bar
    Then the open table views should be exactly "browse-import"
    When user drops the "fixtures/cars-small.csv" file of the project onto status bar
    Then the open table views should be exactly "browse-import, cars-small"
    When user opens smiles dataset
    And user opens curves dataset
    Then the open table views should be exactly "browse-import, cars-small, smiles, curves"

  Scenario: Four datasets, two of them dropped, open as four tabs in that order
    Then 4th view tab should have text "browse-import"
    And 5th view tab should have text "cars-small"
    And 6th view tab should have text "smiles"
    And 7th view tab should have text "curves"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A tab dragged onto another moves next to it, and back again
    When user drags "curves" view tab to "browse-import" view tab
    Then 4th view tab should have text "browse-import"
    And 5th view tab should have text "curves"
    And 6th view tab should have text "cars-small"
    And 7th view tab should have text "smiles"
    When user drags "curves" view tab to "smiles" view tab
    Then 4th view tab should have text "browse-import"
    And 5th view tab should have text "cars-small"
    And 6th view tab should have text "smiles"
    And 7th view tab should have text "curves"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A dataset opened later goes to the end and leaves the new order alone
    When user drags "smiles" view tab to "browse-import" view tab
    Then 5th view tab should have text "smiles"
    And 6th view tab should have text "cars-small"
    When user opens demog dataset
    Then 4th view tab should have text "browse-import"
    And 5th view tab should have text "smiles"
    And 6th view tab should have text "cars-small"
    And 7th view tab should have text "curves"
    And 8th view tab should have text "demog"
    And no errors should have been logged

  Scenario: The reordered tabs come back in that order from the saved project
    When user drags "smiles" view tab to "browse-import" view tab
    Then 5th view tab should have text "smiles"
    And 6th view tab should have text "cars-small"
    When user opens the Save project dialog from the ribbon
    And user types "BDD-Tabs-{time}" into text input in "Save project" dialog
    And user clicks on OK in the Save project dialog and the project uploads
    Then the "Save project" dialog should close
    When user presses Escape
    Then 1 project named "BDD-Tabs-{time}" should be on the server
    And the "browse-import" table of the "BDD-Tabs-{time}" project should be saved as a snapshot
    And the "cars-small" table of the "BDD-Tabs-{time}" project should be saved as a snapshot
    And the "smiles" table of the "BDD-Tabs-{time}" project should be saved as a snapshot
    And the "curves" table of the "BDD-Tabs-{time}" project should be saved as a snapshot
    When user closes all views
    Then no table should be left in the workspace
    When user reloads the page
    And user opens the "BDD-Tabs-{time}" project and waits for its table
    Then 4th view tab should have text "browse-import"
    And 5th view tab should have text "smiles"
    And 6th view tab should have text "cars-small"
    And 7th view tab should have text "curves"
    And no errors should have been logged
