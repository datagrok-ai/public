@journey
Feature: Dataset tabs reordered by drag and drop keep their order through a project
  Two files dropped onto the window and two demo tables open as four view tabs in that order; a tab
  dragged onto another moves next to it and can be dragged back, a dataset opened later goes to the end,
  and a project saved through the ribbon's Save dialog and reopened brings the tabs back in the order
  they had when it was saved. Translated from the manual-only TestTrack case
  General/tabs-reordering-ui.md (tried at the user's request).

  The order is read from the tab strip of the full shell (simple mode off), counting only the tabs of
  the views named, so the Home, Browse and Toolbox tabs do not shift it. A tab dropped on the middle of
  another lands right after it. "Smooth", "snaps into place" and "no glitches" are not observable
  claims; the error floor stands for them. How the Save dialog stores the tables is the projects
  features' claim, not this one's.

  The project is named with the run's time and removed, with its tables, views and picture, at
  feature end and swept at its start.

  Background:
    Given user is logged in
    And simple mode is off
    And the user's own project "BDD-Tabs-{time}" is removed now and at feature end
    When user drops the "fixtures/browse-import.csv" file of the project onto status bar
    Then the open table views should be exactly "browse-import"
    When user drops the "fixtures/cars-small.csv" file of the project onto status bar
    Then the open table views should be exactly "browse-import, cars-small"
    When user opens smiles dataset
    And user opens curves dataset
    Then the open table views should be exactly "browse-import, cars-small, smiles, curves"

  Scenario: Four datasets, two of them dropped, open as four tabs in that order
    Then the view tabs should be in the order "browse-import, cars-small, smiles, curves"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A tab dragged onto another moves next to it, and back again
    When user drags "curves" view tab to "browse-import" view tab
    Then the view tabs should be in the order "browse-import, curves, cars-small, smiles"
    When user drags "curves" view tab to "smiles" view tab
    Then the view tabs should be in the order "browse-import, cars-small, smiles, curves"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A dataset opened later goes to the end and leaves the new order alone
    When user drags "smiles" view tab to "browse-import" view tab
    Then the view tabs should be in the order "browse-import, smiles, cars-small, curves"
    When user opens demog dataset
    Then the view tabs should be in the order "browse-import, smiles, cars-small, curves, demog"
    When user closes the current view
    Then the open table views should be exactly "browse-import, cars-small, smiles, curves"
    And the view tabs should be in the order "browse-import, smiles, cars-small, curves"
    And no errors should have been logged

  Scenario: The reordered tabs come back in that order from the saved project
    When user opens the Save project dialog from the ribbon
    And user types "BDD-Tabs-{time}" into text input in "Save project" dialog
    And user clicks on OK in the Save project dialog and the project uploads
    Then the "Save project" dialog should close
    When user presses Escape
    Then 1 project named "BDD-Tabs-{time}" should be on the server
    When user closes all views
    Then no table should be left in the workspace
    When user opens the "BDD-Tabs-{time}" project and waits for its table
    Then the view tabs should be in the order "browse-import, smiles, cars-small, curves"
    And no errors should have been logged
