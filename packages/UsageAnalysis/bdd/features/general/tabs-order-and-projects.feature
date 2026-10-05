@serial
Feature: Datasets opened together keep their order through a project
  Two files dropped onto the window and two demo tables open as four table views without an error,
  and a project saved through the ribbon's Save dialog opens them again in the same order, with a
  table opened afterwards going to the end. Translated from the manual-only TestTrack case
  General/tabs-reordering-ui.md (tried at the user's request).

  Not translated, and why: reordering the tabs by drag and drop, returning a tab to its place, and a
  reordered strip surviving a new table, a closed table and the project round trip — no step reads
  the order of the view tabs (MISSING.md); the order claimed here is the order the views opened in,
  which a drag does not change. "Smooth", "snaps into place" and "no glitches" are not observable
  claims; the error floor stands for them.

  The project is named with the run's time and removed, with its tables, views and picture, at
  feature end and swept at its start.

  Background:
    Given user is logged in
    And the user's own project "BDD-Tabs-{time}" is removed now and at feature end

  Scenario: Four datasets, two of them dropped, open in four views without an error
    When user drops the "fixtures/browse-import.csv" file of the project onto status bar
    And user drops the "fixtures/cars-small.csv" file of the project onto status bar
    And user opens smiles dataset
    And user opens curves dataset
    Then the open table views should be exactly "browse-import, cars-small, smiles, curves"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The project saved through the ribbon opens its views in the same order
    When user drops the "fixtures/browse-import.csv" file of the project onto status bar
    And user drops the "fixtures/cars-small.csv" file of the project onto status bar
    And user opens smiles dataset
    And user opens curves dataset
    Then the open table views should be exactly "browse-import, cars-small, smiles, curves"
    When user opens the Save project dialog from the ribbon
    And user types "BDD-Tabs-{time}" into text input in "Save project" dialog
    And user clicks on OK in the Save project dialog and the project uploads
    Then the "Save project" dialog should close
    When user presses Escape
    Then 1 project named "BDD-Tabs-{time}" should be on the server
    When user closes all views
    Then no table should be left in the workspace
    When user opens the "BDD-Tabs-{time}" project and waits for its table
    Then the open table views should be exactly "browse-import, cars-small, smiles, curves"
    When user opens demog dataset
    Then the open table views should be exactly "browse-import, cars-small, smiles, curves, demog"
    And no errors should have been logged
