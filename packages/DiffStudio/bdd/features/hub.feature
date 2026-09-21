@journey @diffstudio @realizes:diffstudio.app.diff-studio
Feature: The app's hub, and the Open model menu
  Apps > Diff Studio lands on the app's hub — templates and the library as cards — and every
  manual case begins by opening a model from the Open model icon of the ribbon, under Library.
  The other features go to a model by its address; this one walks the two ways in that a person
  uses, translated from step 1 of files/TestTrack/DiffStudio/open-model.md.

  Scenario: The app opens on its hub, and a library card opens a model
    Given user opens the Diff Studio app
    Then the "Diff Studio" view should be current
    And Create button should be visible
    When user double-clicks on Bioreactor hub card
    Then the "Bioreactor" view should be current
    And "Process mode" input should be visible

  Scenario: The Open model icon switches to another model of the library
    When user clicks on open model button
    And user picks "Library > PK-PD" from the open menu
    Then the "PK-PD" view should be current
    And dose input should be visible
    And "Process mode" input should be absent
    And no errors should have been logged
