@journey @diffstudio @realizes:diffstudio.app.diff-studio
Feature: A model saved to the library
  Saving a model from Diff Studio into the shared library and finding it there again. Translated
  from files/TestTrack/DiffStudio/catalog.md and the spec beside it.

  Saving writes a real file into System:AppData/DiffStudio/library, so the feature deletes exactly
  what it added when it ends — the old suite did not, which is why the stand's library carries
  PK-PD(1)…PK-PD(14) from earlier runs.

  The platform announces the change on its own event bus (diff-studio:library-changed), so the save
  is claimed by that event rather than by the balloon the old spec watched for and had to treat as
  best-effort because it auto-dismisses.

  Background:
    Given user is logged in
    And user opens the "PK-PD" model of the Diff Studio library

  Scenario: The model is the one the library serves
    Then the "PK-PD" view should be current
    And dose input should be visible

  Scenario: The model answers to its inputs
    When user clicks on Multiaxis tab
    And user takes a snapshot of line chart viewer
    And user enters "5000" into dose input
    Then dose input should have value "5000"
    And line chart viewer should have repainted

  Scenario: Saving it to the library is announced
    Given user listens for "diff-studio:library-changed" custom event
    When user saves the model to the Diff Studio library
    Then the "diff-studio:library-changed" custom event should have fired

  Scenario: The Model Hub lists the saved model
    Given user opens the Model Hub
    Then PK-PD link in gallery should be visible

  Scenario: The model runs from the catalog and its chart follows its inputs
    When user double-clicks on PK-PD link in gallery
    Then dose input should be visible
    When user takes a picture of viewer
    And user enters "5000" into dose input
    Then dose input should have value "5000"
    And viewer should look different
    And no errors should have been logged
