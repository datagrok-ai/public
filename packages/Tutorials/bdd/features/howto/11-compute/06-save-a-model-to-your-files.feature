@guide @tutorial
Feature: How-to · Compute
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). A model saved to the account's files and opened again: the
  Lotka-Volterra model of Diff Studio's library (Browse > Apps > Diff Studio), the ribbon's SAVE (Save to My files),
  the name typed, then the file found in Browse > Files > My files and clicked, which opens the model in Diff
  Studio's preview with its inputs (a double-click opens it in two views at once, "Lotka-Volterra (2)" and "(3)": a
  Diff Studio bug, reported, so the clip clicks once). Saving writes
  Predator-prey.ivp into the account's home folder and adds it to Diff Studio's recent models (a file in the same
  folder); the feature deletes the model file and puts the recent models back at its end (and sweeps a leftover
  model file at its start).

  Scenario: Save a model to your files
    Save a model, and it is a file in your Datagrok files that opens in Diff Studio again.

    Given user is logged in
    And simple mode is off
    And the "DiffStudio" package is installed
    And the package autostarts have completed
    And the file "diff-studio-recent.d42" of the user's home is put back at feature end
    And the file "Predator-prey.ivp" of the user's home is put back at feature end
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    # caption: Diff Studio is in Browse > Apps
    When user clicks on browse tab
    And user clicks on Apps tree node inside browse tree
    And user double-clicks on Diff-Studio gallery card
    Then the "Diff Studio" view should be current
    When user double-clicks on "Lotka-Volterra" model card
    Then "Prey" input should be visible
    # caption: SAVE keeps a copy of the model in your files
    When user clicks on Save button
    Then "Save to My files" dialog should be visible
    # caption: Give it a name
    When user enters "Predator-prey" into "Name" input in "Save to My files" dialog
    And user clicks on Save button in "Save to My files" dialog
    Then the "Save to My files" dialog should close
    # caption: It is in Browse > Files > My files
    When user expands Files tree node inside browse tree
    And user clicks on "Files > My files" tree node inside browse tree
    Then Predator-prey.ivp link in gallery should be visible
    # caption: Click it, and the model opens with its inputs
    When user clicks on Predator-prey.ivp link in gallery
    Then "Prey" input should be visible
    And no errors should have been logged
