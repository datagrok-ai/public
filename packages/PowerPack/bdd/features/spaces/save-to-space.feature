@journey @realizes:file.menu.save.tables-as-project @realizes:views.space
Feature: The Space row of the Save project dialog puts a project into a space
  The Space row (core#1507): a new project goes into the picked space, a saved one moves or is linked,
  a name the space holds keeps OK off, and the "+" buttons make nested spaces. Fixtures go at start and end.

  Background:
    Given user is logged in
    And simple mode is off
    And the browse panel is closed
    And a space named "BDDSaveSpace{time}" is on the server
    And a space named "BDDSaveLink{time}" is on the server
    And a space named "BDDSaveMove{time}" is on the server
    And no project named "BDDSpaceProj{time}, BDDNestedProj{time}" is on the server
    And no space named "BDDSaveInner{time}, BDDSaveTop{time}" is on the server

  Scenario: The Space row says what an empty choice means
    Given user opens demog-1000 dataset
    When user opens the Save project dialog from the ribbon
    Then "Save project" dialog should be visible
    And space input in "Save project" dialog should contain text "Save to a space (optional)"
    When user hovers over space input in "Save project" dialog
    Then tooltip should contain text "Empty: saved to your own projects"

  Scenario: A new project picked into a space is saved straight there
    When user enters "BDDSpaceProj{time}" into Name text input in "Save project" dialog
    And user clicks on space input in "Save project" dialog
    And user clicks on "BDDSaveSpace{time}" tree node
    Then space input in "Save project" dialog should contain text "BDDSaveSpace{time}"
    And "Save project" dialog should contain text "Shared with whoever can see \"BDDSaveSpace{time}\"."
    And "Move" choice chip in "Save project" dialog should be hidden
    When user clicks on OK in the Save project dialog and the project uploads
    Then the project named "BDDSpaceProj{time}" should be in the space "BDDSaveSpace{time}" on the server

  Scenario: Re-saved with Link, the project stays and the other space links it
    When user opens the Save project dialog from the ribbon
    Then space input in "Save project" dialog should contain text "BDDSaveSpace{time}"
    And "Save project" dialog should contain text "Lives in \"BDDSaveSpace{time}\"."
    When user clicks on space input in "Save project" dialog
    And user clicks on "BDDSaveLink{time}" tree node
    Then "Save project" dialog should contain text "Moves from \"BDDSaveSpace{time}\" to \"BDDSaveLink{time}\"."
    When user clicks on "Move" choice chip in "Save project" dialog
    And user clicks on "Link" menu item
    Then "Save project" dialog should contain text "Stays in \"BDDSaveSpace{time}\"; \"BDDSaveLink{time}\" gets a link to it."
    When user clicks on OK in the Save project dialog and the project uploads
    Then the project named "BDDSpaceProj{time}" should be in the space "BDDSaveSpace{time}" on the server
    And the project named "BDDSpaceProj{time}" should be linked from the space "BDDSaveLink{time}" on the server

  Scenario: Moved into the space that links it, the project makes the link its home
    When user opens the Save project dialog from the ribbon
    And user clicks on space input in "Save project" dialog
    And user clicks on "BDDSaveLink{time}" tree node
    Then "Save project" dialog should contain text "Moves from \"BDDSaveSpace{time}\" to \"BDDSaveLink{time}\"."
    And "Save project" dialog should not contain text "already has a project named"
    And OK button in "Save project" dialog should be enabled
    When user clicks on OK in the Save project dialog and the project uploads
    Then the project named "BDDSpaceProj{time}" should be in the space "BDDSaveLink{time}" on the server
    And the project named "BDDSpaceProj{time}" should not be in the space "BDDSaveSpace{time}" on the server

  Scenario: Re-saved with Move, the project goes to the other space
    When user opens the Save project dialog from the ribbon
    And user clicks on space input in "Save project" dialog
    And user clicks on "BDDSaveMove{time}" tree node
    Then "Save project" dialog should contain text "Moves from \"BDDSaveLink{time}\" to \"BDDSaveMove{time}\"."
    When user clicks on OK in the Save project dialog and the project uploads
    Then the project named "BDDSpaceProj{time}" should be in the space "BDDSaveMove{time}" on the server
    And the project named "BDDSpaceProj{time}" should not be in the space "BDDSaveLink{time}" on the server

  Scenario: A space that holds a project of the name holds OK off
    When user opens the Save project dialog from the ribbon
    And user selects "Save a copy" in radio input in "Save project" dialog
    And user enters "BDDSpaceProj{time}" into Name text input in "Save project" dialog
    And user clicks on space input in "Save project" dialog
    And user clicks on "BDDSaveMove{time}" tree node
    Then "Save project" dialog should contain text "\"BDDSaveMove{time}\" already has a project named \"BDDSpaceProj{time}\""
    And OK button in "Save project" dialog should be disabled
    When user enters "BDDSpaceProjCopy{time}" into Name text input in "Save project" dialog
    Then "Save project" dialog should not contain text "already has a project named"
    And OK button in "Save project" dialog should be enabled
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close

  Scenario: Spaces made from the "+" buttons nest, and the project is saved into the inner one
    When user opens the Save project dialog from the ribbon
    And user selects "Save a copy" in radio input in "Save project" dialog
    And user enters "BDDNestedProj{time}" into Name text input in "Save project" dialog
    And user clicks on space input in "Save project" dialog
    And user clicks on "New top-level space" button
    And user types "BDDSaveTop{time}" into "Rename New space" text input
    And user presses Enter
    Then space input in "Save project" dialog should contain text "BDDSaveTop{time}"
    When user clicks on "New space inside BDDSaveTop{time}" button
    And user types "BDDSaveInner{time}" into "Rename New space" text input
    And user presses Enter
    Then space input in "Save project" dialog should contain text "BDDSaveTop{time} / BDDSaveInner{time}"
    When user clicks on space input in "Save project" dialog
    And user clicks on OK in the Save project dialog and the project uploads
    Then the project named "BDDNestedProj{time}" should be in the space "BDDSaveInner{time}" on the server
