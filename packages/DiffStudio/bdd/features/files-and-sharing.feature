@journey @diffstudio @realizes:diffstudio.app.diff-studio
Feature: A model file previewed from Browse
  An .ivp file of the shared library opened straight from the Files tree, and the inputs it brings
  with it. Translated from files/TestTrack/DiffStudio/files-and-sharing.md and the spec beside it.

  The manual case shares the model by copying its address into a second tab; a feature is given one
  page, so it loads the address again from scratch instead — which is what pasting the link does,
  minus the second tab.

  Background:
    Given user is logged in
    And the browse panel is open

  Scenario: The library folder is reachable from the Files tree
    Given Files tree node inside browse tree is expanded
    And "Files > App Data" tree node inside browse tree is expanded
    And "Files > App Data > DiffStudio" tree node inside browse tree is expanded
    Then "Files > App Data > DiffStudio > library" tree node inside browse tree should be visible

  Scenario: A model file opens as a preview with its inputs
    When user clicks on "Files > App Data > DiffStudio > library" tree node inside browse tree
    And user clicks on pk.ivp link in gallery
    Then step input should be visible
    And count input should be visible
    And Multiaxis tab should be absent
    And Facet tab should be absent

  Scenario: The inputs the preview brings can be set
    When user enters "0.1" into step input
    Then step input should have value "0.1"
    When user hovers over count input
    And user clicks on plus icon in count input
    Then count input should not have value "1"

  Scenario: The address carries the inputs, and loading it again brings them back
    Then the page address should contain "step"
    When user opens the page address of the current view
    Then step input should have value "0.10"
    And no errors should have been logged
