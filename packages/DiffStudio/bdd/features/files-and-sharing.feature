@journey @diffstudio @realizes:diffstudio.app.diff-studio
Feature: A model file previewed from Browse
  An .ivp file of the shared library opened straight from the Files tree, and the inputs it brings
  with it. Translated from files/TestTrack/DiffStudio/files-and-sharing.md and the spec beside it.

  The manual case shares the model by copying its address into a second tab; a feature is given one
  page, so it loads the address again from scratch instead — which is what pasting the link does,
  minus the second tab.

  Step is set the way the case sets it, by dragging its slider: the track spans 0.01 to 0.1 in steps
  of 0.0009, so the drag lands within a step of the value asked for rather than on it exactly.

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

  Scenario: The slider sets Step, as a reader would set it
    When user takes a snapshot of line chart viewer
    And user drags the slider of step input to 0.1
    Then step input should have a value between 0.09 and 0.1
    And line chart viewer should have repainted

  Scenario: The clicker counts Count up to four
    When user hovers over count input
    And user clicks on plus icon in count input
    And user clicks on plus icon in count input
    And user clicks on plus icon in count input
    Then count input should have value "4"

  Scenario: The inputs take a typed value too
    When user enters "0.1" into step input
    Then step input should have value "0.1"

  Scenario: The address carries the inputs, and loading it again brings them back
    Then the page address should contain "step"
    When user opens the model at the page address
    Then step input should have value "0.10"
    And no errors should have been logged
