@demo @realizes:u2.range-slider @realizes:u2.multi-select @realizes:u2.button-group @realizes:u2.combobox
Feature: Sliders, multi-select, button groups and comboboxes
  The inputs whose value is not one text: handles, tags, segmented buttons and popup lists.

  Scenario: Range slider handles answer the keyboard
    Given user opens the "Range slider" demo page
    Then value of range readout should have text "20 … 70"
    When user presses ArrowRight in Minimum slider handle in first range slider
    Then value of range readout should have text "21 … 70"
    When user presses ArrowLeft in Maximum slider handle in first range slider
    Then value of range readout should have text "21 … 69"
    When user presses ArrowRight in Minimum slider handle in second range slider
    Then value of dose readout should have text "10 … 1000 mg"

  Scenario: A multi-select shows its picks as tags
    Given user opens the "Multi-select" demo page
    Then value of targets readout should have text "Kinase"
    When user selects "Protease" in Targets input
    And user closes Targets input
    Then value of targets readout should have text "Kinase, Protease"
    When user clicks on "Remove Kinase" button in Targets input
    Then value of targets readout should have text "Protease"

  Scenario: Button groups as actions, a single toggle and a multi toggle
    Given user opens the "Multi-select" demo page
    Then Grid button should be selected
    When user clicks on Cards button
    Then value of layout readout should have text "cards"
    And Cards button should be selected
    And Grid button should not be selected
    When user clicks on Bold button
    And user clicks on Italic button
    Then value of style readout should have text "bold, italic"
    And Bold button should be selected
    When user clicks on Bold button
    Then value of style readout should have text "italic"

  Scenario: Comboboxes over local and async items
    Given user opens the "Async" demo page
    When user selects "Kyiv" in first combobox
    Then value of local readout should have text "Kyiv"
    When user types "Osl" into second combobox
    And user selects "Oslo" in second combobox
    Then value of remote readout should have text "Oslo"

  Scenario: An async view refreshes on its filter
    Given user opens the "Async" demo page
    Then async view should contain text "Berlin"
    When user types "Par" into "Filter the async view…" input
    Then async view should contain text "Paris"
    And async view should not contain text "Berlin"
