@guide @tutorial @viewport:1280x720
Feature: How-to · Select and filter
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). Filter presets: with two cards filtering, Save or Apply >
  Save... in the filter panel's menu stores the panel's filters under a name, the panel's reset icon clears them,
  and the name under Save or Apply puts them back. demog (5850 rows: RACE Asian 72, of them F 37) is opened off
  camera with a filter panel holding the RACE and SEX cards. The preset is a server entity: it is deleted before
  the feature starts and when it ends.

  Scenario: Save filters and apply them later
    Save the filter panel's filters under a name, and bring them back in one click.

    Given user is logged in
    And simple mode is off
    And no saved filter state "Asian women" is kept, now or when the feature ends
    And user opens demog dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    And user opens an empty filter panel
    And user hovers over filter panel
    And user adds a card for "SEX" to the filter panel
    And user adds a card for "RACE" to the filter panel
    Then all rows should pass the filter
    # caption: Set the filters you want to keep
    When user clicks on the "category Asian of RACE" area of filter panel
    Then 72 rows should pass the filter
    When user clicks on the "category F of SEX" area of filter panel
    Then 37 rows should pass the filter
    # caption: Save or Apply > Save... keeps these filters under a name
    When user picks "Save or Apply | Save..." from the filter panel menu
    And user types "Asian women" into Name input in "Save filter preset" dialog
    And user clicks on OK button in "Save filter preset" dialog
    Then the "Save filter preset" dialog should close
    # caption: Reset clears the filters
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    # caption: Pick the name to bring them back
    When user picks "Save or Apply > Asian women" from the filter panel menu
    Then 37 rows should pass the filter
    And the "selected categories of RACE" reading of filter panel should be "Asian"
    And the "selected categories of SEX" reading of filter panel should be "F"
    When the video says "The preset is saved, so it is there next time too"
    Then no errors should have been logged
