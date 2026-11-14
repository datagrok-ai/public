@journey @serial @viewers @realizes:viewers.filters
Feature: Filter presets in the Save or Apply menu
  A filter preset is a named entity on the server that can be shared. Save or Apply lists the ones that
  fit the table, newest first, with who saved each beside its name: 20 at once, the next ones on Show
  more, any of them through the menu's search box, in the order Sort by picks. Hovering a preset offers
  to share it and to delete it when one may; a second account it was shared with to view finds it in
  the same menu under its author's name, applies it, and is offered neither. One journey on demog-1000, where RACE is Asian in 15 rows, as the
  running account and the sharing user. The 24 presets "bdd bulk #" on SEX are made through the server
  one after another, so the menu holds more than 20; the one saved through the dialog is the newest.
  Not translated here: which tables are offered a preset (filter-panel-saved-states.feature, and the
  datlas suite domain_saved_filters_test with the server side of the sort keys and of sharing).

  Serial: it signs in as the sharing user on its page.

  Background:
    Given user is logged in
    And no saved filter state "bdd shared preset" is kept, now or when the feature ends
    And 24 saved filter presets "bdd bulk #" on the text column "SEX" are kept until the feature ends
    And user opens demog-1000 dataset
    And user opens an empty filter panel

  Scenario: A preset saved with a description is listed first, as one's own
    When user adds a card for "RACE" to the filter panel
    And user clicks on the "category Asian of RACE" area of filter panel
    Then 15 rows should pass the filter
    When user picks "Save or Apply | Save..." from the filter panel menu
    And user types "bdd shared preset" into Name input in "Save filter preset" dialog
    And user types "Asian patients only" into Description input in "Save filter preset" dialog
    And user clicks on OK button in "Save filter preset" dialog
    Then the "Save filter preset" dialog should close
    When user opens the viewer menu of filter panel
    Then the open menu should list "Save or Apply > bdd shared preset" above "bdd bulk 24"
    And the filter preset "bdd shared preset" in the open menu should be by "me"
    And the filter preset "bdd shared preset" in the open menu should offer "share and delete"
    When user closes the context menu

  Scenario: The menu shows the first 20 presets, the next ones on scrolling to the end, and finds any by search
    When user opens the viewer menu of filter panel
    Then the open menu should list "Save or Apply > bdd bulk 6"
    And the open menu should not list "Save or Apply > bdd bulk 5"
    When user scrolls the presets under Save or Apply to the end in the open menu
    Then the open menu should list "Save or Apply > bdd bulk 5"
    And the open menu should list "Save or Apply > bdd bulk 1"
    When user closes the context menu
    And user opens the viewer menu of filter panel
    Then the open menu should not list "Save or Apply > bdd bulk 4"
    When user types "bulk 4" into the search box of "Save or Apply" in the open menu
    Then the open menu should list "Save or Apply > bdd bulk 4"
    When user closes the context menu

  Scenario: An older preset applied comes first under Sort by Used by me
    When user opens the viewer menu of filter panel
    And user types "bulk 3" into the search box of "Save or Apply" in the open menu
    And user picks "Save or Apply > bdd bulk 3" from the open menu
    Then the filter panel should have 2 filters
    When user opens the viewer menu of filter panel
    Then the open menu should list "Save or Apply > bdd shared preset" above "bdd bulk 24"
    When user picks "Save or Apply > Sort by > Used by me" from the open menu
    Then the open menu should list "Save or Apply > bdd bulk 3" above "bdd shared preset"
    When user picks "Save or Apply > Sort by > Updated" from the open menu
    Then the open menu should not list "Save or Apply > bdd bulk 3"
    When user closes the context menu

  Scenario: A preset shared to view is applied by the second account, which can neither share nor delete it
    When user opens the viewer menu of filter panel
    And user clicks the share icon of the filter preset "bdd shared preset" in the open menu
    Then "Share bdd shared preset" dialog should be visible
    When user picks the sharing user in "User, group, or email" input in "Share bdd shared preset" dialog
    And user unchecks "Send notifications" input in "Share bdd shared preset" dialog
    And user clicks on OK button in "Share bdd shared preset" dialog
    Then the "Share bdd shared preset" dialog should close
    When user signs in as the sharing user
    Then the sharing user should be signed in
    Given user opens demog-1000 dataset
    And user opens an empty filter panel
    When user opens the viewer menu of filter panel
    Then the open menu should not list "Save or Apply > bdd bulk 24"
    And the filter preset "bdd shared preset" in the open menu should name its author
    And the filter preset "bdd shared preset" in the open menu should offer "nothing"
    When user picks "Save or Apply > bdd shared preset" from the open menu
    Then 15 rows should pass the filter
    When user signs in as themselves again
    Then the running account should be signed in

  Scenario: Deleting one's own preset asks first, then takes it out of the menu
    Given user opens demog-1000 dataset
    And user opens an empty filter panel
    When user opens the viewer menu of filter panel
    And user clicks the delete icon of the filter preset "bdd shared preset" in the open menu
    Then "Are you sure?" dialog should be visible
    And "Are you sure?" dialog should contain text "deleted for everyone it is shared with"
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    When user opens the viewer menu of filter panel
    Then the open menu should list "Save or Apply > bdd bulk 24"
    And the open menu should not list "Save or Apply > bdd shared preset"
    When user closes the context menu
    Then no errors should have been logged
    And no error or warning balloon should have been shown
