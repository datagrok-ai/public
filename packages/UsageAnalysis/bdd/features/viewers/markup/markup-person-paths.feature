@journey @viewers @realizes:viewers.markup
Feature: Markup — the keyboard in the grid moves it, Show Title hides its title, Edit content applies
  Three things a person does to a markup viewer that the property-driven features do not: walk the
  grid with the arrow keys and watch the `${COLUMN}` references follow the current row, switch
  Show Title off while a title is set, and change the content through **Edit content...** and OK.
  `markup-row-binding.feature` moves the current row through the API; here the row moves because a
  grid cell was clicked and ArrowDown pressed. demog-1000: row 1 is AGE 26, SEX F, WEIGHT 74.10;
  row 2 is AGE 30, SEX F.

  The Edit dialog's OK is clicked here: `markup-chrome.feature` dismisses the dialog with Escape
  because on a stand with no internet the help pane's preloader covered OK; if the preloader comes
  back on a stand, this scenario is where it shows.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a markup viewer with:
      | content | Age: ${AGE} Sex: ${SEX} |
    Then markup viewer should report no error

  Scenario: ArrowDown in the grid moves the current row, and the markup follows it
    When user clicks on the "cell 1 of AGE" area of grid
    Then the "current row" reading of markup viewer should be 1
    And the "text" reading of markup viewer should be "Age: 26 Sex: F"
    When user presses ArrowDown
    Then the "current row" reading of markup viewer should be 2
    And the "text" reading of markup viewer should be "Age: 30 Sex: F"
    When user presses ArrowUp
    Then the "current row" reading of markup viewer should be 1
    And the "text" reading of markup viewer should be "Age: 26 Sex: F"
    And no errors should have been logged

  Scenario: Show Title off hides a title that is still set
    When user sets properties of markup viewer:
      | Show Title | true         |
      | Title      | Patient card |
    Then title of markup viewer should have text "Patient card"
    When user sets "Show Title" property of markup viewer to "false"
    Then title of markup viewer should be hidden
    When user sets "Show Title" property of markup viewer to "true"
    Then title of markup viewer should have text "Patient card"
    When user sets properties of markup viewer:
      | Title      |       |
      | Show Title | false |
    Then no errors should have been logged

  Scenario: Edit content... and OK put the new content on screen
    Given user makes row 1 current
    Then the "text" reading of markup viewer should be "Age: 26 Sex: F"
    When user picks "Edit content..." from the context menu of markup viewer
    Then Edit dialog should be visible
    When user enters "Weight: ${WEIGHT}" into input in Edit dialog
    And user clicks on "OK" button in Edit dialog
    Then Edit dialog should be absent
    And "content" property of markup viewer should be "Weight: ${WEIGHT}"
    And the "text" reading of markup viewer should be "Weight: 74.10"
    And no errors should have been logged
