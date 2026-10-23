@viewers @realizes:viewers.grid
Feature: Pin if editable pins a column for the users who may edit it
  `grid-ui.md` "Context Panel - Permissions": next to Editable by, Advanced > Permissions has a Pin if
  editable switch, which `grid-editing.feature` leaves out. With the switch on, the grid pins the
  column for a user Editable by names, and only for them: a column editable by someone else stays
  where it is, and the switch turned off unpins it again (`GridCore.applyPinIfEditable`). On
  demog-1000's SEVERITY, with the account the run signs in with and the library's second account
  (`bddsecond`, or DATAGROK_SHARING_LOGIN) as the other user; the tags are on the table opened for
  the feature, so nothing is left on the server.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And the context panel is open
    And the grid should pin the columns ""
    When user clicks on the "header SEVERITY" area of grid
    Then the context panel should show "SEVERITY"
    Given Advanced accordion header in context panel is expanded
    And Permissions accordion header in context panel is expanded

  Scenario: The switch pins the column for its editor, not for anyone else, and unpins it when off
    When user enters the sharing user's login into "Editable by" input in context panel
    And user checks "Pin if editable" input in context panel
    Then the grid should pin the columns ""
    When user enters the current user's login into "Editable by" input in context panel
    Then the grid should pin the columns "SEVERITY"
    When user unchecks "Pin if editable" input in context panel
    Then the grid should pin the columns ""
    And the "column order" reading of grid should contain "SEVERITY"
    When user clears "Editable by" input in context panel
    And user presses Enter in "Editable by" input in context panel
    Then no errors should have been logged
