@demo @realizes:u2.dialog
Feature: Dialogs and confirm
  confirmDialog() asks with CANCEL and the action as a verb, and the work it starts keeps it open until done; a dialog
  of a fixed height keeps its footer at the bottom; the first focus passes read-only fields by; a required
  field says so.

  Background:
    Given user opens the "Dialogs" demo page

  Scenario: Delete asks first; CANCEL deletes nothing, DELETE deletes once its work is done
    When user clicks on "Delete campaign..." button
    Then "Delete campaign" dialog should be visible
    And "Delete campaign" dialog should contain text "This cannot be undone."
    When user clicks on CANCEL button in "Delete campaign" dialog
    Then "Delete campaign" dialog should be hidden
    And status bar should contain text "Delete cancelled"
    And value of deleted readout should have text "0"
    When user clicks on "Delete campaign..." button
    And user clicks on Delete button in "Delete campaign" dialog
    Then "Delete campaign" dialog should be hidden
    And status bar should contain text "Campaign KHL-3 deleted"
    And value of deleted readout should have text "1"

  Scenario: The first focus is CANCEL, never the action, so Enter deletes nothing
    When user clicks on "Delete campaign..." button
    Then CANCEL button in "Delete campaign" dialog should be focused
    When user presses Enter
    Then "Delete campaign" dialog should be hidden
    And status bar should contain text "Delete cancelled"
    And value of deleted readout should have text "0"

  Scenario: A fixed height: the footer stays at the bottom
    When user clicks on "Sized dialog" button
    Then "Progress tracker" dialog should be 420 px tall
    And body of "Progress tracker" dialog should be 339 px tall
    When user clicks on CANCEL button in "Progress tracker" dialog
    Then "Progress tracker" dialog should be hidden

  Scenario: The first focus goes to the first field that takes edits
    When user clicks on "Submit..." button
    Then Status input in Submit dialog should be focused
    And Name input in Submit dialog should be read-only
    When user clicks on CANCEL button in Submit dialog
    Then Submit dialog should be hidden

  Scenario: Required fields say so
    Then "template-name" text input should be required
    And "template-key" text input should be required
    And "template-description" text input should not be required
