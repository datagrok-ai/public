@journey
Feature: Hit Design template picker: choosing, adding, cloning and deleting templates
  The Template picker of the New Campaign form lists the app's templates by file name and carries three icons beside
  it: New Template (+), Clone template and Delete template. A choice rebuilds the form below it from the template;
  "+" and Clone open the New Template editor; Delete asks before it removes the template file, then moves the picker
  to another template. Every scenario starts on a freshly opened app.

  Two templates are made for the feature through the JS API (bindings/landing.ts): "BDD Picker Source", with a
  required field Owner, the stages Design and Make, the Demo File Submit function and a layout with a bar chart; and
  "BDD Picker Doomed", which the last scenario deletes. A clone saved as "BDD Picker Clone" goes at feature end.
  Clone takes the template the form was last built for, so it is pressed once the form shows Owner.

  A clone proposes "<name> (copy)" and the key with the next free number, so CREATE is enabled at once; its Submit
  function shows the function by its friendly name; it carries the template's layout. The editor opens while the
  descriptors cannot be listed (template-editor.feature).

  Background:
    Given user is logged in
    And the descriptors cannot be listed, put back at feature end
    And the "Hit Design" templates are on the server:
      | name                     | key          | fields                 | stages       | submit           | layout             |
      | BDD Picker Source {time} | BDDPKS{time} | Owner: String required | Design, Make | Demo File Submit | bar chart of Stage |
      | BDD Picker Doomed {time} | BDDPKD{time} |                        | Design       |                  |                    |
    And the "Hit Design" templates named "BDD Picker Clone {time}" are removed now and at feature end
    And user opens the "Hit Design" app

  Scenario: The picker lists the templates, offers New, Clone and Delete, and a choice builds its form
    Then Template input should offer the choice "BDD Picker Source {time}"
    And Template input should offer the choice "BDD Picker Doomed {time}"
    And "New Template" icon should be visible
    And "Clone template" icon should be visible
    And "Delete template" icon should be visible
    When user selects "BDD Picker Source {time}" in Template input
    Then "Campaign Name" input should be visible
    And "Owner" input should be visible
    And Start button should be visible

  Scenario: The "+" icon opens an empty New Template editor, and CANCEL comes back
    When user closes all views
    And user opens the "Hit Design" app
    And user clicks on "New Template" icon
    Then the "New Template" view should be current
    And the page address should contain "/new-template"
    And template name input should have the value ""
    And template key input should have the value ""
    When user clicks on Cancel button
    Then the "Hit Design" view should be current

  Scenario: Clone fills the editor with the selected template's fields and stages
    When user closes all views
    And user opens the "Hit Design" app
    And user selects "BDD Picker Source {time}" in Template input
    Then "Owner" input should be visible
    When user clicks on "Clone template" icon
    Then the "New Template" view should be current
    And the template editor should list the fields "Owner (String, required)"
    And the template editor should list the stages "Design, Make"
    When user clicks on Cancel button
    Then the "Hit Design" view should be current

  Scenario: A clone proposes a name and a key no template has
    When user closes all views
    And user opens the "Hit Design" app
    And user selects "BDD Picker Source {time}" in Template input
    Then "Owner" input should be visible
    When user clicks on "Clone template" icon
    Then the "New Template" view should be current
    And template name input should have the value "BDD Picker Source {time} (copy)"
    And template key input should not have the value "BDDPKS{time}"
    And Create button should be enabled
    When user clicks on Cancel button
    Then the "Hit Design" view should be current

  Scenario: A clone keeps the template's submit function
    When user closes all views
    And user opens the "Hit Design" app
    And user selects "BDD Picker Source {time}" in Template input
    Then "Owner" input should be visible
    When user clicks on "Clone template" icon
    Then the "New Template" view should be current
    And "Submit function" input should have the value "Demo File Submit"
    When user clicks on Cancel button
    Then the "Hit Design" view should be current

  Scenario: A clone saved under a new name carries the template's layout, and CREATE comes back with it selected
    When user closes all views
    And user opens the "Hit Design" app
    And user selects "BDD Picker Source {time}" in Template input
    Then "Owner" input should be visible
    When user clicks on "Clone template" icon
    Then the "New Template" view should be current
    When user types "BDD Picker Clone {time}" into template name input
    And user types "BDDPKC{time}" into template key input
    And user clicks on Create button
    Then the "Hit Design" view should be current
    And Template input should have the value "BDD Picker Clone {time}"
    And the "Hit Design" template "BDD Picker Clone {time}" on the server should have the fields "Owner (String, required)"
    And the "Hit Design" template "BDD Picker Clone {time}" on the server should have the layout "a layout"

  Scenario: Delete asks for a confirmation, removes the template and moves the picker to another one, with no error
    When user closes all views
    And user opens the "Hit Design" app
    And user selects "BDD Picker Doomed {time}" in Template input
    And user clicks on "Delete template" icon
    Then "Delete template" dialog should contain text "Are you sure you want to delete template"
    And "Delete template" dialog should contain text "BDD Picker Doomed {time}"
    When user confirms the "Delete template" dialog
    Then the "Hit Design" template "BDD Picker Doomed {time}" should not be on the server
    And an info balloon should have been shown
    And no error or warning balloon should have been shown
    And Template input should not offer the choice "BDD Picker Doomed {time}"
    And Template input should not have the value ""
