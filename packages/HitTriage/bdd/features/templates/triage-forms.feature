@journey
Feature: Hit Triage and PepTriage forms: START and CREATE refuse what they cannot use
  Hit Triage and PepTriage keep their own new-campaign form and template editor. START and CREATE are disabled while
  the form is invalid, with the reason shown: START while a required field is empty; the editor shows its name error
  under the input, and its CREATE is disabled while the name or the key is empty or taken (PepTriage: or the sequence
  column is empty), the reason as its tooltip, so an empty name writes no templates/.json and an existing template's
  name overwrites nothing.

  A Hit Triage template "BDD HT Form" (a required String field Owner, data from a file) is made for the feature
  through the JS API (bindings/landing.ts). No scenario clicks START or CREATE: they write. The editor of Hit
  Triage is the parent of its landing (parentView), and grok.shell.closeAll then ends on the landing rather than Home,
  so the scenarios do not close every view between them. The editors open while the descriptors cannot be listed
  (template-editor.feature); they log that they could not (bdd.config.json allows that one error).

  Background:
    Given user is logged in
    And the descriptors cannot be listed, put back at feature end
    And the "Hit Triage" templates are on the server:
      | name               | key          | fields                 |
      | BDD HT Form {time} | BDDHTF{time} | Owner: String required |
    And user opens the "Hit Triage" app

  Scenario: While a required field is empty, START is disabled
    When user selects "BDD HT Form {time}" in Template input
    Then "Owner" input should be visible
    And Start button should be disabled

  Scenario: The editor shows its error while the name is empty or an existing template's
    When user clicks on "New Template" icon
    Then the "New Template" view should be current
    When user types "BDD HT Form {time}" into template name input
    Then the template editor should show the error "Template name is empty or already exists"
    When user types "BDD HT Free {time}" into template name input
    Then the template editor should not show the error "Template name is empty or already exists"
    When user clears template name input
    Then the template editor should show the error "Template name is empty or already exists"
    When user clicks on Cancel button
    Then the "Hit Triage" view should be current

  Scenario: CREATE is disabled while the name is empty or an existing template's
    When user clicks on "New Template" icon
    Then the "New Template" view should be current
    And template name input should have the value ""
    And Create button should be disabled
    When user types "BDD HT Form {time}" into template name input
    Then the template editor should show the error "Template name is empty or already exists"
    And Create button should be disabled
    When user clicks on Cancel button
    Then the "Hit Triage" view should be current

  Scenario: In PepTriage, CREATE is disabled while the name is empty
    Given the "Bio" package is installed
    When user opens the "PepTriage" app
    And user clicks on "New Template" icon
    Then the "New Template" view should be current
    And template name input should have the value ""
    And Create button should be disabled
    When user clicks on Cancel button
    Then the "PepTriage" view should be current
