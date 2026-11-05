@journey
Feature: Hit Design new campaign: the form a template builds, its required fields, and Start
  Choosing a template in the New Campaign section builds the form below it: Campaign Name, one input per field of
  the template, and START. START with every required field filled makes the campaign <KEY>-<n>, puts its code in the
  address, opens it on one row (Molecule, Stage at the template's first stage, V-iD) with the molecule sketcher on
  that row, and saves it at once under the typed name. Every scenario but the first starts on a freshly opened app.

  Two templates are made for the feature through the JS API (bindings/landing.ts): "BDD New Campaign", with a
  required String field Owner, a Number field Budget and the stages Idea and Made; and "BDD New Layout", whose layout
  has a bar chart next to the grid. The campaigns their keys start go with them at feature end, with their rows in
  the hitdesign database. Needs Chem on the stand (the sketcher).

  START is disabled while a required field is empty, the reason its tooltip. A new campaign from a template with a
  layout opens with the template's layout, not the grid alone.

  Background:
    Given user is logged in
    And the "Chem" package is installed
    And the "Hit Design" templates are on the server:
      | name                    | key          | fields                                 | stages    | layout             |
      | BDD New Campaign {time} | BDDNEW{time} | Owner: String required; Budget: Number | Idea, Made |                    |
      | BDD New Layout {time}   | BDDNLY{time} | Layout note: String                    | Idea      | bar chart of Stage |
    And user opens the "Hit Design" app

  Scenario: A template's form has Campaign Name, the template's fields and START
    When user selects "BDD New Campaign {time}" in Template input
    Then "Campaign Name" input should be visible
    And "Owner" input should be visible
    And "Budget" input should be visible
    And Start button should be visible

  Scenario: While a required field is empty, START is disabled
    When user closes all views
    And user opens the "Hit Design" app
    And user selects "BDD New Campaign {time}" in Template input
    Then "Owner" input should be visible
    And Start button should be disabled
    When user types "Ann" into "Owner" input
    Then Start button should be enabled

  Scenario: START makes the campaign under its code and name and opens it on one row with the sketcher
    When user closes all views
    And user opens the "Hit Design" app
    And user selects "BDD New Campaign {time}" in Template input
    And user types "Ann" into "Owner" input
    And user types "BDD New Run {time}" into "Campaign Name" input
    And user clicks on Start button
    Then the "BDD New Run {time}" view should be current
    And the page address should contain "campaignId=BDDNEW{time}-1"
    And sketcher dialog should be visible
    And the table should have the columns "Molecule, Stage, V-iD"
    And the table should have 1 row
    And the value of "Stage" column in row 1 should be "Idea"
    And the "Hit Design" campaign "BDDNEW{time}-1" should be on the server
    And the "Hit Design" campaign "BDDNEW{time}-1" on the server should have the name "BDD New Run {time}"
    And the "Hit Design" campaign "BDDNEW{time}-1" on the server should have the status "In Progress"
    When user clicks on CANCEL button in sketcher dialog
    Then sketcher dialog should be absent

  Scenario: A new campaign from a template with a layout opens with that layout
    When user closes all views
    And user opens the "Hit Design" app
    And user selects "BDD New Layout {time}" in Template input
    Then "Layout note" input should be visible
    When user clicks on Start button
    Then the "BDDNLY{time}-1" view should be current
    And sketcher dialog should be visible
    When user clicks on CANCEL button in sketcher dialog
    Then bar chart viewer should be visible
