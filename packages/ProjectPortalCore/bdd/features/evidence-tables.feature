@portal
Feature: Evidence tables generated into documents
  For a document criterion, the evidence behind it (the run selected for it with its inputs and
  outputs, the newest snapshot, the gate's data conditions, the decisions citing it) is generated
  into a table: a CSV and a Markdown file in the program's folder, linked as a new draft version of
  a document over the criterion, which records the evidence it was generated from.

  Background:
    Given user opens the Project Portal app

  Scenario: The process description's evidence is generated into a draft document
    Given a scratch project "BDD-ET" at stage "D4" in the "EX-BIOC" portfolio
    And the workspace of the scratch project "BDD-ET" in the "EX-BIOC" portfolio is built
    And a run of the example model at feed rate 2 titled "Feed 2" saved to the scratch project "BDD-ET" in the "EX-BIOC" portfolio
    And a run of the example model at feed rate 4 titled "Feed 4" saved to the scratch project "BDD-ET" in the "EX-BIOC" portfolio
    When user clicks on biologics development portfolio
    And user clicks on scratch evidence project card
    And user clicks on choose process description run button
    Then choose run dialog should be visible
    When user selects "Feed 4" in selected run input
    And user enters "Higher titer at the same temperature" into selection rationale input
    And user clicks on select run button
    Then choose run dialog should not be visible
    When user clicks on generate process description evidence button
    Then process description document card should contain text "Evidence table, Locked process description · v1 · draft"
    And process description document card should contain text "Table"

  Scenario: Another run is selected, the table goes stale, and is generated again
    Given a scratch project "BDD-ET2" at stage "D4" in the "EX-BIOC" portfolio
    And the workspace of the scratch project "BDD-ET2" in the "EX-BIOC" portfolio is built
    And a run of the example model at feed rate 2 titled "Feed 2" saved to the scratch project "BDD-ET2" in the "EX-BIOC" portfolio
    And a run of the example model at feed rate 4 titled "Feed 4" saved to the scratch project "BDD-ET2" in the "EX-BIOC" portfolio
    When user clicks on biologics development portfolio
    And user clicks on scratch stale evidence project card
    And user clicks on choose process description run button
    And user selects "Feed 4" in selected run input
    And user enters "Higher titer at the same temperature" into selection rationale input
    And user clicks on select run button
    Then choose run dialog should not be visible
    When user clicks on generate process description evidence button
    Then process description document card should contain text "Evidence table, Locked process description · v1 · draft"
    And process description document card should not contain text "stale"
    When user clicks on choose process description run button
    And user selects "Feed 2" in selected run input
    And user enters "Lower feed holds the titer" into selection rationale input
    And user clicks on select run button
    Then choose run dialog should not be visible
    And process description document card should contain text "Evidence changed since written: another run selected: Feed 2 (was Feed 4)"
    When user clicks on regenerate process description evidence link
    Then process description document card should contain text "Evidence table, Locked process description · v2 · draft"
    And process description document card should contain text "Evidence table, Locked process description · v1 · superseded"
    And process description document card should not contain text "Evidence changed since written"
