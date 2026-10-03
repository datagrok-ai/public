@portal
Feature: Data package manifest
  The program page's Data package lists what the program's gate case rests on, one line each:
  the criteria up to the next gate with their state, documents, evidence tables, saved runs, the
  program's datasets and snapshots, each with its status and where it comes from. It opens as a
  table and downloads as CSV.

  Background:
    Given user opens the Project Portal app

  Scenario: A project's data package lists its criteria and its saved run
    Given a scratch project "BDD-PKG" at stage "D5" in the "EX-BIOC" portfolio
    And a run of the example model saved to the scratch project "BDD-PKG" in the "EX-BIOC" portfolio
    When user clicks on biologics development portfolio
    And user clicks on scratch package project card
    And user clicks on data package button
    Then data package dialog should be visible
    And data package summary should contain text "1 run"
    And data package lines should contain text "Fed-batch titer (example) · BDD-PKG"
    And data package lines should contain text "Production bioreactor"
    And data package lines should contain text "Process description"
    And download data package csv button should be visible
    When user clicks on open data package as table button
    Then data package dialog should not be visible
    And the open table should have 6 rows and a name starting with "BDD-PKG data package"
