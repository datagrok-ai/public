@guide @help:access/files
Feature: Open a file from your computer
  A guide: the answer to "how do I open a CSV file from my computer in Datagrok?", written as a
  scenario so that it is also a test of that answer. `grok-bdd guide features/guides/import-local-csv.feature`
  renders it into a video with every step's element lit and the pointer moving to it.

  Scenario: Open a local CSV file as a table
    Given user is logged in
    And simple mode is off
    And the browse panel is open
    When user uploads "fixtures/browse-import.csv" through "Open local file" icon inside browse toolbar
    Then the "browse-import" view should be current
    And the table should have 5 rows
    And the table should have 3 columns
