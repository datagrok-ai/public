@portal
Feature: Documents over files
  A document is a file with a kind, title, version and status, linked to the gate criterion it
  answers. The program page has a card per document criterion of the current and next gate:
  what exists, what is missing, and the evidence behind it. Approving a document meets its
  criterion.

  Background:
    Given user opens the Project Portal app

  Scenario: A file of the project folder is linked to the process description and approved
    Given a scratch project "BDD-DOCS" at stage "D4" in the "EX-BIOC" portfolio
    And the workspace of the scratch project "BDD-DOCS" in the "EX-BIOC" portfolio is built
    And a file "process-description-v1.md" in the folder of the scratch project "BDD-DOCS" in the "EX-BIOC" portfolio
    When user clicks on biologics development portfolio
    And user clicks on scratch docs project card
    Then process description document card should contain text "no document yet"
    When user clicks on link process description document button
    Then link document dialog should be visible
    When user selects "process-description-v1.md" in document file input
    And user enters "Process description, BDD" into document title input
    And user enters "1.0" into document version input
    And user clicks on link document button
    Then link document dialog should not be visible
    And process description document card should contain text "Process description, BDD · v1.0 · draft"
    And process description document card should contain text "in progress"
    When user clicks on approve process description document link
    Then process description document card should contain text "Process description, BDD · v1.0 · approved"
    And process description document card should contain text "met"
