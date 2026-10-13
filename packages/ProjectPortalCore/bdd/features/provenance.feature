@portal
Feature: Provenance
  From any record of the portal (a decision, a run, a document, a snapshot, a program) the
  provenance panel walks what it rests on, down to the data, and what rests on it, up to the
  program. Each record in it opens.

  Background:
    Given user opens the Project Portal app

  Scenario: The provenance of a decision: its readiness and snapshot, the data, and the project it moved
    Given a scratch project "BDD-PROV" at stage "TV" in the "EX-SMD" portfolio
    When user clicks on small molecule discovery portfolio
    And user clicks on provenance project card
    And user clicks on open program button
    And user clicks on assay criterion met button
    And user checks evidence snapshot checkbox
    And user enters "Assay validated" into rationale input
    And user clicks on advance button
    And user clicks on latest decision provenance link
    Then provenance dialog should be visible
    And provenance start should contain text "Decision: advance"
    And provenance rests on should contain text "Criterion: Primary assay validated"
    And provenance rests on should contain text "Snapshot: BDD-PROV"
    And provenance rests on should contain text "Data: SAR table of BDD-PROV"
    And provenance used by should contain text "Program: BDD-PROV"
    When user clicks on first used-by record
    Then provenance dialog should not be visible
    And program title should contain text "BDD-PROV"
