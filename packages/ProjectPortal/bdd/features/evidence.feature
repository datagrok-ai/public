@portal
Feature: Decisions with evidence
  A gate decision rests on evidence: the gate's readiness, the runs selected for its criteria, and
  documents or pages the approver attaches. Each piece keeps its name as it was when attached, so
  the program's history reads the same after a run is renamed or deleted.

  Background:
    Given user opens the Project Portal app

  Scenario: An approver advances a project with its readiness and a link as evidence
    Given a scratch project "BDD-EVIDENCE" at stage "TV" in the "EX-SMD" portfolio
    When user clicks on small molecule discovery portfolio
    And user clicks on evidence project card
    And user clicks on open program button
    And user clicks on assay criterion met button
    Then decision evidence should contain text "Readiness of Hit Finding: 1 criterion"
    When user enters "https://example.org/reviews/hf-gate" into evidence link input
    And user clicks on attach evidence button
    Then decision evidence should contain text "Link: example.org/…/hf-gate"
    When user enters "Assay validated, Z' 0.7" into rationale input
    And user clicks on advance button
    Then latest decision evidence should contain text "Link: example.org/…/hf-gate"
    And latest decision evidence should contain text "Criterion: Primary assay validated"

  Scenario: A decision takes a snapshot of the project now, which opens read-only from the history
    Given a scratch project "BDD-EVIDENCE" at stage "TV" in the "EX-SMD" portfolio
    And the workspace of the scratch project "BDD-EVIDENCE" in the "EX-SMD" portfolio is built
    When user clicks on small molecule discovery portfolio
    And user clicks on evidence project card
    And user clicks on open program button
    And user clicks on assay criterion met button
    And user checks evidence snapshot checkbox
    And user enters "Assay validated" into rationale input
    And user clicks on advance button
    Then latest decision evidence should contain text "Snapshot: BDD-EVIDENCE"
    And program snapshots should contain text "Target Validation"
    When user clicks on latest decision snapshot link
    Then snapshot dialog should be visible
    And snapshot stage should contain text "Target Validation"
    And snapshot criteria should contain text "Primary assay validated"
    And snapshot table should contain text "the file holds what was frozen"
