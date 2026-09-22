@journey
Feature: Peptides landing view
  The landing view offers examples in FASTA, separator and HELM notation.

  Not translated, and why: opening it from the app browser, as the manual case does — Peptides
  registers no app (its Browse > Apps "Peptides" group is other packages' browsePath), so the view
  is opened through the package's entry function, Peptides:Peptides.

  Background:
    Given user is logged in
    And the Peptides package is initialized
    And user opens the Peptides landing view
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The landing view offers three demos and folds the side panels
    Then the "Peptides" view should be current
    And there should be 3 visible button in Peptides landing view
    And "Simple demo" button should be visible
    And "Complex demo" button should be visible
    And "HELM demo" button should be visible
    And toolbox should be hidden
    And context panel should be hidden
    And help panel should be hidden
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario Outline: A demo button opens the dataset in its original notation
    When user clicks on "<demo> demo" button
    Then the "PeptidesView" view should be current
    And the table should have <rows> rows
    And "<column>" column should have semantic type "Macromolecule"
    And "<column>" column should have units "<notation>"
    And context panel should be visible
    When user closes the current view
    And user switches to the "Peptides" view
    Then no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | demo    | rows | column          | notation  |
      | Simple  | 647  | AlignedSequence | fasta     |
      | Complex | 540  | MSA             | separator |
      | HELM    | 334  | HELM            | helm      |
