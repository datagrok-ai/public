@journey @realizes:bio.rendering @realizes:bio.rendering.biln @realizes:bio.rendering.separator @realizes:bio.rendering.monomer @realizes:GROK-12164
Feature: Sequence cell renderers
  A Macromolecule column is painted by the renderer its units tag selects, and the grid reports
  that renderer as the column's cell type: a HELM column by the HELM renderer, a separator column
  by the sequence renderer, each with its monomers in their own colors; a column a transform
  derives takes the renderer of its own notation, not its source's (GROK-12164); the per-position
  columns of Split to Monomers take the monomer renderer.

  Background:
    Given user is logged in
    And the Bio package is initialized

  Scenario: A HELM column renders with the HELM renderer in monomer colors
    Given user opens filter_HELM dataset
    Then "HELM string" column should have semantic type "Macromolecule"
    And "HELM string" column should have units "helm"
    And the "cell type of HELM string" reading of grid should be "helm"
    And the "cell 1 of HELM string" area of grid should be painted in at least 3 colors
    And no error or warning balloon should have been shown

  Scenario: A separator column renders with the sequence renderer in monomer colors
    Given user opens filter_MSA dataset
    Then "MSA" column should have semantic type "Macromolecule"
    And "MSA" column should have units "separator"
    And "MSA" column should have tag "separator" equal to "/"
    And the "cell type of MSA" reading of grid should be "sequence"
    And the "cell 1 of MSA" area of grid should be painted in at least 3 colors
    And no error or warning balloon should have been shown

  Scenario: Converting HELM to separator gives the new column the separator renderer
    Given user opens filter_HELM dataset
    When user picks "Bio > Transform > Convert Sequence Notation..." from the top menu
    And user selects "separator" in "Convert to" input in "Convert Sequence Notation" dialog
    And user clicks on OK button in "Convert Sequence Notation" dialog
    Then the top menu command should have completed
    And 1 new column should have been added
    And a new column matching "^separator\(HELM string\)" should have been added
    And "separator(HELM string)" column should have units "separator"
    And the "cell type of separator(HELM string)" reading of grid should be "sequence"
    And "HELM string" column should have units "helm"
    And the "cell type of HELM string" reading of grid should be "helm"
    And no error or warning balloon should have been shown

  Scenario: Split to Monomers columns render with the monomer renderer
    Given user opens filter_MSA dataset
    When user picks "Bio > Transform > Split to Monomers..." from the top menu
    And user clicks on OK button in "Split to Monomers" dialog
    Then the top menu command should have completed
    And 17 new columns should have been added
    And "1" column should have semantic type "Monomer"
    And the "cell type of 1" reading of grid should be "Monomer"
    And the "cell 1 of 1" area of grid should be painted
    And no error or warning balloon should have been shown
