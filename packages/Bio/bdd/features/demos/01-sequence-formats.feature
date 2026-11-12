@guide @demo
Feature: Bioinformatics · Sequence formats
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag
  leaves the set-up and the checks off camera). Bio's sample files in FASTA, separator, HELM, BILN
  and aligned notation open from Browse as Macromolecule columns, drawn by one renderer: monomers in
  their colors, the monomer under the pointer in a tooltip, aligned sequences position by position.

  Scenario: Every notation opens as a Macromolecule column
    FASTA, separator, HELM, BILN and aligned files open as one column type: monomers drawn in color, details under the pointer and the same tools for every notation.

    Given user is logged in
    And simple mode is off
    And the Bio package is initialized
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---Bio tree node inside browse tree is expanded
    And Files---App-Data---Bio---samples tree node inside browse tree is expanded
    When user double-clicks on Files---App-Data---Bio---samples---FASTA.fasta tree node inside browse tree
    Then the "FASTA" view should be current
    And "sequence" column should have semantic type "Macromolecule"
    And "sequence" column should have units "fasta"
    Given the toolbox pane is hidden
    And the context panel is closed
    And the "sequence" column of the grid is 420 px wide
    When the video says "FASTA: one letter per residue"
    And user hovers over monomer 12 of the "cell 2 of sequence" area of grid
    And user double-clicks on Files---App-Data---Bio---samples---SEPARATOR_PT.csv tree node inside browse tree
    Then the "SEPARATOR_PT" view should be current
    And "sequence" column should have semantic type "Macromolecule"
    And "sequence" column should have units "separator"
    Given the toolbox pane is hidden
    And the context panel is closed
    And the "sequence" column of the grid is 620 px wide
    When the video says "Separator: monomers between dashes"
    And user double-clicks on Files---App-Data---Bio---samples---HELM.csv tree node inside browse tree
    Then the "HELM" view should be current
    And "HELM" column should have semantic type "Macromolecule"
    And "HELM" column should have units "helm"
    Given the toolbox pane is hidden
    And the context panel is closed
    And the "HELM" column of the grid is 620 px wide
    When the video says "HELM: non-natural monomers, colored like their natural analogs"
    And user hovers over monomer 9 of the "cell 3 of HELM" area of grid
    And user double-clicks on Files---App-Data---Bio---samples---BILN.csv tree node inside browse tree
    Then the "BILN" view should be current
    And "biln" column should have semantic type "Macromolecule"
    And "biln" column should have units "biln"
    Given the toolbox pane is hidden
    And the context panel is closed
    And the "biln" column of the grid is 620 px wide
    When the video says "BILN: chains and ring bonds in one line"
    And user double-clicks on Files---App-Data---Bio---samples---MSA.csv tree node inside browse tree
    Then the "MSA" view should be current
    And "MSA" column should have semantic type "Macromolecule"
    Given the toolbox pane is hidden
    And the context panel is closed
    And the "MSA" column of the grid is 620 px wide
    When the video says "Aligned: monomers line up position by position"
    And user hovers over monomer 10 of the "cell 4 of MSA" area of grid
    Then no errors should have been logged
