@guide @demo
Feature: 1 · Any notation · Sequence formats
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves the
  set-up and the checks off camera). Bio's sample files in FASTA, separator, HELM and BILN open from Browse as
  Macromolecule columns, drawn by one renderer: monomers in their colors, the monomer under the pointer in a
  tooltip, a BILN line with its disulfide drawn as a bracket; the clip ends with the pointer off the grid, the
  brackets in view. One note per notation says what differs. Set up off camera: the Browse panel narrowed so the
  sequences get the room (docked afresh at its default width at feature end), the toolbox and the context panel
  closed, the sequence columns widened and their letters enlarged, HELM's activities shown to two decimals,
  BILN_complex shown without its notes, its monomer names up to eight letters (Thr-ol, D-2Nal) and in rows tall
  enough for its brackets. Left out: aligned (MSA) files, which the alignment clips show.

  Scenario: Every notation opens as a Macromolecule column
    FASTA, separator, HELM and BILN open as one column type, drawn and explained alike.

    Given user is logged in
    And simple mode is off
    And the Bio package is initialized
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---Bio tree node inside browse tree is expanded
    And Files---App-Data---Bio---samples tree node inside browse tree is expanded
    And the browse panel is 50 px narrower
    When user double-clicks on Files---App-Data---Bio---samples---FASTA.fasta tree node inside browse tree
    Then the "FASTA" view should be current
    And "sequence" column should have semantic type "Macromolecule"
    And "sequence" column should have units "fasta"
    Given the toolbox pane is hidden
    And the context panel is closed
    And the "sequence" column of the grid is 480 px wide
    And the monomers of the "sequence" column are drawn at 14 px
    When the video says "FASTA: one letter per residue"
    And user hovers over monomer 12 of the "cell 2 of sequence" area of grid
    And user double-clicks on Files---App-Data---Bio---samples---SEPARATOR_PT.csv tree node inside browse tree
    Then the "SEPARATOR_PT" view should be current
    And "sequence" column should have semantic type "Macromolecule"
    And "sequence" column should have units "separator"
    Given the toolbox pane is hidden
    And the context panel is closed
    And the "sequence" column of the grid is 700 px wide
    And the monomers of the "sequence" column are drawn at 14 px
    When the video says "Separator: monomers between dashes"
    And user double-clicks on Files---App-Data---Bio---samples---HELM.csv tree node inside browse tree
    Then the "HELM" view should be current
    And "HELM" column should have semantic type "Macromolecule"
    And "HELM" column should have units "helm"
    Given the toolbox pane is hidden
    And the context panel is closed
    And the "HELM" column of the grid is 620 px wide
    And user sets the format of "Activity" column in the grid to "0.00"
    When the video says "HELM: non-natural monomers, colored like their natural analogs"
    And user hovers over monomer 9 of the "cell 3 of HELM" area of grid
    And user double-clicks on Files---App-Data---Bio---samples---BILN_complex.csv tree node inside browse tree
    Then the "BILN_complex" view should be current
    And "biln" column should have semantic type "Macromolecule"
    And "biln" column should have units "biln"
    Given the toolbox pane is hidden
    And the context panel is closed
    And the grid shows only the columns "name, category, biln"
    And the "biln" column of the grid is 620 px wide
    And the monomers of the "biln" column are drawn at 14 px
    And the monomer names of the "biln" column are cut at 8 characters
    And the rows of the grid are 56 px tall
    When the video says "BILN: rings and cross-links written in the line"
    And user hovers over monomer 2 of the "cell 2 of biln" area of grid
    And user moves the pointer away from grid
    Then no errors should have been logged
