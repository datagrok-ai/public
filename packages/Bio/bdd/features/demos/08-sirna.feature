@guide @demo
Feature: Bioinformatics · siRNA
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves the
  set-up and the checks off camera). Bio's siRNA demo (Bioinformatics | siRNA: demo-files/SIRNA.csv with its
  layout, opened off camera): HELM duplexes drawn by SequenceTranslator's oligo renderer, sense over antisense,
  each sugar modification a color band (2'-OMe, 2'-F), each linkage an arc (phosphorothioate and others), the
  GalNAc conjugate a pill. The pointer on a nucleotide shows its sugar, base and 3' linkage with their structures,
  on a linkage the linkage's structure, on the conjugate the conjugate's. The context panel then lists the
  modifications of the duplex with their colors, and draws its sense strand atom by atom (Bio's Molecular
  Structure, which the demo table also keeps as molfile columns), opened full screen under the magnifier. The
  hover points come from the oligo renderer (`getMonomerPoint` / `getLinkagePoint`, SequenceTranslator).

  Scenario: siRNA duplexes, every modification in view
    Modified sugars, linkages and conjugates of siRNA duplexes, explained under the pointer and drawn atom by atom.

    Given user is logged in
    And simple mode is off
    And the Bio package is initialized
    And user calls "Bio:demoBioSiRNA" function
    And the notification balloons are closed
    And the toolbox pane is hidden
    And the context panel is closed
    And the "oligo_helm" column of the grid is 640 px wide
    When user hovers over monomer 1 of the duplex in the "cell 1 of oligo_helm" area of grid
    And user hovers over monomer 2 of the duplex in the "cell 1 of oligo_helm" area of grid
    And user hovers over the linkage after monomer 2 of the duplex in the "cell 1 of oligo_helm" area of grid
    And user hovers over monomer 20 of the duplex in the "cell 1 of oligo_helm" area of grid
    And the video says "Every modification counted, each strand drawn atom by atom"
    Given the context panel is open
    Then the context panel should show the current cell
    When user expands "Oligo-Nucleotide" section
    And user expands "Oligo Structures" section
    And user expands the "Sense" pane inside the "Oligo Structures" section
    Then the "Sense" section should show a drawn molecule
    When user clicks on the molecule in the "Sense" section
    And user hovers over the molecule in the "Molecule" dialog at 30% across and 45% down
    And user hovers over the molecule in the "Molecule" dialog at 73% across and 38% down
    And user hovers over the molecule in the "Molecule" dialog at 69% across and 59% down
    Then no error or warning balloon should have been shown
