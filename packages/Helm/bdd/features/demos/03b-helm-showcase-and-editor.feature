@guide @demo
Feature: Bioinformatics · HELM rendering and editor
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves
  the set-up and the checks off camera). Helm's own showcase sample: HELM of every shape — disulfide-rich
  and lasso peptides, siRNA duplexes, a tRNA cloverleaf, peptide–oligo and BLOB conjugates — drawn in the
  grid as their graphs, a monomer's structure under the pointer, and the HELM editor on a two-chain,
  disulfide-bridged peptide with its notation, computed properties and the whole molecule's structure. Cancel leaves the cell as it was.

  Scenario: Any HELM drawn in the grid, edited in one editor
    Peptides, cyclic and branched, oligonucleotide duplexes, folded RNA and conjugates are drawn from their HELM; a double-click opens the editor.

    Given user is logged in
    And simple mode is off
    And the Helm package is initialized
    And user opens helm-showcase dataset
    Then "HELM" column should have units "helm"
    Given user filters rows where "Name" is one of "Disulfide-rich (insulin-like), Lasso peptide (macrolactam + tail), siRNA with ss/as strand annotations, Cloverleaf — four-arm tRNA, Peptide – SMCC – RNA conjugate"
    And the grid shows only the columns "Name, Category, HELM"
    And the "Name" column of the grid is 230 px wide
    And the "HELM" column of the grid is 900 px wide
    And the rows of the grid are 150 px tall
    When user hovers over monomer 20 of the "cell 33 of HELM" area of grid
    And the video says "One renderer for peptides, duplexes, folded RNA and conjugates"
    And user hovers over monomer 6 of the "cell 14 of HELM" area of grid
    And user double-clicks on the "cell 14 of HELM" area of grid
    Then HELM editor should be visible
    When user clicks on HELM tab
    Then notation pane should contain text "PEPTIDE1"
    When user clicks on Properties tab
    Then formula field should be visible
    When user clicks on "Molecular Structure" tab in HELM editor
    And user clicks on CANCEL button in HELM editor
    Then HELM editor should be absent
    And no errors should have been logged
