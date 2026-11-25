@guide @demo
Feature: 1 · Any notation · HELM rendering and editor
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves
  the set-up and the checks off camera). Helm's own showcase sample: HELM of every shape — disulfide-rich
  and lasso peptides, siRNA duplexes, a tRNA cloverleaf, peptide–oligo and BLOB conjugates — drawn in the
  grid as their graphs, a monomer's structure under the pointer, and the HELM editor on a two-chain,
  disulfide-bridged peptide: its BILN (the editor opens on its HELM tab), then the whole molecule's structure,
  whose Explore pane holds Chem's panels for it; the clip ends inside the editor on the molecule's Properties.
  The editor's palette is switched to the peptides off camera, so it does not open on an empty Favorites tab.
  Nothing is saved: the scenario's end closes the editor the way Cancel does. Left out: the Properties tab (its
  extinction coefficient carries no unit).

  Scenario: Any HELM drawn in the grid, edited in one editor
    Peptides, duplexes, folded RNA and conjugates drawn from HELM; a double-click opens the editor.

    Given user is logged in
    And simple mode is off
    And the Helm package is initialized
    And user opens helm-showcase dataset
    Then "HELM" column should have units "helm"
    Given the browse panel is closed
    And the toolbox pane is hidden
    And the context panel is closed
    And user filters rows where "Name" is one of "Disulfide-rich (insulin-like), Lasso peptide (macrolactam + tail), siRNA with ss/as strand annotations, Cloverleaf — four-arm tRNA, Peptide – SMCC – RNA conjugate"
    And the grid shows only the columns "Name, Category, HELM"
    And the "Name" column of the grid is 230 px wide
    And the "HELM" column of the grid is 900 px wide
    And the rows of the grid are 150 px tall
    When user hovers over monomer 20 of the "cell 33 of HELM" area of grid
    And the video says "One renderer for peptides, duplexes, folded RNA and conjugates"
    And user hovers over monomer 6 of the "cell 14 of HELM" area of grid
    And user double-clicks on the "cell 14 of HELM" area of grid
    Then HELM editor should be visible
    Given user clicks on Peptides palette tab
    When user clicks on "BILN" tab in HELM editor
    And user clicks on "Molecular Structure" tab in HELM editor
    And the video says "The whole molecule, with every chemistry panel"
    And user expands "Explore" section in HELM editor
    And user expands the "Chemistry" pane inside the "Explore" section
    And user expands the "Properties" pane inside the "Chemistry" section
    Then no errors should have been logged
