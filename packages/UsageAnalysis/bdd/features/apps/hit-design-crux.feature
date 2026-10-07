@crux @sketcher-controls
Feature: A row added in Hit Design, drawn in Crux
  Hit Design opens the molecule sketcher on the empty Molecule cell of a row the user adds ("Add new row" on the design
  view's ribbon), and OK writes the drawn molecule into the cell, which registers its V-iD. Here the sketcher is Crux,
  Chem's own, pinned as the session's sketcher: it opens empty and ready, and the molecule is drawn on its own controls.

  The campaign is a fixed fixture made once per stand and reused, since a campaign cannot be taken back from the
  database (bindings/hit-design.ts); its files are put back at feature end as they were. Needs Chem on the stand. The
  feature drives Crux's own controls (@sketcher-controls): a run that pins another sketcher skips it.

  Background:
    Given user is logged in
    And the "Chem" package is installed
    And the molecule sketcher is "Crux"
    And the package autostarts have completed
    And the Hit Design campaign "BDD Crux" is on the server, and is as it was again when the feature ends

  # HOST-052
  Scenario: A new row opens Crux empty and ready, and a molecule drawn there is written to the row with its V-iD
    Given user opens the "Hit Design" app
    When user clicks on "BDD Crux" link
    Then the "BDD Crux" view should be current
    And the table should have 1 row
    When user clicks on "Add new row" icon in toolbar
    Then sketcher dialog should be visible
    And the "ready" reading of crux sketcher widget should be "true"
    And the "atoms" reading of crux sketcher widget should be 0
    When user clicks on crux benzene tool
    And user clicks on crux canvas
    And user clicks on crux nitrogen tool
    And user clicks on the "atom 0" area of crux sketcher widget
    Then the "smiles" reading of crux sketcher widget should be the molecule "c1ccncc1"
    When user clicks on OK button in sketcher dialog
    Then sketcher dialog should be absent
    And the table should have 2 rows
    And the molecule in row 2 of "Molecule" column should be "c1ccncc1"
    And the Hit Design campaign "BDD Crux" should be saved with 2 rows, row 2 the molecule "c1ccncc1" with its V-iD
