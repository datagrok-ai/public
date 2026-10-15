@sketcher-controls
Feature: A query drawn in the substructure filter on SPGI: its rows, the sketcher reopened, and Ketcher
  The product owner's report of 2026-10-07 (crux-sketch spike query-roundtrip, H47): "Open mol1k. Open filters panel. In
  structure filter draw toluene (with crux sketcher). Mark the bond of benzene to methyl as aromatic. Click ok in dialog.
  Problems: Not all things passing the filter are correct ... if I click on query molecule again to open sketcher back,
  crux sketcher does not recognize the aromatic bond that was marked and draws just normal toluene. From there, if I
  switch to Ketcher sketcher, that one recognizes the aromatic bond on toluene." The query's rows are claimed against
  RDKit's matches of its SMARTS, `[#6]:[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1` (the methyl bond aromatic); what Crux
  holds through its "smarts" reading (its own SMARTS), what Ketcher shows through its own molfile, each matched against
  the column the same way. The other query features, on mol1K, SPGI and smiles.csv, with Crux and with Ketcher, are the
  package test category "query round trip"; this one checks the report's flow on SPGI (spgi-100), its molblock column.
  The feature drives Crux's own controls (@sketcher-controls): a run that pins another sketcher skips it.

  Background:
    Given user is logged in
    And the molecule sketcher is "Crux"
    And the package autostarts have completed

  Scenario: Toluene drawn in Crux, its methyl bond then marked aromatic, filters by the aromatic bond, and the sketcher reopened, Ketcher and Crux again hold that bond
    Given user opens spgi-100 dataset
    When user clicks on filter icon in toolbar
    And user clicks on "Sketch" text in "Structure" filter card
    And user clicks on crux benzene tool
    And user clicks on crux canvas
    And user clicks on crux single bond tool
    And user clicks on the "atom 0" area of crux sketcher widget
    Then the "smiles" reading of crux sketcher widget should be the molecule "Cc1ccccc1"
    And the filter should pass exactly the molecules of "Structure" column containing "[#6]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1"
    When user clicks on crux aromatic bond tool
    And user clicks on the "bond 6" area of crux sketcher widget
    And user clicks on OK button in sketcher dialog
    Then the filter should pass exactly the molecules of "Structure" column containing "[#6]:[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1"
    When user clicks on sketcher thumbnail in "Structure" filter card
    Then the "smarts" reading of crux sketcher widget should find exactly the molecules of "Structure" column containing "[#6]:[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1"
    When user clicks on "Options" icon in sketcher dialog
    And user picks "Ketcher" from the open menu
    Then the query Ketcher shows in sketcher dialog should find exactly the molecules of "Structure" column containing "[#6]:[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1"
    When user clicks on "Options" icon in sketcher dialog
    And user picks "Crux" from the open menu
    Then the "smarts" reading of crux sketcher widget should find exactly the molecules of "Structure" column containing "[#6]:[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1"
