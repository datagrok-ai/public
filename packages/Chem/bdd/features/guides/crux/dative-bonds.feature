@guide @sketcher-controls
Feature: Dative bonds in Crux
  A guide: how do I draw a dative (coordinate) bond in Crux, Chem's molecule sketcher? The dative bond tool draws one
  from its donor, where the drag starts, to its acceptor, drawn as an arrow onto it: here pyridine's nitrogen donates
  to a platinum atom. Crux writes the bond as `->` in the SMILES the sketcher holds, read as written. The feature
  drives Crux's own controls (@sketcher-controls): a run that pins another sketcher skips it.

  Scenario: Bind pyridine's nitrogen to a platinum atom with a dative bond
    Given user is logged in
    And simple mode is off
    And the molecule sketcher is "Crux"
    And the package autostarts have completed
    And the Crux sketcher is open on "c1ccncc1.[Pt]"
    # caption: Pick the dative bond tool
    When user clicks on Crux dative bond tool
    # caption: Drag from the nitrogen, the donor, to the platinum, the acceptor
    And user drags the "atom 3" area of Crux sketcher widget to the "atom 6" area
    # caption: Clean Up tidies the drawing
    And user clicks on Crux clean up button
    # caption: An N→Pt dative bond: pyridine bound to platinum
    Then the "smiles" reading of Crux sketcher widget should be "[Pt]<-[n]1ccccc1"
