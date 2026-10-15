@guide @sketcher-controls
Feature: Rings and templates in Crux
  A guide: how do I draw ring systems in Crux, Chem's molecule sketcher? The ring bar below the canvas places a ring on
  the empty canvas, fuses it onto a bond it is clicked on, and puts a spiro ring on an atom. The Structure Library (SL,
  the ring bar's last button) holds Ketcher's templates: a search shows the groups that match, opened, and a template
  picked there is placed with a click. Each result is the molecule the sketcher holds, compared by RDKit. The feature
  drives Crux's own controls (@sketcher-controls): a run that pins another sketcher skips it.

  Scenario: Place a ring, fuse another, add a spiro ring and a template from the Structure Library
    Given user is logged in
    And simple mode is off
    And the molecule sketcher is "Crux"
    And the package autostarts have completed
    And the Crux sketcher is open
    # caption: Pick benzene from the ring bar
    When user clicks on Crux benzene tool
    # caption: Click the canvas to place it
    And user clicks on Crux canvas
    # caption: Pick cyclohexane
    And user clicks on Crux cyclohexane tool
    # caption: Click a bond of the benzene: the new ring is fused onto it
    And user clicks on the "bond 0" area of Crux sketcher widget
    # caption: Tetralin: a benzene fused with a cyclohexane
    Then the sketcher in sketcher dialog should hold the molecule "c1ccc2c(c1)CCCC2"
    # caption: Pick cyclopropane
    When user clicks on Crux cyclopropane tool
    # caption: Click a CH2 of the new ring: a spiro ring grows there
    And user clicks on the "atom 7" area of Crux sketcher widget
    # caption: A spiro compound: the two rings share one carbon
    Then the sketcher in sketcher dialog should hold the molecule "c1ccc2c(c1)CCC1(CC1)C2"
    # caption: Open the Structure Library
    When user clicks on Crux structure library button
    # caption: Search for azulene: its group opens with the match
    And user types "azulene" into Crux structure library search
    # caption: Pick azulene
    And user clicks on "Azulene" button inside Crux structure library
    # caption: Click an empty spot of the canvas to place it
    And user clicks on Crux canvas 80% across and 25% down
    # caption: Azulene sits beside the spiro compound
    Then the sketcher in sketcher dialog should hold the molecule "c1ccc2c(c1)CCC1(CC1)C2.c1ccc2cccc-2cc1"
