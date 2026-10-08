@guide @sketcher-controls
Feature: Atropisomers in Crux
  A guide: how do I set the configuration of a biaryl's axis (an atropisomer) in Crux, Chem's molecule sketcher? A
  wedge at an end of the single bond between two rings that cannot turn makes that bond an axis, P or M, labelled
  beside it (R, S, E and Z labels shown, as Settings shows them). The axis's own menu offers P and M, the current one
  checked, and choosing the other redraws the wedge that gives it. The biaryl is the product owner's (2026-10-07), as a
  molfile. What the axis is, is read where its menu shows it checked. The feature drives Crux's own controls
  (@sketcher-controls): a run that pins another sketcher skips it.

  Scenario: Turn a biaryl's axis from P to M and back from its menu
    Given user is logged in
    And simple mode is off
    And the molecule sketcher is "Crux"
    And the package autostarts have completed
    And the Crux sketcher is open on this molfile, showing R, S, E and Z labels:
      """
      
           RDKit          2D
      
       16 17  0  0  0  0  0  0  0  0999 V2000
          0.6495    2.6250    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
         -0.6495    1.8750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
         -1.9486    2.6250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
         -3.2476    1.8750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
         -3.2476    0.3750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
         -1.9486   -0.3750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
         -1.9486   -1.8750    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
         -0.6495    0.3750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
          0.6495   -0.3750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
          1.9486    0.3750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
          1.9486    1.8750    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
          3.2476   -0.3750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
          3.2476   -1.8750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
          1.9486   -2.6250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
          0.6495   -1.8750    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
         -0.6495   -2.6250    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
        1  2  1  0
        2  3  1  0
        3  4  2  0
        4  5  1  0
        5  6  2  0
        6  7  1  0
        6  8  1  0
        8  9  1  0
        9 10  2  0
       10 11  1  0
       10 12  1  0
       12 13  2  0
       13 14  1  0
       14 15  2  0
       15 16  1  0
        8  2  2  0
        9 15  1  6
      M  END
      """
    # caption: Right-click the axis, the bond between the two rings
    When user opens the Crux context menu on the "bond 7" area
    # caption: The axis is P: its menu has P checked
    Then Crux atropisomer P item should be checked
    # caption: Choose M
    When user clicks on Crux atropisomer M item
    # caption: Right-click the axis again
    And user opens the Crux context menu on the "bond 7" area
    # caption: Now M is checked, and the label beside the axis reads (M)
    Then Crux atropisomer M item should be checked
    # caption: Choose P to turn it back
    When user clicks on Crux atropisomer P item
    # caption: Right-click the axis once more
    And user opens the Crux context menu on the "bond 7" area
    # caption: P again
    Then Crux atropisomer P item should be checked
