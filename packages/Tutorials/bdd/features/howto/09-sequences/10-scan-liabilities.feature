@guide @tutorial
Feature: How-to · Sequences
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). Bio's antibodies sample: Bio > Annotate > Scan Liabilities scans
  the heavy chains for sequence liabilities (deamidation, isomerization, oxidation and the like), in the browser, and
  marks every hit on the chains; the pointer on a marked residue names the liability. The heavy chains are numbered
  in IMGT off camera first, as the antibody numbering how-to shows, so the scan runs on the CDR3 the context menu
  extracts on camera here. Off camera, before the first step: the numbering, the grid's columns chosen and
  sized, the rows made taller. Nothing is made on the server.

  Scenario: Scan sequences for liabilities
    Scan Liabilities marks deamidation, isomerization, oxidation and other motifs on every sequence.

    Given user is logged in
    And simple mode is off
    And the "HMM" package is installed
    And user opens antibodies dataset
    And the Bio package is initialized
    And the toolbox pane is hidden
    And the context panel is closed
    And user picks "Bio > Annotate > Apply Numbering Scheme..." from the top menu
    And user clicks on OK button in "Apply Antibody Numbering" dialog
    Then a new column "AntibodyHC (aligned)" should have been added
    Given the notification balloons are closed
    And the grid shows only the columns "Antibody_ID, AntibodyHC (aligned)"
    And the header of the "AntibodyHC (aligned)" column starts at position "95"
    And the "AntibodyHC (aligned)" column of the grid stays 560 px wide
    And the rows of the grid are 30 px tall
    When the video says "Antibody heavy chains, numbered in IMGT"
    # caption: Right-click a chain: extract CDR3 as a column
    And user picks "Annotations > Extract CDR3 as Column" from the context menu of the "cell 12 of AntibodyHC (aligned)" area of grid
    Then a new column "AntibodyHC (aligned)(CDR3)" should have been added
    # caption: Drag its edge to widen it
    When user drags the "column resizer AntibodyHC (aligned)(CDR3)" area of grid by 300 pixels to the right
    # caption: Bio > Annotate > Scan Liabilities
    When user picks "Bio > Annotate > Scan Liabilities..." from the top menu
    Then "Scan Sequence Liabilities" dialog should be visible
    # caption: Pick the column to scan: the CDR3s
    When user selects "AntibodyHC (aligned)(CDR3)" in Sequence input in "Scan Sequence Liabilities" dialog
    And user clicks on OK button in "Scan Sequence Liabilities" dialog
    Then a new column "~AntibodyHC (aligned)(CDR3)_annotations" should have been added
    When the video says "Each liability found is underlined on its CDR3"
    And user hovers over monomer 3 of the "cell 12 of AntibodyHC (aligned)(CDR3)" area of grid
    Then tooltip should contain text "Isomerization (DG)"
    And no error or warning balloon should have been shown
    And no errors should have been logged
