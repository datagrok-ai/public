@guide @tutorial
Feature: How-to · Sequences
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). Bio's antibodies sample (493 heavy and light chains): Bio >
  Annotate > Apply Numbering Scheme numbers the heavy chains in IMGT, in the browser (the HMM package's ANARCI,
  the default engine), colors the frameworks and CDRs on every chain (the pointer on a residue names its region),
  and adds the chains aligned in IMGT positions, whose header shows the regions; the context menu of a numbered
  sequence extracts CDR3 as a column of its own, keeping its IMGT positions, and a drag on its edge widens it. The
  wheel with Shift held over the aligned column's header scrolls it three positions a notch, to CDR3. Off
  camera, before the first step: the grid's columns chosen and sized. Nothing is made on the server.

  Scenario: Number antibody chains
    Number antibody chains in IMGT: frameworks and CDRs marked on every chain, and CDR3 as a column of its own.

    Given user is logged in
    And simple mode is off
    And the "HMM" package is installed
    And user opens antibodies dataset
    And the Bio package is initialized
    And the toolbox pane is hidden
    And the context panel is closed
    And the grid shows only the columns "Antibody_ID, AntibodyHC"
    And the "AntibodyHC" column of the grid is 600 px wide
    When the video says "493 antibodies, their heavy chains here"
    # caption: Bio > Annotate > Apply Numbering Scheme
    And user picks "Bio > Annotate > Apply Numbering Scheme..." from the top menu
    Then "Apply Antibody Numbering" dialog should be visible
    And Scheme input in "Apply Antibody Numbering" dialog should have value "imgt"
    # caption: Pick a scheme: IMGT, Kabat, Chothia and more
    When user selects "imgt" in Scheme input in "Apply Antibody Numbering" dialog
    And user clicks on OK button in "Apply Antibody Numbering" dialog
    Then a new column "AntibodyHC (aligned)" should have been added
    And "AntibodyHC (aligned)" column should have tag ".numberingScheme" equal to "imgt"
    When the video says "Frameworks and CDRs, colored on every chain"
    And user hovers over monomer 30 of the "cell 3 of AntibodyHC" area of grid
    Then tooltip should contain text "CDR1"
    When user moves the pointer away from grid
    And the video says "The new column: the chains aligned in IMGT positions"
    # caption: Shift + wheel over its header scrolls the positions
    And user scrolls the mouse wheel down 43 times over the "header AntibodyHC (aligned)" area of grid holding Shift
    And user hovers over the region track under the "position 111A of AntibodyHC (aligned) header" area of grid
    Then tooltip should contain text "CDR3"
    # caption: Right-click a chain: extract CDR3 as a column
    When user picks "Annotations > Extract CDR3 as Column" from the context menu of the "cell 12 of AntibodyHC (aligned)" area of grid
    Then a new column "AntibodyHC (aligned)(CDR3)" should have been added
    # caption: Drag its edge to widen it
    When user drags the "column resizer AntibodyHC (aligned)(CDR3)" area of grid by 200 pixels to the right
    And the video says "CDR3 alone, with its IMGT positions: 105 to 117"
    And user hovers over the "position 111A of AntibodyHC (aligned)(CDR3) header" area of grid
    Then tooltip should contain text "Position: 111A"
    And no error or warning balloon should have been shown
    And no errors should have been logged
