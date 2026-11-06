@guide
Feature: Keep a molecule's panels in view while clicking through molecules
  A guide: the answer to "how do I see the 3D structure, the Gasteiger charges and the properties of each molecule I
  click, without reopening them in the Context Panel?". A panel a function computes for a molecule, dragged out of the
  Context Panel, becomes a viewer of its own: it recalculates for each molecule that becomes current, and docks into
  the table view like any other viewer, so it is saved with the layout. It holds no table of its own, and clicking
  anything other than a molecule leaves it as it is.

  In Crux mode the panels compute in the browser; in RDKit mode Gasteiger Partial Charges and 3D Structure are server-side
  Python, so this guide skips there. The rules behind it are in panels/floating-pane.feature and the platform's tests
  `xamgle | Context panel | Dragged out pane`.

  Scenario: Dock Properties, Gasteiger and 3D under the grid
    Given user is logged in
    And simple mode is off
    And the package autostarts have completed
    And the Chemical Engine is "Crux"
    And user opens smiles-50 dataset
    And the context panel is open
    # caption: Click a molecule
    When user clicks on the "cell 1 of canonical_smiles" area of grid
    # caption: Its panels are in the Context Panel
    Then the context panel should show the current cell
    # caption: Open the Chemistry group
    When user expands Chemistry accordion header in context panel
    # caption: Drag Properties out of the Context Panel: it becomes a viewer
    And user drags Properties accordion header in context panel out of the accordion
    # caption: It shows the properties of this molecule
    Then "Properties" floating panel should contain the text "C24H23NO3"
    # caption: Dock it under the grid
    When user docks info panel viewer in "Properties" floating panel to the bottom edge of the view
    # caption: Drag out Gasteiger Partial Charges
    And user drags "Gasteiger Partial Charges" accordion header in context panel out of the accordion
    # caption: Dock it next to Properties
    And user docks info panel viewer in "Gasteiger Partial Charges" floating panel to the right side of info panel viewer in "Properties" dock panel
    # caption: Open the Structure group
    And user expands Structure accordion header in context panel
    # caption: Drag out 3D Structure
    And user drags "3D Structure" accordion header in context panel out of the accordion
    # caption: Dock it next to Gasteiger Partial Charges
    And user docks info panel viewer in "3D Structure" floating panel to the right side of info panel viewer in "Gasteiger Partial Charges" dock panel
    # caption: Make the docked panels taller
    And user drags the top border of info panel viewer in "Properties" dock panel by 210 pixels up
    # caption: Share the width between them
    And user drags the right border of info panel viewer in "Properties" dock panel by 132 pixels left
    # caption: And with 3D Structure
    And user drags the right border of info panel viewer in "Gasteiger Partial Charges" dock panel by 200 pixels left
    # caption: Click another molecule
    And user clicks on the "cell 3 of canonical_smiles" area of grid
    # caption: Properties recalculates for it
    Then "Properties" dock panel should contain the text "C22H23NO3"
    # caption: So do the charges and the 3D structure
    And no loading indicator should be visible
    # caption: Click another one
    When user clicks on the "cell 5 of canonical_smiles" area of grid
    # caption: Its formula in Properties
    Then "Properties" dock panel should contain the text "C12H14F3NO"
    # caption: The panels stay with the molecule you click
    And no loading indicator should be visible
    And no errors should have been logged
