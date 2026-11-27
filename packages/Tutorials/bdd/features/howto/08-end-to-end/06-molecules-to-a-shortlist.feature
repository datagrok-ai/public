@guide @tutorial
Feature: How-to · End to end
  An end-to-end how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films
  it as a demo clip that outlines what each step presses). One task from start to finish: from 100 compounds (spgi,
  opened off camera with the toolbox open and the OpenChemLib sketcher) to a shortlist. Chem > Search > Substructure
  Search with the pyrrolidine ring (C1CCNC1) keeps 33; a scatter plot of NIBR logP against TPSA colored by
  Competition assay shows them; a histogram of Competition assay puts the most potent (lowest values) on the left, and
  a selection box over its first two bins selects 7 of the 33; Select > Extract Selected Rows makes the shortlist a
  table of its own. Chem > Search > Substructure Search asks first which molecule column to search (spgi has more
  than one), then opens the sketcher. Only browser-side Chem functions (the substructure search runs on RDKit in the page). Every move after the
  start is on camera. Nothing is made on the server.

  Scenario: From molecules to a shortlist
    Narrow a set of compounds by structure, look at their properties and activity, and pull the best into a table.

    Given user is logged in
    And simple mode is off
    And the "Chem" package is installed
    And the molecule sketcher is "OpenChemLib"
    And the package autostarts have completed
    And user opens spgi dataset
    And the toolbox pane is shown
    And the context panel is closed
    And the browse panel is closed
    When the video says "The task: a shortlist of potent compounds with a pyrrolidine ring"
    # caption: Keep the molecules with the ring: Chem > Search > Substructure Search
    And user picks "Chem > Search > Substructure Search..." from the top menu
    Then "Substructure search" dialog should be visible
    When user clicks on OK button in "Substructure search" dialog
    Then sketcher dialog should be visible
    When user clicks on molecule input of sketcher dialog
    And user types "C1CCNC1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then 33 rows should pass the filter
    When the video says "33 of 100 have the ring"
    # caption: Plot two properties: logP against polar surface area
    And user clicks on "scatter plot" icon in toolbox
    Then scatter plot viewer should be visible
    When user picks "NIBR logP" in the "x" column selector of scatter plot viewer
    And user picks "TPSA" in the "y" column selector of scatter plot viewer
    # caption: Color by activity. Lower Competition assay values are more potent
    And user picks "Competition assay" in the "color" column selector of scatter plot viewer
    # caption: A histogram of the activity: the most potent are on the left
    And user clicks on "histogram" icon in toolbox
    Then histogram viewer should be visible
    When user picks "Competition assay" in the "value" column selector of histogram viewer
    # caption: Shift+drag over the first bars to select the most potent
    And user drags a selection box from the "bin 1" area to the "bin 2" area of histogram viewer
    Then 7 rows should be selected
    And every selected row should pass the filter
    When the video says "The plot marks where they sit"
    And user hovers over scatter plot viewer
    # caption: Select > Extract Selected Rows: the shortlist as a table of its own
    And user picks "Select > Extract Selected Rows" from the top menu
    Then the table should have 7 rows
    When the video says "The shortlist, with structures and every property"
    And user hovers over the "cell 1 of Structure" area of grid
    Then no errors should have been logged
