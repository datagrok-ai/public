@guide @tutorial
Feature: How-to · Basics
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). The context panel on the right shows whatever was clicked last,
  in panes that open with a click: the table's name in the status bar brings the table (General opened, on demog), a
  column header the column (Stats and Plots opened), and a molecule cell the molecule (smiles.csv: Chemistry >
  Properties opened, computed in the browser). A plain cell or a row number is no object of its own: the panel stays where it
  was, so the clip does not click one. The table comes first: with a column current, a click on the table's name
  left the panel on the column (2026-11-22). demog is opened off camera with the context panel open;
  the table of molecules is opened on camera from Browse > Files > Demo > chem (no jump); nothing is made on the
  server. The molecule's other panes stay
  closed: Descriptors, Gasteiger, 3D Structure and Identifiers run server scripts or reach outside services.

  Scenario: The context panel
    The panel on the right shows what you click: a column, the table, a molecule.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And the toolbox pane is hidden
    And the context panel is open
    # caption: The panel on the right shows what you click: here, the table's name
    When user clicks on "demog" text in status bar
    Then the context panel should show "demog"
    # caption: Its panes open with a click
    When user clicks on "General" accordion header in context panel
    # caption: A column: click its header
    And user clicks on the "header AGE" area of grid
    Then the context panel should show "AGE"
    # caption: Stats: the column's numbers at a glance
    When user clicks on "Stats" accordion header in context panel
    # caption: Plots: how the values are spread
    And user clicks on "Plots" accordion header in context panel
    # caption: Now a table of molecules, opened from Browse
    When user clicks on browse tab
    And user expands Files tree node inside browse tree
    And user expands Files---Demo tree node inside browse tree
    And user expands Files---Demo---chem tree node inside browse tree
    And user double-clicks on Files---Demo---chem---smiles.csv tree node inside browse tree
    Then the "smiles" view should be current
    # caption: A molecule: click its cell
    When user clicks on the "cell 3 of canonical_smiles" area of grid
    Then the context panel should show the current cell
    # caption: Chemistry > Properties: computed for this molecule
    When user clicks on "Chemistry" accordion header in context panel
    And user clicks on "Properties" accordion header in context panel
    Then no errors should have been logged
