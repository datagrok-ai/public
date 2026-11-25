@guide @demo
Feature: Chemistry · Reactions in tables
  Demo clips of reaction cells (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves the set-up and the
  checks off camera). Crux draws a reaction cell scaled to the cell, and a ChemDraw page cut to its reaction (its key
  box and the rest of the page left out), legible as the cell grows. A double-click on a reaction cell opens the Crux
  reaction editor, in Crux mode and in RDKit mode alike; it goes full screen, follows the dialog when it is resized,
  and OK writes the edit into the cell. A ChemDraw page opened in the editor keeps its labels (OTf, TBS, Boc) through
  Crux's CDXML copy and paste (Ctrl+A selects the drawing only, not the rows of the table behind the dialog), and the
  Reaction pane of a page says how many molecules on no step it leaves out. The tables are Chem's own, opened from
  the package's files in the project (nothing written to the stand): files/tests/reactions-crux.csv (the reaction
  corpus), files/reactions-cdxml.csv (eight ChemDraw pages) and the showcase table
  fixtures/BDD-reactions-showcase.csv. What the cells draw and what the editor writes is Chem's package tests
  ('Reactions: CDXML', 'Reactions: editor').

  Scenario: Reaction cells, from SMILES and from ChemDraw
    Reactions and routes drawn to fit their cells; ChemDraw pages cut to their reaction, legible as the cells grow.

    Given user is logged in
    And simple mode is off
    And the Chemical Engine is "Crux"
    And user opens the "../files/tests/reactions-crux.csv" file of the project
    And the toolbox pane is hidden
    And the context panel is closed
    Given user filters rows where "kind" is one of "route: enumeration, route: linear, route: convergent, route: one-pot"
    And user makes row 1 current
    And the grid shows only the columns "id, reaction"
    And the "id" column of the grid is 200 px wide
    And the "reaction" column of the grid is 900 px wide
    And the rows of the grid are 150 px tall
    When the video says "Reactions and routes, each drawn to fit its cell"
    And user hovers over the "cell 32 of reaction" area of grid
    And user hovers over the "cell 35 of reaction" area of grid
    And user hovers over the "cell 38 of reaction" area of grid
    Given user opens the "../files/reactions-cdxml.csv" file of the project
    And the grid shows only the columns "name, cdxml"
    And the "name" column of the grid is 220 px wide
    And the "cdxml" column of the grid is 420 px wide
    And the rows of the grid are 150 px tall
    When the video says "ChemDraw pages, cut to their reaction"
    And user hovers over the "cell 3 of cdxml" area of grid
    And user hovers over the "cell 4 of cdxml" area of grid
    And user drags the "column resizer cdxml" area of grid by 560 pixels to the right
    And user hovers over the "cell 3 of cdxml" area of grid
    And user hovers over the "cell 4 of cdxml" area of grid
    And user hovers over the "cell 5 of cdxml" area of grid
    Then no errors should have been logged

  Scenario: The reaction editor, from a cell
    A double-click opens the Crux reaction editor, full screen or resized; OK writes the edit. In RDKit mode too.

    Given user is logged in
    And simple mode is off
    And the Chemical Engine is "Crux"
    And user opens the "fixtures/BDD-reactions-showcase.csv" file of the project
    And the toolbox pane is hidden
    And the context panel is closed
    Given user filters rows where "kind" is "reaction"
    And user makes row 20 current
    And the grid shows only the columns "name, reaction"
    And the "name" column of the grid is 190 px wide
    And the "reaction" column of the grid is 760 px wide
    And the rows of the grid are 130 px tall
    When the video says "Double-click a reaction: the Crux reaction editor"
    And user double-clicks on the "cell 1 of reaction" area of grid
    Then reaction editor dialog should be visible
    When user clicks on "Full screen" icon in reaction editor dialog
    And user clicks on "Exit full screen" icon in reaction editor dialog
    And the video says "The drawing follows the dialog"
    And user resizes reaction editor dialog by 260 and 0 pixels from its corner
    And the video says "A pyridine for the phenyl; OK writes it back"
    And user clicks on crux nitrogen tool
    And user clicks on the "atom 16" area of crux sketcher widget
    And user clicks on the "atom 43" area of crux sketcher widget
    And user clicks on OK button in reaction editor dialog
    Then reaction editor dialog should be absent
    And the reaction in row 1 of "reaction" column should have the reactants "CC(=O)c1ccc(Br)cc1.OB(O)c1ccncc1" and the products "CC(=O)c1ccc(-c2ccncc2)cc1"
    Given the Chemical Engine of this page is "RDKit"
    When the video says "RDKit mode: the same editor"
    And user double-clicks on the "cell 2 of reaction" area of grid
    Then reaction editor dialog should be visible
    When user clicks on CANCEL button in reaction editor dialog
    Then reaction editor dialog should be absent
    And no errors should have been logged

  Scenario: A ChemDraw page, through the editor and back
    A ChemDraw page keeps its labels through Crux's copy and paste; its Reaction pane says what lies on no step.

    Given user is logged in
    And simple mode is off
    And the Chemical Engine is "Crux"
    And user opens the "../files/reactions-cdxml.csv" file of the project
    And the toolbox pane is hidden
    And the context panel is open
    And the grid shows only the columns "name, cdxml"
    And the "name" column of the grid is 220 px wide
    And the "cdxml" column of the grid is 520 px wide
    And the rows of the grid are 150 px tall
    When the video says "Molecules on no step are counted, not dropped silently"
    And user clicks on the "cell 5 of cdxml" area of grid
    Then the context panel should show the current cell
    When user expands "Reaction" section
    Then the Reaction pane should say "2 molecules on no step were left out"
    When the video says "Copied and pasted, the page keeps OTf and TBS"
    And user double-clicks on the "cell 5 of cdxml" area of grid
    Then reaction editor dialog should be visible
    When user presses Control+A in crux canvas
    Then no rows should be selected
    When user presses Control+C in crux canvas
    And user clicks on crux clear button
    And user presses Control+V in crux canvas
    And user points at the label "OTf" in the sketcher
    And user points at the label "TBS" in the sketcher
    Then the sketcher should draw the labels "OTf" and "TBS"
    When user clicks on CANCEL button in reaction editor dialog
    Then no errors should have been logged
