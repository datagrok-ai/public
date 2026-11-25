@guide @demo
Feature: Chemistry · Reactions drawn by Crux
  Demo clips of Chem's reactions in Crux mode (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves the
  set-up and the checks off camera), over the showcase table fixtures/BDD-reactions-showcase.csv, opened from the
  project (nothing written to the stand): named reactions with their agents and the product beside each, four
  synthetic routes with a ring scaffold, and seven of Chem's enumeration templates. Crux draws a reaction with the
  molecule options Chem's molecule cells use, a route as one chain with numbered steps, and a template as a generic
  scheme whose notes the Reaction pane lists; the Reaction pane draws the current reaction large, its steps, and maps
  a route's atoms step by step. The last clip switches the Chemical Engine to RDKit and back in Chem's settings, for
  the account's own group (the server gives an account its group's settings over All users'); the settings are the
  stand's, so the step that remembers them puts them back at feature end. What Crux draws, atom by atom, is Chem's
  package tests ('Reactions: Crux').

  Scenario: Reactions drawn like the molecules beside them
    Crux draws each reaction in Chem's molecule style, a thin arrow with the agents over it; widen the column and the scheme grows.

    Given user is logged in
    And simple mode is off
    And the Chemical Engine is "Crux"
    And user opens the "fixtures/BDD-reactions-showcase.csv" file of the project
    And the toolbox pane is hidden
    And the context panel is closed
    Then "reaction" column should have semantic type "ChemicalReaction"
    And "product" column should have semantic type "Molecule"
    # a row the filter hides is current: no row of the clip is highlighted
    Given user filters rows where "kind" is "reaction"
    And user makes row 20 current
    And the grid shows only the columns "name, reaction, product"
    And the "name" column of the grid is 190 px wide
    And the "reaction" column of the grid is 300 px wide
    And the "product" column of the grid is 220 px wide
    And the rows of the grid are 130 px tall
    When the video says "Reactions in the style of the molecule cells"
    And user hovers over the "cell 1 of reaction" area of grid
    And user hovers over the "cell 1 of product" area of grid
    And user drags the "column resizer reaction" area of grid by 640 pixels to the right
    And user hovers over the "cell 2 of reaction" area of grid
    And user hovers over the "cell 2 of product" area of grid
    And user scrolls the mouse wheel down 3 times over the "cell 4 of reaction" area of grid
    And user hovers over the "cell 7 of reaction" area of grid
    And user hovers over the "cell 8 of reaction" area of grid
    And user hovers over the "cell 8 of product" area of grid
    And user hovers over the "cell 9 of reaction" area of grid
    And user hovers over the "cell 9 of product" area of grid
    Then no errors should have been logged

  Scenario: Routes drawn as whole chains
    A multi-step route is one chain, its steps numbered and its scaffold kept in place from step to step.

    Given user is logged in
    And simple mode is off
    And the Chemical Engine is "Crux"
    And user opens the "fixtures/BDD-reactions-showcase.csv" file of the project
    And the toolbox pane is hidden
    And the context panel is closed
    Given user filters rows where "kind" is "route"
    And user makes row 20 current
    And the grid shows only the columns "name, reaction"
    And the "name" column of the grid is 190 px wide
    And the "reaction" column of the grid is 760 px wide
    And the rows of the grid are 180 px tall
    When the video says "One chain: the steps numbered, the ring kept in place"
    And user hovers over the "cell 10 of reaction" area of grid
    And user hovers over the "cell 12 of reaction" area of grid
    And user drags the "column resizer reaction" area of grid by 540 pixels to the right
    And user hovers over the "cell 10 of reaction" area of grid
    And user hovers over the "cell 11 of reaction" area of grid
    And user hovers over the "cell 12 of reaction" area of grid
    And user hovers over the "cell 13 of reaction" area of grid
    Then no errors should have been logged

  Scenario: The Reaction pane maps a route step by step
    Click a route: the context panel draws it large with its steps and maps each step by itself (it has no button), and full screen shows them mapped.

    Given user is logged in
    And simple mode is off
    And the Chemical Engine is "Crux"
    And user opens the "fixtures/BDD-reactions-showcase.csv" file of the project
    And the toolbox pane is hidden
    Given user filters rows where "kind" is "route"
    And user makes row 20 current
    And the context panel is open
    And user drags context panel splitter by -300 and 0 pixels
    And the grid shows only the columns "name, reaction"
    And the "name" column of the grid is 170 px wide
    And the "reaction" column of the grid is 560 px wide
    And the rows of the grid are 170 px tall
    When user clicks on the "cell 10 of reaction" area of grid
    Then the context panel should show the current cell
    When user expands "Reaction" section
    And the video says "Atoms mapped by themselves, one step at a time"
    Then the Reaction pane should have mapped 3 steps
    When user clicks on "Full screen" icon inside "Reaction" section
    Then "Reaction" dialog should be visible
    When user scrolls the Reaction dialog to its last step
    And user closes "Reaction" dialog
    Then no errors should have been logged

  Scenario: Templates drawn as generic schemes
    Chem's reaction templates are drawn with R, Ar and X; click one and the context panel lists what its atoms stand for.

    Given user is logged in
    And simple mode is off
    And the Chemical Engine is "Crux"
    And user opens the "fixtures/BDD-reactions-showcase.csv" file of the project
    And the toolbox pane is hidden
    Given user filters rows where "kind" is "template"
    And user makes row 1 current
    And the context panel is open
    And user drags context panel splitter by -250 and 0 pixels
    And the grid shows only the columns "name, reaction"
    And the "name" column of the grid is 200 px wide
    And the "reaction" column of the grid is 560 px wide
    And the rows of the grid are 110 px tall
    When the video says "Templates as generic schemes: R, Ar, X"
    And user hovers over the "cell 14 of reaction" area of grid
    And user hovers over the "cell 16 of reaction" area of grid
    And user clicks on the "cell 15 of reaction" area of grid
    Then the context panel should show the current cell
    When user expands "Reaction" section
    And user scrolls reaction template notes into view
    Then reaction template notes should be visible
    When user clicks on the "cell 17 of reaction" area of grid
    And user scrolls reaction template notes into view
    Then reaction template notes should be visible
    When user clicks on the "cell 16 of reaction" area of grid
    And user scrolls reaction template notes into view
    Then reaction template notes should be visible
    And no errors should have been logged

  Scenario: Switching the Chemical Engine
    In RDKit mode Chem draws reactions as it always has; switched back to Crux, the same cells are drawn by Crux.

    Given user is logged in
    And simple mode is off
    And the Chemical Engine is "Crux"
    And the "Chem" package settings are put back at feature end
    And user opens the "fixtures/BDD-reactions-showcase.csv" file of the project
    Given user filters rows where "name" is one of "Suzuki coupling, Amide coupling with HATU, Paracetamol from phenol"
    And user makes row 20 current
    And the grid shows only the columns "name, reaction"
    And the "name" column of the grid is 170 px wide
    And the "reaction" column of the grid is 620 px wide
    And the rows of the grid are 150 px tall
    And the context panel is open
    And the browse panel is open
    And Platform tree node inside browse tree is expanded
    When user clicks on Platform---Plugins tree node inside browse tree
    And user clicks on "Chem" gallery card
    And user expands Settings accordion header in context panel
    And user clicks on package settings group button in "Settings" pane in context panel
    And user picks their own group from the package settings groups
    And user selects "RDKit" in "Chemical Engine" input in "Settings" pane in context panel
    And user clicks on SAVE button in "Settings" pane in context panel
    And user clicks on BDD-reactions-showcase tab
    Then the Chemical Engine should be "RDKit"
    When the video says "RDKit mode: drawn as Chem always has"
    And user hovers over the "cell 1 of reaction" area of grid
    And user hovers over the "cell 10 of reaction" area of grid
    And user clicks on Plugins tab
    And user clicks on "Chem" gallery card
    And user expands Settings accordion header in context panel
    And user clicks on package settings group button in "Settings" pane in context panel
    And user picks their own group from the package settings groups
    And user selects "Crux" in "Chemical Engine" input in "Settings" pane in context panel
    And user clicks on SAVE button in "Settings" pane in context panel
    And user clicks on BDD-reactions-showcase tab
    Then the Chemical Engine should be "Crux"
    When the video says "Crux mode: the same cells, drawn by Crux"
    And user hovers over the "cell 1 of reaction" area of grid
    And user hovers over the "cell 10 of reaction" area of grid
    Then no errors should have been logged
