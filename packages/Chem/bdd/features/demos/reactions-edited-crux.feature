@guide @demo
Feature: Chemistry · Reactions edited in Crux Sketch
  Demo clips of the reaction editor (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves the set-up and
  the checks off camera). The Reaction Sketcher edits its reaction in Crux Sketch, the text under it and the publication
  picture beside it following each other. The clips of a reaction cell edited in Crux Sketch and of a route the editor
  did not open are retired: reaction-tables.feature's "The reaction editor, from a cell" opens, edits and writes back a
  cell in both engine modes, and the editor opens convergent routes too since crux-core 15b1bbd65 (Chem's package tests
  keep the refusal of a route it cannot read). What the editor writes for the corpus and Chem's templates is Chem's
  package tests ('Reactions: editor', 'Reactions: app').

  Scenario: The Reactions app draws what you type, and writes what you draw
    Type a reaction and Crux Sketch draws it; draw on it, and the text and the publication picture follow.

    Given user is logged in
    And simple mode is off
    And the Chemical Engine is "Crux"
    And user opens Chem's Reactions app
    And the toolbox pane is hidden
    And the context panel is closed
    When the video says "Type a reaction: Crux Sketch draws it"
    And user enters "Benzylamine acetylated" into "Title" input
    And user types "CC(=O)Cl.NCc1ccccc1>>CC(=O)NCc1ccccc1" into reactions app text
    Then the "atoms" reading of crux sketcher widget should be 23
    When the video says "Draw on it: the text and the picture follow"
    And user clicks on crux single bond tool
    # a methyl on the amine's nitrogen, then on the amide's (the atoms in the order the text writes them)
    And user clicks on the "atom 4" area of crux sketcher widget
    And user clicks on the "atom 15" area of crux sketcher widget
    Then the Reactions app's text should have the reactants "CC(=O)Cl.CNCc1ccccc1" and the products "CC(=O)N(C)Cc1ccccc1"
    When user hovers over reactions app drawing
    Then no errors should have been logged
