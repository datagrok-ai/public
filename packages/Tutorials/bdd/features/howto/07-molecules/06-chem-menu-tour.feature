@guide @tutorial @viewport:1280x720
Feature: How-to · Molecules
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). A walk down the Chem menu, which a table view shows once
  a column holds molecules, with one line per group: Search, Calculate, Analyze, Transform, then Admetica. The
  pointer rests on a few items of each group; nothing is picked, and nothing that runs Python, a container or an
  outside service is opened (Descriptors, Map Identifiers, IUPAC Name, Curate, Mutate, Generate Conformers and
  Admetica stay closed). The pointer goes the way a hand does (`user rests the pointer on {string} in the top
  menu`, bindings/howto-data-2.ts): the Dart menu hides a flyout left at a steep angle. Docking's package is on
  the stand but adds nothing to this menu. smiles is opened off camera; nothing is made on the server.

  Scenario: What the Chem menu does
    The Chem menu holds the tools for a table of molecules: search, calculate, analyze and transform.

    Given user is logged in
    And simple mode is off
    And the "Chem" package is installed
    And the package autostarts have completed
    And user opens smiles dataset
    And the grid shows only the columns "molregno, canonical_smiles, RingCount, NumHDonors, NumHAcceptors"
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    When the video says "Chem: the tools for a table of molecules" at the top
    # caption: Search: by a fragment, by similarity, or for diversity
    And user rests the pointer on "Chem > Search" in the top menu
    And user rests the pointer on "Chem > Search > Similarity Search..." in the top menu
    And user rests the pointer on "Chem > Search > Substructure Search..." in the top menu
    Then "Substructure Search..." menu item should be enabled
    # caption: Calculate: properties, scores and alerts as new columns
    When user rests the pointer on "Chem > Calculate" in the top menu
    And user rests the pointer on "Chem > Calculate > Toxicity Risks..." in the top menu
    And user rests the pointer on "Chem > Calculate > MPO Score..." in the top menu
    And user rests the pointer on "Chem > Calculate > Chemical Properties..." in the top menu
    Then "Chemical Properties..." menu item should be enabled
    # caption: Analyze: how structure relates to activity
    When user rests the pointer on "Chem > Analyze" in the top menu
    And user rests the pointer on "Chem > Analyze > Chemical Space..." in the top menu
    And user rests the pointer on "Chem > Analyze > R-Groups Analysis..." in the top menu
    And user rests the pointer on "Chem > Analyze > Activity Cliffs..." in the top menu
    And user rests the pointer on "Chem > Analyze > Scaffold Tree" in the top menu
    Then "Activity Cliffs..." menu item should be enabled
    # caption: Transform: change the molecules: notation, cleanup, reactions
    When user rests the pointer on "Chem > Transform" in the top menu
    And user rests the pointer on "Chem > Transform > Convert Notation..." in the top menu
    And user rests the pointer on "Chem > Transform > Reactions" in the top menu
    Then "Convert Notation..." menu item should be enabled
    # caption: Admetica: ADMET predictions, computed on the server
    When user rests the pointer on "Chem > Admetica" in the top menu
    Then "Admetica" menu item should be enabled
    And no errors should have been logged
